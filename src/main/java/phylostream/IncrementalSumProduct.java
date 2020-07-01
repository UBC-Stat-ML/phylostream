package phylostream;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import bayonet.graphs.GraphUtils;
import bayonet.marginal.BinaryFactor;
import bayonet.marginal.DiscreteFactorGraph;
import bayonet.marginal.FactorGraph;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;
import briefj.BriefCollections;
import briefj.collections.UnorderedPair;

/**
 * Manages a cache suitable for incremental computation of logNormalization() using the sum product algorithm. 
 * 
 * Keeps track of a distinguished leaf, called latestTip, which orients all message passing computations 
 * towards it. 
 * 
 * This object should be notified of all updates on the topology and factors. See notifyFactorUpdated().
 */
public class IncrementalSumProduct<N> {
  
  final FactorGraph<N> graph;
  
  N latestTip;
  
  final Map<Pair<N, N>, UnaryFactor<N>> cachedMessages; 
  
  public IncrementalSumProduct(DiscreteFactorGraph<N> graph, N latestTip) {
    // TODO: check graph.getTopology is a tree
    this.graph = graph;
    this.latestTip = latestTip;
    if (topology().degreeOf(latestTip) != 1) 
      throw new RuntimeException("" + latestTip + " should be a leaf.");
    // This implementation is based on memoization because updates are assumed small.
    // However, for initialization, use the sum product machinery (based on loops) to avoid potentially 
    // costly recursion for the initialization phase.
    SumProduct<N> sumProduct = new SumProduct<>(graph);
    sumProduct.computeMarginal(latestTip); 
    cachedMessages = sumProduct.cachedMessages;
  }
  
  public double logNormalization() {
    UnaryFactor<N> result = message(latestEdge()); //cachedMessages.get(latestEdge());
    UnaryFactor<N> modelFactor = graph.getUnary(latestTip);
    if (modelFactor != null) {
      List<UnaryFactor<N>> toMultiply = new ArrayList<UnaryFactor<N>>(2);
      toMultiply.add(modelFactor);
      toMultiply.add(result);
      result = graph.factorOperations().pointwiseProduct(toMultiply);
    }
    return result.logNormalization();
  }
  
  /**
   * Signal that a new tip called 'freshLatestTip' is now the reference tip.
   * 
   * Call this AFTER performing the topological change.
   * 
   * @param freshLatestTip
   */
  public void updateTip(N freshLatestTip) {  
    latestTip = freshLatestTip;
    // trigger message computation; this will re-orient messages towards new tip, see (*) in message(Pair)
    messageToLatestTip();
  }
  
  /**
   * Delete the cache corresponding to an edge in the graph, along with the 
   * other messages that depended on it. 
   * 
   * When performing NNI, call this method on the edge that latestTip will be exchanged with; 
   * doing so BEFORE performing the topological change, to allow proper tracking of message dependencies.
   * 
   * @param edge
   */
  public void notifyFactorUpdated(UnorderedPair<N, N> edge) {
    notifyFactorUpdated(orient(edge));
  }
  
  // Support methods below
  
  UnaryFactor<N> messageToLatestTip() {
    UnaryFactor<N> result = message(latestEdge());
    if (cachedMessages.size() != topology().edgeSet().size()) 
      throw new RuntimeException("Potential memory leak detected. Cache after recomputation should be of the same size as the number of edges."); 
    return result;
  }
  
  UnaryFactor<N> message(Pair<N,N> edge) {
    if (cachedMessages.containsKey(edge)) return cachedMessages.get(edge);
    cachedMessages.remove(reverse(edge)); // rationale: see (*) in updateTip()
    N source = edge.getLeft(),
      destination = edge.getRight();

    List<UnaryFactor<N>> toMultiply = new ArrayList<>();
    for (Pair<N,N> incomingEdge : GraphUtils.distinctIncoming(topology(), edge))
      toMultiply.add(message(incomingEdge));
    UnaryFactor<N> modelFactor = graph.getUnary(source);
    if (modelFactor != null)
      toMultiply.add(modelFactor);
    
    // marginalize one node
    BinaryFactor<N> binaryFactor = graph.getBinary(source, destination);
    UnaryFactor<N> result = graph.factorOperations().marginalize(binaryFactor, toMultiply);
    cachedMessages.put(edge, result);
    return result;
  }
  
  void notifyFactorUpdated(Pair<N, N> edge) {
    // orient using the cache
    if (edge == null) return; 
    cachedMessages.remove(edge);
    for (N neighbour : Graphs.neighborListOf(topology(), edge.getRight()))
      if (neighbour != edge.getLeft())
        notifyFactorUpdated(Pair.of(edge.getRight(), neighbour));
  }
  
  /** Use the cached message to orient the edge towards the latestTip */
  Pair<N,N> orient(UnorderedPair<N, N> edge) {
    Pair<N,N> o1 = Pair.of(edge.getFirst(), edge.getSecond()); if (!cachedMessages.containsKey(o1)) o1 = null;
    Pair<N,N> o2 = Pair.of(edge.getSecond(), edge.getFirst()); if (!cachedMessages.containsKey(o2)) o2 = null;
    if (o1 == null && o2 == null) return null;
    if (o1 != null && o2 != null) throw new RuntimeException();
    return o1 != null ? o1 : o2;
  }
  
  static <N> Pair<N,N> reverse(Pair<N,N> edge) { return Pair.of(edge.getRight(), edge.getLeft()); }
  
  /** The edge pointing towards latestTip */
  Pair<N,N> latestEdge() {
    List<N> neighbours = Graphs.neighborListOf(topology(), latestTip);
    if (neighbours.size() != 1) throw new RuntimeException("" + latestTip + " should be a leaf.");
    return Pair.of(BriefCollections.pick(neighbours), latestTip);
  }
  
  UndirectedGraph<N, ?> topology() {
    return graph.getTopology();
  }
}
