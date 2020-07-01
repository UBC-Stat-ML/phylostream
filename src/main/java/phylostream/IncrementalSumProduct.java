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
  
  boolean messagesComputed;
  
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
    messagesComputed = true;
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
    messagesComputed = true;
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
    if (!messagesComputed) throw new RuntimeException();
    UnorderedPair<N, N> deletedEdge = deletedEdge(latestEdge(freshLatestTip));
    notifyFactorUpdated(deletedEdge); 
    latestTip = freshLatestTip; 
    messageToLatestTip();
  }
  
  /**
   * Delete the cache corresponding to an edge in the graph, along with the 
   * other messages that depended on it. 
   * 
   * When performing NNI, call this method on the edge that latestTip will be exchanged with; 
   * doing so BEFORE performing the topological change, to allow proper tracking of message dependencies.
   * 
   * Several such notifications can be done in a row, however we assume 
   * logNormalization will be called between the last notification and the next 
   * call to updateTip().
   * 
   * @param edge
   */
  public void notifyFactorUpdated(UnorderedPair<N, N> edge) {
    messagesComputed = false;
    notifyFactorUpdated(orient(edge));
  }
  
  // Support methods below
  
  /**
   * Construct the pair {v, w} which used to be an edge an is no
   * longer one after the tip update.
   * 
   *         o 
   *         ^
   *         | newLatestEdge
   *  \      |      /
   * --o-----o-----o--
   *   v           w
   */
  UnorderedPair<N, N> deletedEdge(Pair<N,N> newLatestEdge) {
    List<N> nodes = new ArrayList<>(2);
    for (N neighbour : Graphs.neighborListOf(topology(), newLatestEdge.getLeft()))
      if (!neighbour.equals(newLatestEdge.getRight()))
        nodes.add(neighbour);
    return UnorderedPair.of(nodes.get(0), nodes.get(1));
  }
  
  UnaryFactor<N> messageToLatestTip() {
    UnaryFactor<N> result = message(latestEdge());
    if (cachedMessages.size() != topology().edgeSet().size()) 
      throw new RuntimeException("Potential memory leak detected. Cache after recomputation should be of the same size as the number of edges in the factor graph."); 
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
    if (edge == null) return; 
    cachedMessages.remove(edge);
    for (N neighbour : Graphs.neighborListOf(topology(), edge.getRight()))
      if (!neighbour.equals(edge.getLeft()))
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
  Pair<N,N> latestEdge() { return latestEdge(latestTip); }
  Pair<N,N> latestEdge(N latestTip) {
    List<N> neighbours = Graphs.neighborListOf(topology(), latestTip);
    if (neighbours.size() != 1) throw new RuntimeException("" + latestTip + " should be a leaf.");
    return Pair.of(BriefCollections.pick(neighbours), latestTip);
  }
  
  UndirectedGraph<N, ?> topology() {
    return graph.getTopology();
  }
}
