package phylostream;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import org.jgrapht.UndirectedGraph;

import bayonet.marginal.BinaryFactor;
import bayonet.marginal.DiscreteFactorGraph;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;
import briefj.BriefCollections;

/**
 * Incremental computation of sum product. 
 * 
 * Keeps track of the latestTip, and allow local changes around that tip without full re-computation, 
 * via incremental updates of the messages.
 * 
 * (*) Invariant: keep the messages computed pointing towards latestTip.
 */
public class IncrementalSumProduct<N> {
  
  final DiscreteFactorGraph<N> graph;
  
  N latestTip;
  
  final Map<Pair<N, N>, UnaryFactor<N>> cachedMessages; 
  
  public double logNormalization() {
    // By (*) this message should already be available:
    UnaryFactor<N> result = cachedMessages.get(latestEdge());
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
   * Engraft new tip called 'freshLatestTip' on the edge 
   * engraftmentEdge, with end points (v, w) shown below in the ASCII figure. 
   * 
   * Use an array newFactors on the 3 introduced edges, 
   * shown as {f0, f1, f2} in the ASCII figure.
   * 
   *         o freshLatestTip
   *         |
   *         | f0
   *  \      |      /
   * --o-----o-----o--
   *   v f1    f2  w
   * 
   * @param freshLatestTip
   * @param engraftmentEdge
   * @param newFactors
   */
  public void addTip(N freshLatestTip, Pair<N,N> engraftmentEdge, BinaryFactor<N> [] newFactors) {  
    throw new RuntimeException("not implemented");
    // TODO: check freshLatestTip is not in there
    // TODO: update topology
    // TODO: add binary factor
    // TODO: update cache
    // TODO: update latestTip
  }
  

  
//  /**
//   * 
//   * Involution property: to faciliate MCMC rejection, we have that the following code should bring back the configuration to the same state as the start state
//   * 
//   * result = localUpdate(x, y)
//   * localUpdate(x, result)
//   * 
//   * @param topologyCode 0: no topology change; 1,2,3,4: the 4 possible NNI moves; should be an involution, i.e. such that applying the same code twice go back to original state
//   * @param updatedFactors length 7; contains the 7 neighbour factors that can be updated
//   * @return an array of the previous binary factor that can be passed in to get the involution
//   */
//  public BinaryFactor<N> [] localUpdate(int topologyCode, BinaryFactor<N> [] updatedFactors) {
//    throw new RuntimeException("not implemented");
//  }
  
  /**
   * 
   *         x   freshLatestTip
   *         |
   *         |   f0
   *  \      |      /
   * --x-----x-----x--
   *   left        right
   *   
   *     f1     f2
   * 
   * @param updates
   * @return
   */
  public BinaryFactor<N> updateLocalBranch(N neighbour, BinaryFactor<N> update) {
    throw new RuntimeException("not implemented");
  }
  
  public void interchange(Pair<N,N> otherEdge) {
    throw new RuntimeException("not implemented");
  }
  
  public List<Pair<N,N>> interchangeNeighbours() {
    throw new RuntimeException("not implemented");
  }
  
  /** The edge pointing towards latestTip */
  Pair<N,N> latestEdge() {
    List<N> neighbours = Graphs.neighborListOf(topology(), latestTip);
    if (neighbours.size() != 1) throw new RuntimeException("" + latestTip + " should be a leaf.");
    return Pair.of(BriefCollections.pick(neighbours), latestTip);
  }
  
  IncrementalSumProduct(DiscreteFactorGraph<N> graph, N latestTip) {
    // TODO: check graph.getTopology is a tree
    this.graph = graph;
    this.latestTip = latestTip;
    if (topology().degreeOf(latestTip) != 1) 
      throw new RuntimeException("" + latestTip + " should be a leaf.");
    SumProduct<N> sumProduct = new SumProduct<>(graph);
    sumProduct.computeMarginal(latestTip); 
    cachedMessages = sumProduct.cachedMessages;
  }
  
  UndirectedGraph<N, ?> topology() {
    return graph.getTopology();
  }
}
