package phylostream;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.UndirectedGraph;

import bayonet.marginal.DiscreteFactorGraph;
import bayonet.marginal.FactorGraph;
import bayonet.marginal.UnaryFactor;
import briefj.collections.UnorderedPair;
import conifer.EvolutionaryModel;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import conifer.models.LikelihoodFactoryContext;

/**
 * A rooted phylogenetic tree supporting 'growth' i.e. adding tips, as well as basic 
 * tree operators such as changing a branch length or performing NNI operators.
 * 
 * See TestGrowingTree for usage.
 * . 
 * Keeps track of dynamic programming messages to ensure quick likelihood evaluation,
 * via IncrementalSumProduct (one per rate category). 
 */
public class GrowingTree {
  
  // TODO: to support quick resampling, will need to switch to VAVR's Hash array mapped trie (https://en.wikipedia.org/wiki/Hash_array_mapped_trie)
  //       and check everything is working fine with cloning (or write custom cloning)
  //       annoying part with this will be to migrate the jgrapht SimpleGraph, UnrootedTree, and bayonet FactorGraphs to use VAVR 
  //       best way to do this is probably write a java.util.Collection with an inner VAVR object.. not sure-issue is then changes via iterators, etc
  
  // NB: both 'tree' and 'likelihood.graph' share the same underlying topology object
  final UnrootedTree unrootedTree; 
  final List<IncrementalSumProduct<TreeNode>> likelihoods; // one for each rate category
  final TreeNode root;
  final Map<TreeNode, TreeNode> rootingParentPointers; // node -> parent 
  
  final EvolutionaryModel model;
  final TreeObservations observations;
  
  private double annealingParameter;
  
  public GrowingTree(UnrootedTree tree, TreeNode root, EvolutionaryModel model, TreeObservations observations) {
    this.root = root;
    this.model = model;
    this.observations = observations;
    this.unrootedTree = tree;
    TreeNode latestTip = tree.leaves().stream().filter((TreeNode it) -> it != root).findFirst().get();
    this.annealingParameter = 1.0;
    List<FactorGraph<TreeNode>> graphs = EvolutionaryModelUtils.buildFactorGraphs(model, tree, root, observations);
    this.likelihoods = new ArrayList<IncrementalSumProduct<TreeNode>>(graphs.size());
    for (int category = 0; category < graphs.size(); category++) 
      likelihoods.add(new IncrementalSumProduct<TreeNode>(graphs.get(category), latestTip));
    this.rootingParentPointers = parentsMap(tree, root);
  }
  
  public double logLikelihood() {
    List<UnaryFactor<TreeNode>> marginals = new ArrayList<UnaryFactor<TreeNode>>(nCategories());
    for (IncrementalSumProduct<TreeNode> sumProduct : likelihoods)
      marginals.add(sumProduct.tipMarginals());
    LikelihoodComputationContext context = new LikelihoodComputationContext(marginals);
    return model.computeLogLikelihood(context);
  }

  /**
   * Setup the factor graph and tree object to accommodate one new tip. 
   * 
   * Initially, the tip is set to have an anneal parameter of 0.0. 
   * 
   * Notation used for the input argument setting up initial branch lengths l0 and l1
   * and placement point:
   * 
   *         o freshLatestTip
   *         |
   *         | l0
   *  \      |      /
   * --o-----o-----o--
   *   v l1  x     w
   */
  public void addTip(TreeNode freshLatestTip, TreeNode v, TreeNode w, double l0, double l1) {  
    if (unrootedTree.getTopology().containsVertex(freshLatestTip))
      throw new RuntimeException();
    double oldLength = getBranchLength(v, w);
    TreeNode x = TreeNode.nextUnlabelled();
    if (l1 > oldLength) throw new RuntimeException();
    removeBranch(v, w);
    topology().addVertex(freshLatestTip);
    topology().addVertex(x);
    updateRooting(v, w, x, freshLatestTip);
    addBranch(v, x, l1);
    addBranch(x, w, oldLength - l1);
    addBranch(x, freshLatestTip, l0);
    // NB: NOT adding the tip unary since it is completely annealed.
    this.annealingParameter = 0.0;
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.updateTip(freshLatestTip);
  }
  
  public double getBranchLength(TreeNode v, TreeNode w) {
    return unrootedTree.getBranchLength(v, w);
  }
  
  /**
   * Do not edit the result. Use updateBranchLength() instead.
   */
  public Map<UnorderedPair<TreeNode, TreeNode>, Double> getBranchLengths() {
    return unrootedTree.getBranchLengths(); 
  }
  
  public void updateBranchLength(TreeNode x, TreeNode y, double length) {
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.notifyFactorUpdated(UnorderedPair.of(x, y)); 
    removeBranch(x, y);
    addBranch(x, y, length); 
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.recomputeMessages();
  }
  
  /**
   * Perform an interchange between the latestTip and the edge 
   * (rootOfDisconnectedSubtree, otherEndPointOfCutEdge). 
   * 
   * The move creates the edges (latestTip, otherEndPointOfCutEdge) and
   * (rootOfDisconnectedSubtree, x), 
   * where x = unique neighbour of latestTip before the move
   * 
   * Schematic of the notation before the move:
   * 
   *         o latestTip                v [possibly a subtree]  
   *         |                          o rootOfDisconnectedSubtree
   *         | latestBranchLen          | otherBranchLen
   *  \      |                          |     
   * --o-----o--------------------------o-- ...
   *         x                          otherEndPointOfCutEdge
   *         
   * Schematic of the notation after the move:
   * 
   *         v [possibly a subtree]      o latestTip       
   *         o rootOfDisconnectedSubtree |    
   *         | otherBranchLen            | latestBranchLen
   *  \      |                           |     
   * --o-----o---------------------------o-- ...
   *         x                           otherEndPointOfCutEdge
   * 
   * Return the arguments used to undo this move. 
   */
  public Pair<TreeNode, TreeNode> interchange(TreeNode rootOfDisconnectedSubtree, TreeNode otherEndPointOfCutEdge) {
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.notifyFactorUpdated(UnorderedPair.of(rootOfDisconnectedSubtree, otherEndPointOfCutEdge)); 
    TreeNode x = latestEdge().getLeft();
    if (!topology().containsEdge(x, otherEndPointOfCutEdge))
      throw new RuntimeException();
    double latestBranchLen = unrootedTree.getBranchLength(latestTip(), x);
    double otherBranchLen = unrootedTree.getBranchLength(rootOfDisconnectedSubtree, otherEndPointOfCutEdge);
    removeBranch(x, latestTip());
    removeBranch(rootOfDisconnectedSubtree, otherEndPointOfCutEdge);
    
    // Update the oriention before adding the new branches
    if (rootingParentPointers.get(rootOfDisconnectedSubtree) == otherEndPointOfCutEdge) {
      /*
       * This means the config before was 
       * 
       *         o latestTip                v  
       *         ^                          o rootOfDisconnectedSubtree
       *         |                          ^ 
       *  \      |                          |     
       * --o-----o?------------------------?o-- ...
       *         x                          otherEndPointOfCutEdge
       *         
       * and hence after the move we should have
       * 
       *         v [possibly a subtree]      o latestTip       
       *         o rootOfDisconnectedSubtree ^    
       *         ^                           | 
       *  \      |                           |     
       * --o-----o?-------------------------?o-- ...
       *         x                           otherEndPointOfCutEdge
       *         
       * where '?' denotes an orientation that do not matter
       * since it will be the same before and after.
       */
      rootingParentPointers.put(rootOfDisconnectedSubtree, x);
    } else { 
      /*
       * This means the config before was 
       * 
       *         o latestTip                v [possibly a subtree]  
       *         ^                          o rootOfDisconnectedSubtree
       *         |                          | 
       *  \      |                          v     
       * --o<----o<-------------------------o--> ...
       *         x                          otherEndPointOfCutEdge
       *         
       * and hence after the move we should have
       * 
       *         v                           o latestTip       
       *         o rootOfDisconnectedSubtree ^    
       *         |                           | 
       *  \      v                           |     
       * --o<----o-------------------------->o--> ...
       *         x                           otherEndPointOfCutEdge
       *         
       * I.e. in this sub-case the root was within the displaced 
       * subtree, hence the edge between x and otherEndPointOfCutEdge 
       * needs to be inverted.      
       */
      rootingParentPointers.put(x, rootOfDisconnectedSubtree);
      rootingParentPointers.put(otherEndPointOfCutEdge, x);
    }
    rootingParentPointers.put(latestTip(), otherEndPointOfCutEdge);
    
    addBranch(otherEndPointOfCutEdge, latestTip(), latestBranchLen);
    addBranch(rootOfDisconnectedSubtree, x, otherBranchLen);
    
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.recomputeMessages();
    
    return Pair.of(rootOfDisconnectedSubtree, x);
  }
  
  public double getLatestTipAnnealingParameter() { return annealingParameter; }
  
  public void setLatestTipAnnealingParameter(double annealingParameter) {
    if (this.annealingParameter == annealingParameter) return;
    this.annealingParameter = annealingParameter;
    for (int category = 0; category < likelihoods.size(); category++) {
      DiscreteFactorGraph<TreeNode> graph = graph(category);
      if (graph.getUnary(latestTip()) != null)
        graph.removeUnary(latestTip());
      if (annealingParameter != 0.0)
        model.buildObservation(latestTip(), context(category, annealingParameter));
    }
  }
  
  void updateRooting(TreeNode _v, TreeNode _w, TreeNode x, TreeNode freshLatestTip) {
    Pair<TreeNode, TreeNode> directed = orient(_v, _w);
    TreeNode v = directed.getLeft();
    TreeNode w = directed.getRight();
    /*
     * After calling orient, we can assume the following directed structure:
     * 
     *         o freshLatestTip
     *         ^
     *         | l0
     *  \      |      /
     * --o---->o---->o--
     *   v l1  x     w
     */
    rootingParentPointers.put(w, x);
    rootingParentPointers.put(x, v);
    rootingParentPointers.put(freshLatestTip, x);
  }
  
  void addBranch(TreeNode x, TreeNode y, double length) {
    unrootedTree.addEdge(x, y, length);
    Pair<TreeNode, TreeNode> directedEdge = orient(x, y);
    for (int category = 0; category < nCategories(); category++) 
      model.buildTransition(directedEdge.getLeft(), directedEdge.getRight(), context(category)); 
  }
  
  LikelihoodFactoryContext context(int category) { return context(category, Double.NaN); } // when the context will not make use of annealing parameter
  LikelihoodFactoryContext context(int category, double annealingParameter) {
    return new LikelihoodFactoryContext(likelihoods.get(category).graph, unrootedTree, observations, category, annealingParameter);
  }
  
  UndirectedGraph<TreeNode, ?> topology() {
    return unrootedTree.getTopology();
  }
  
  DiscreteFactorGraph<TreeNode> graph(int category) {
    return (DiscreteFactorGraph<TreeNode>) (likelihoods.get(category).graph);
  }
  
  TreeNode latestTip() { return likelihoods.get(0).latestTip; }
  Pair<TreeNode, TreeNode> latestEdge() { return likelihoods.get(0).latestEdge(); } 
  
  void removeBranch(TreeNode v, TreeNode w) {
    unrootedTree.removeEdge(v, w);
    for (int category = 0; category < nCategories(); category++) {
      DiscreteFactorGraph<TreeNode> graph = graph(category);
      graph.removeBinary(v, w);
      graph.removeBinary(w, v);
    }
  }

  int nCategories() { return likelihoods.size(); }
  
  /**
   * From an undirected edge {x, y}, return a directed edge 
   * (result.getLeft(), result.getRight()) pointing away from 
   * the root.
   */
  Pair<TreeNode, TreeNode> orient(TreeNode x, TreeNode y) {
         if (rootingParentPointers.get(x) == y) return Pair.of(y, x);
    else if (rootingParentPointers.get(y) == x) return Pair.of(x, y);
    else throw new RuntimeException();
  }
  
  public static Map<TreeNode, TreeNode> parentsMap(UnrootedTree tree,
      TreeNode root) {
    Map<TreeNode, TreeNode> result = new LinkedHashMap<TreeNode, TreeNode>();
    for (Pair<TreeNode,TreeNode> directed : tree.getRootedEdges(root))
      result.put(directed.getRight(), directed.getLeft());
    return result;
  }
  
  @Override
  public String toString() {
    return unrootedTree.toString();
  }
}
