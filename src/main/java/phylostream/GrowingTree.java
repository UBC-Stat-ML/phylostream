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

public class GrowingTree {
  
  // TODO: to support quick resampling, will need to switch to VAVR's Hash array mapped trie (https://en.wikipedia.org/wiki/Hash_array_mapped_trie)
  //       and check everything is working fine with cloning (or write custom cloning)
  //       annoying part with this will be to migrate the jgrapht SimpleGraph, UnrootedTree, and bayonet FactorGraphs to use VAVR 
  //       best way to do this is probably write a java.util.Collection with an inner VAVR object.. not sure-issue is then changes via iterators, etc
  
  // NB: both 'tree' and 'likelihood.graph' share the same underlying topology object
  final UnrootedTree tree;
  final List<IncrementalSumProduct<TreeNode>> likelihoods; // one for each rate category
  final Map<TreeNode, TreeNode> rootingParentPointers; // node -> parent 
  
  final EvolutionaryModel model;
  final TreeObservations observations;
  
  public GrowingTree(UnrootedTree tree, TreeNode root, EvolutionaryModel model, TreeObservations observations) {
    this.model = model;
    this.observations = observations;
    this.tree = tree;
    TreeNode latestTip = tree.leaves().get(0);
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
    double oldLength = tree.getBranchLength(v, w);
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
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.updateTip(freshLatestTip);
  }
  
  public void updateBranchLength(TreeNode x, TreeNode y, double length) {
    removeBranch(x, y);
    addBranch(x, y, length); 
    for (IncrementalSumProduct<TreeNode> likelihood : likelihoods)
      likelihood.notifyFactorUpdated(UnorderedPair.of(x, y)); 
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
  
  public void setLatestTipAnnealingParameter(double annealingParameter) {
    for (int category = 0; category < likelihoods.size(); category++) {
      graph(category).removeUnary(latestTip());
      model.buildObservation(latestTip(), context(category, annealingParameter));
    }
  }
  
  void addBranch(TreeNode x, TreeNode y, double length) {
    tree.addEdge(x, y, length);
    Pair<TreeNode, TreeNode> directedEdge = orient(x, y);
    for (int category = 0; category < nCategories(); category++) 
      model.buildTransition(directedEdge.getLeft(), directedEdge.getRight(), context(category, 1.0)); 
  }
  
  LikelihoodFactoryContext context(int category, double annealingParameter) {
    return new LikelihoodFactoryContext(likelihoods.get(category).graph, tree, observations, category, annealingParameter);
  }
  
  UndirectedGraph<TreeNode, ?> topology() {
    return tree.getTopology();
  }
  
  DiscreteFactorGraph<TreeNode> graph(int category) {
    return (DiscreteFactorGraph<TreeNode>) (likelihoods.get(category).graph);
  }
  
  TreeNode latestTip() { return likelihoods.get(0).latestTip; }
  
  void removeBranch(TreeNode v, TreeNode w) {
    tree.removeEdge(v, w);
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
}
