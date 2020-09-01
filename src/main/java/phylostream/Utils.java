package phylostream;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.eclipse.xtext.xbase.lib.Functions.Function1;
import org.jgrapht.Graphs;
import org.junit.Assert;

import briefj.collections.Counter;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.FastaUtils;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import viz.components.TreeViz;
import viz.core.Viz;

public class Utils {
  
  /**
   * Order in each pair coincides with order of arguments in GrowingTree.interchange
   */
  public static List<Pair<TreeNode,TreeNode>> possibleNNIs(GrowingTree tree) {
    List<Pair<TreeNode,TreeNode>> result = new ArrayList<>();
    TreeNode latestTipNeighbour = tree.latestEdge().getLeft();
    for (TreeNode other : Graphs.neighborListOf(tree.topology(), latestTipNeighbour))
      if (other != tree.latestTip())
        for (TreeNode root : Graphs.neighborListOf(tree.topology(), other))
          if (root != latestTipNeighbour)
            result.add(Pair.of(root, other)); 
    return result;
  }
  
  public static List<Pair<TreeNode,TreeNode>> possibleBranchMoves(GrowingTree tree) {
    return possibleBranchMoves(tree, 3);
  }
  
  public static List<Pair<TreeNode,TreeNode>> possibleBranchMoves(GrowingTree tree, int depth) {
    List<Pair<TreeNode,TreeNode>> result = new ArrayList<>();
    Pair<TreeNode, TreeNode> latestEdge = tree.latestEdge();
    possibleBranchMoves(tree, depth, result, latestEdge.getLeft(), latestEdge.getRight());
    return result;
  }
  
  private static void possibleBranchMoves(GrowingTree tree, int depth, List<Pair<TreeNode,TreeNode>> result, TreeNode parent, TreeNode current) {
    if (depth == 0) return;
    result.add(Pair.of(current, parent));
    for (TreeNode other : Graphs.neighborListOf(tree.topology(), current))
      if (other != parent)
        possibleBranchMoves(tree, depth - 1, result, current, other);
  }
  
  // Adapted from VasiliNovikov's in https://stackoverflow.com/questions/4965335/how-to-print-binary-tree-diagram
  public static String toString(UnrootedTree tree, TreeNode root) {
    StringBuilder buffer = new StringBuilder(50);
    toString(tree, root, buffer, "", "", null);
    return buffer.toString();
  }
  
  private static void toString(UnrootedTree tree, TreeNode node, StringBuilder buffer, String prefix, String childrenPrefix, TreeNode parent) {
    buffer.append(prefix);
    buffer.append(node.toString() + (parent == null ? "" : " [" + tree.getBranchLength(parent, node) + "]"));
    buffer.append('\n');
    for (Iterator<TreeNode> it = Graphs.neighborListOf(tree.getTopology(), node).iterator(); it.hasNext();) {
      TreeNode next = it.next();
      if (next != parent) {
        if (it.hasNext()) {
          toString(tree, next, buffer, childrenPrefix + "├── ", childrenPrefix + "│   ", node);
        } else {
          toString(tree, next, buffer, childrenPrefix + "└── ", childrenPrefix + "    ", node);
        }
      }
    }
  }
  
  public static void show(UnrootedTree tree) {
    show(tree, TopologyUtils.arbitraryNode(tree));
  }
  
  public static void show(UnrootedTree tree, TreeNode root) {
    Map<TreeNode,TreeNode> parents = GrowingTree.parentsMap(tree, root);
    Function1<? super TreeNode, ? extends List<TreeNode>> children = (TreeNode node) -> {
      List<TreeNode> result = new ArrayList<>();
      for (TreeNode neighbour : Graphs.neighborListOf(tree.getTopology(), node))
        if (parents.get(neighbour) == node)
          result.add(neighbour);
      return result;
    };
    Function1<? super TreeNode, ? extends Number> branchLengths = (TreeNode node) -> {
      TreeNode parent = parents.get(node);
      if (parent == null) return 0.0;
      return 1; //100.0 * tree.getBranchLength(UnorderedPair.of(node, parent));
    };
    TreeViz<TreeNode> viz = new TreeViz<TreeNode>(root, children, Viz.fixHeight(500), branchLengths);
    viz.show();
  }
  
  /**
   * Compare the logLikelihood computed from scratch from the one computed incrementally.
   */
  public static void check(GrowingTree tree) {
    if (tree.getLatestTipAnnealingParameter() != 1.0)
      throw new RuntimeException("This check only works with un-annealed tip.");
    TreeNode arbitraryRoot = TopologyUtils.arbitraryNode(tree.unrootedTree);
    LikelihoodComputationContext context = new LikelihoodComputationContext(
          EvolutionaryModelUtils.buildFactorGraphs(
            tree.model, 
            tree.unrootedTree, 
            arbitraryRoot, 
            tree.observations), 
          arbitraryRoot);
    final double reference = tree.model.computeLogLikelihood(context);
    final double incremental = tree.logLikelihood();
    Assert.assertEquals(reference, incremental, 1e-6); 
  }
  
  public static void main(String [] args) {
    LinkedHashMap<TreeNode, String> data = FastaUtils.readFasta(new File("data/covid19_BC_data-July-22.fasta.gz"));
    System.out.println("nTaxa = " + data.size());
    Counter<Character> letters = new Counter<Character>();
    List<LinkedHashSet<Character>> sets = new ArrayList<>();
    int size = -1;
    for (String str : data.values()) {
      if (size == -1) {
        size = str.length();
        System.out.println("nSites = " + size);
        for (int i = 0; i < size; i++)
          sets.add(new LinkedHashSet<>());
      }
      if (size != str.length()) throw new RuntimeException();
      for (int i = 0; i < size; i++) {
        letters.incrementCount(str.charAt(i), 1.0);
        if (str.charAt(i) != 'N')
          sets.get(i).add(str.charAt(i));
      }
    }
    Counter<Integer> sizeFreqs = new Counter<>();
    for (LinkedHashSet<Character> set : sets)
      sizeFreqs.incrementCount(set.size(), 1.0);
    System.out.println(sizeFreqs);
    System.out.println(letters);
  }
}
