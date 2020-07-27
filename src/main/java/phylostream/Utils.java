package phylostream;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.eclipse.xtext.xbase.lib.Functions.Function1;
import org.jgrapht.Graphs;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
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
}
