package phylostream;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import org.eclipse.xtext.xbase.lib.Functions.Function1;
import org.jgrapht.Graphs;

import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import viz.components.TreeViz;
import viz.core.Viz;

public class Utils {
  
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
