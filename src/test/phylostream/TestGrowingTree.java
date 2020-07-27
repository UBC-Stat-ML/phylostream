package phylostream;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import org.junit.Assert;
import org.junit.Test;

import bayonet.distributions.Random;
import briefj.BriefCollections;
import briefj.collections.UnorderedPair;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import phylostream.Synthetic.Realization;

public class TestGrowingTree {
  
  /**
   * Return a synthetic tree. 
   */
  public GrowingTree getTree() { return getTree(0.005); }
  public GrowingTree getTree(double errorProbability) {
    Random rand = new Random(1);
    Synthetic generator = new Synthetic();
    generator.errorProbability = errorProbability;
    Realization realization = generator.next(rand);
    TreeObservations data = realization.nextDataset(rand);
    return new GrowingTree(realization.trueTree, realization.trueRoot, realization.trueModel, data);
  }

  @Test
  public void testInit() {
    check(getTree());
  }
  
  @Test
  public void testRootInvariance() {
    Random rand = new Random(1);
    Synthetic generator = new Synthetic();
    Realization realization = generator.next(rand);
    TreeObservations data = realization.nextDataset(rand);
    GrowingTree rooting1 = new GrowingTree(realization.trueTree, realization.trueRoot, realization.trueModel, data);
    final double reference = rooting1.logLikelihood();
    
    TreeNode anotherRoot = null;
    loop:for (TreeNode node : rooting1.tree.getTopology().vertexSet())
      if (node != anotherRoot) {
        anotherRoot = node;
        break loop;
      }
    
    GrowingTree rooting2 = new GrowingTree(realization.trueTree, anotherRoot, realization.trueModel, data);
        
    Assert.assertEquals(reference, rooting2.logLikelihood(), 1e-100); 
  }
  
  @Test
  public void testMoves() {
    Random rand = new Random(1);
    GrowingTree tree = getTree();
    
    for (int j = 0; j < 10; j++) {
      
      // perform some branch scalings
      for (UnorderedPair<TreeNode, TreeNode> edge : new ArrayList<>(tree.tree.getTopology().edgeSet())) {
        tree.updateBranchLength(edge.getFirst(), edge.getSecond(), rand.nextDouble());
        if (rand.nextBoolean())
          check(tree); 
      }
      
      // then some NNIs
      for (int i = 0; i < 100; i++) {
        List<Pair<TreeNode,TreeNode>> possibleNNIs = Utils.possibleNNIs(tree);
        Pair<TreeNode,TreeNode> selection = possibleNNIs.get(rand.nextInt(possibleNNIs.size()));
        tree.interchange(selection.getLeft(), selection.getRight());
        if (rand.nextBoolean())
          check(tree); 
      }
    }
    
  }

  
  @Test
  public void testAddTip() {
    GrowingTree tree = getTree();
    
    // full likelihood
    final double reference = tree.logLikelihood();
    
    // pick leaf
    UnrootedTree reducedTree = tree.tree;
    
    // record its position and remove it
    TreeNode leaf = reducedTree.leaves().get(1);
    TreeNode n0 = BriefCollections.pick(Graphs.neighborListOf(reducedTree.getTopology(), leaf));
    if (n0 == tree.root || leaf == tree.root) throw new RuntimeException("Test not designed to cover that case");
    final double pendantBranchLength = reducedTree.getBranchLength(leaf, n0);
    List<TreeNode> n1s = new ArrayList<>();
    for (TreeNode n1 : Graphs.neighborListOf(reducedTree.getTopology(), n0))
      if (n1 != leaf)
        n1s.add(n1); 
    final TreeNode v = n1s.get(0);
    final TreeNode w = n1s.get(1);
    final double v_n0_len = reducedTree.getBranchLength(n0, v);
    
    reducedTree.removeEdge(leaf, n0);
    reducedTree.getTopology().removeVertex(leaf);
    reducedTree.simplify();
    
    // reinsert it
    GrowingTree afterDeletion = new GrowingTree(reducedTree, tree.root, tree.model, tree.observations);
    afterDeletion.addTip(leaf, v, w, pendantBranchLength, v_n0_len);
    afterDeletion.setLatestTipAnnealingParameter(1.0);
    
    final double recomputed = afterDeletion.logLikelihood(); 
    
    Assert.assertEquals(reference, recomputed, 1e-100); 
  }
  
  @Test
  public void testErrorModelConvergence() {
    Assert.assertEquals(getTree(0.0).logLikelihood(), getTree(0.000000001).logLikelihood(), 1e-3); 
    Assert.assertNotEquals(getTree(0.0).logLikelihood(), getTree(0.01).logLikelihood(), 1e-3);
  }
  
  /**
   * Compare the logLikelihood computed from scratch from the one computed incrementally.
   */
  void check(GrowingTree tree) {
    TreeNode arbitraryRoot = TopologyUtils.arbitraryNode(tree.tree);
    LikelihoodComputationContext context = new LikelihoodComputationContext(
          EvolutionaryModelUtils.buildFactorGraphs(
            tree.model, 
            tree.tree, 
            arbitraryRoot, 
            tree.observations), 
          arbitraryRoot);
    final double reference = tree.model.computeLogLikelihood(context);
    final double incremental = tree.logLikelihood();
    Assert.assertEquals(reference, incremental, 1e-100); 
  }
}
