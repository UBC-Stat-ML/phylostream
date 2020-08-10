package phylostream;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import org.junit.Assert;
import org.junit.Test;

import bayonet.distributions.Random;
import blang.core.LogScaleFactor;
import blang.core.WritableRealVar;
import blang.mcmc.RealSliceSampler;
import briefj.BriefCollections;
import briefj.collections.UnorderedPair;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;
import phylostream.io.Synthetic;
import phylostream.io.Synthetic.Realization;

import static phylostream.Utils.check;

public class TestGrowingTree {
  
  /**
   * Return a synthetic tree. 
   */
  public GrowingTree getTree() { return getTree(0.005, 10); }
  public GrowingTree getTree(double errorProbability) { return getTree(errorProbability, 10); }
  public GrowingTree getTree(double errorProbability, int nLeaves) {
    Synthetic generator = new Synthetic();
    generator.nLeaves = nLeaves;
    generator.errorProbability = errorProbability;
    Realization realization = generator.next();
    TreeObservations data = realization.nextDataset();
    return new GrowingTree(realization.trueTree, realization.trueRoot, realization.trueModel, data);
  }

  @Test
  public void testInit() {
    check(getTree());
  }
  
  @Test
  public void testRootInvariance() {
    for (double anneal : Arrays.asList(0.0, 0.345, 1.0)) {
      Synthetic generator = new Synthetic();
      Realization realization = generator.next();
      TreeObservations data = realization.nextDataset();
      GrowingTree rooting1 = new GrowingTree(realization.trueTree, realization.trueRoot, realization.trueModel, data);
      rooting1.setLatestTipAnnealingParameter(anneal);
      final double reference = rooting1.logLikelihood();
      
      TreeNode anotherRoot = null;
      loop:for (TreeNode node : rooting1.unrootedTree.getTopology().vertexSet())
        if (node != anotherRoot) {
          anotherRoot = node;
          break loop;
        }
      
      GrowingTree rooting2 = new GrowingTree(realization.trueTree, anotherRoot, realization.trueModel, data);
      rooting2.setLatestTipAnnealingParameter(anneal);
          
      Assert.assertEquals(reference, rooting2.logLikelihood(), 1e-100); 
    }
  }
  
  @Test
  public void singleBranchTest() {
    GrowingTree tree = getTree(0.005, 2);
    tree.setLatestTipAnnealingParameter(0.0);
//    tree.setLatestTipAnnealingParameter(0.0);
    tree.setLatestTipAnnealingParameter(1.0);
    Random rand = new Random(1);
    Pair<TreeNode, TreeNode> pair = Utils.possibleBranchMoves(tree).get(0);
    for (int i = 0 ; i < 100; i++) {
      tree.updateBranchLength(pair.getLeft(), pair.getRight(), rand.nextDouble()); 
      
      if (rand.nextBoolean()) check(tree);
    }
  }
  
  @Test
  public void testMoves() {
    Random rand = new Random(1);
    GrowingTree tree = getTree();
    
    for (int j = 0; j < 10; j++) {
      // perform some branch scalings
      for (UnorderedPair<TreeNode, TreeNode> edge : new ArrayList<>(tree.unrootedTree.getTopology().edgeSet())) {
        
        for (int i = 0; i < 2; i++) {
          tree.updateBranchLength(edge.getFirst(), edge.getSecond(), rand.nextDouble());
          check(tree);
        }
        
        if (rand.nextBoolean())
          check(tree);
                
        WritableRealVar branchLen = new WritableRealVar() {
          double value = tree.getBranchLength(edge.getFirst(), edge.getSecond());
          @Override
          public double doubleValue() {
            return value;
          }
          @Override
          public void set(double value) {
            this.value = value;
            if (value >= 0.0)
              tree.updateBranchLength(edge.getFirst(), edge.getSecond(), value); 
          }
        };
        LogScaleFactor factor = () -> {
          if (branchLen.doubleValue() <= 0.0)
            return Double.NEGATIVE_INFINITY;
          else return tree.logLikelihood();
        };
        RealSliceSampler sampler = RealSliceSampler.build(branchLen, Arrays.asList(factor), 0.0, 1.0);
        sampler.execute(rand); 
        
      }
      
      // then some NNIs
      loop:for (int i = 0; i < 100; i++) {
        List<Pair<TreeNode,TreeNode>> possibleNNIs = Utils.possibleNNIs(tree);
        if (possibleNNIs.isEmpty()) break loop;
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
    UnrootedTree reducedTree = tree.unrootedTree;
    
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
  

}
