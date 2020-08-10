package phylostream.blang;

import phylostream.GrowingTree;
import phylostream.mcmc.BranchSlicer;
import phylostream.mcmc.NNI;
import blang.types.AnnealingParameter;
import briefj.collections.UnorderedPair;

import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;

import org.apache.commons.lang3.tuple.Pair;

import blang.core.RealDistribution;
import blang.mcmc.Samplers;
import blang.runtime.internals.objectgraph.SkipDependency;
import conifer.UnrootedTree;
import conifer.factors.NonClockTreePriorUtils;
import conifer.EvolutionaryModel;
import conifer.io.TreeObservations;
import conifer.TreeNode;

/**
 * State object used as basis of particles. Compared to GrowingTree (which could be used for max likelihood say), this 
 * has the following specialization:
 * 
 * - contains a prior
 * - keeps track of the tip to add automatically
 * - allows sampling from initial 2-leaves tree to initialize things
 */
@Samplers({BranchSlicer.class, NNI.class}) 
public class BayesianPhylo {
  @SkipDependency(isMutable = true)
  public GrowingTree tree;
  
  public final AnnealingParameter annealingParameter = new AnnealingParameter();
  final RealDistribution branchLengthPrior;
  final EvolutionaryModel model;
  
  @SkipDependency(isMutable = false)
  final TreeObservations observations;
  double logPrior;
  int nLeaves;
  
  /**
   * Create an un-initialized BayesianPhylo
   */
  public BayesianPhylo(EvolutionaryModel model, TreeObservations observations, RealDistribution branchLengthPrior) {
    this.tree = null;
    this.model = model;
    this.observations = observations;
    this.branchLengthPrior = branchLengthPrior;
    //recomputeTreeStatistics(unrootedTree); 
  }
  
  private void recomputeTreeStatistics(UnrootedTree unrootedTree) {
    logPrior = 0.0;
    for (double length : unrootedTree.getBranchLengths().values())
      logPrior += branchLengthPrior.logDensity(length);
    nLeaves = unrootedTree.leaves().size();
  }

  public void sampleInitialCherry(Random rand) {
    List<TreeNode> leaves = observations.getObservedTreeNodes().subList(0, 2);
    UnrootedTree unrootedTree = NonClockTreePriorUtils.sample(rand, branchLengthPrior, leaves);
    TreeNode root = leaves.get(0);
    this.tree = new GrowingTree(unrootedTree, root, model, observations);
    recomputeTreeStatistics(unrootedTree); 
    this.tree.setLatestTipAnnealingParameter(0.0);
  }
  
  public double logPrior() { return logPrior; }
  
  public double logLikelihood() {  
    tree.setLatestTipAnnealingParameter(annealingParameter.doubleValue());
    return tree.logLikelihood();
  }
  
  public double logDensity() { 
    return logPrior() + logLikelihood(); 
  }
  
  /**
   * Sample attachment point uniformly on the tree branches viewed as a continuous space.
   * 
   * Yields constant weight update when using an exponential prior on branch lengths (see pictures Aug 3).
   * Even in such case might need an update on the weights, which can be ignored 
   * if (1) rate of prior on branch length is fixed (2) we do not care about absolute 
   * marginal likelihood numerical values. 
   */
  public void addTipContinuousUniform(TreeNode freshLatestTip, Random rand) {
    // TODO: this makes it quadratic, but can be addressed by keeping these branch lengths in a 
    // binary tree. Will need to carefully update it so wait it becomes a bottleneck. 
    Map<UnorderedPair<TreeNode, TreeNode>, Double> branchLengths = tree.getBranchLengths();
    double totalBranchLength = branchLengths.values().stream().mapToDouble(Double::doubleValue).sum();
    double sum = 0.0;
    double unif = rand.nextDouble() * totalBranchLength;
    for (Entry<UnorderedPair<TreeNode, TreeNode>, Double> entry : branchLengths.entrySet()) {
      sum += entry.getValue();
      if (sum >= unif) {
        addTip(freshLatestTip, entry.getKey().getFirst(), entry.getKey().getSecond(), branchLengthPrior.sample(rand), rand.nextDouble() * entry.getValue());
        return;
      }
    }
    throw new RuntimeException();
  }
    
  public void addTip(TreeNode freshLatestTip, TreeNode v, TreeNode w, double l0, double l1) {
    // update prior:
    final double oldLength = tree.getBranchLength(v, w); // do this first!
    logPrior -= branchLengthPrior.logDensity(oldLength);
    logPrior += branchLengthPrior.logDensity(l1);
    logPrior += branchLengthPrior.logDensity(oldLength - l1);
    logPrior += branchLengthPrior.logDensity(l0);
    // update likelihood:
    tree.addTip(freshLatestTip, v, w, l0, l1);
  }
  
  public void updateBranchLength(TreeNode x, TreeNode y, double length) {
    // update prior
    final double oldLength = tree.getBranchLength(x, y); // do this first!
    logPrior -= branchLengthPrior.logDensity(oldLength);
    logPrior += branchLengthPrior.logDensity(length);
    // update likelihood:
    tree.updateBranchLength(x, y, length);
  }
  
  public Pair<TreeNode,TreeNode> interchange(TreeNode rootOfDisconnectedSubtree, TreeNode otherEndPointOfCutEdge) {
    // no change in prior
    // likelihood:
    return tree.interchange(rootOfDisconnectedSubtree, otherEndPointOfCutEdge);
  }

  @Override
  public String toString() {
    return tree.toString();
  }
}