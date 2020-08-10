package phylostream.mcmc;

import java.util.Arrays;
import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Random;
import blang.core.Constrained;
import blang.core.LogScaleFactor;
import blang.core.WritableRealVar;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.RealSliceSampler;
import blang.mcmc.SampledVariable;
import blang.mcmc.Sampler;
import conifer.TreeNode;
import phylostream.Utils;
import phylostream.blang.BayesianPhylo;

public class BranchSlicer implements Sampler {
  
  @SampledVariable BayesianPhylo phylo;
  
  @ConnectedFactor Constrained ___;
  
  @ConnectedFactor List<LogScaleFactor> __;

  @Override
  public void execute(Random rand) {
    // pick a branch near by to sample
    List<Pair<TreeNode, TreeNode>> possibleBranchMoves = Utils.possibleBranchMoves(phylo.tree);
    Pair<TreeNode, TreeNode> move = possibleBranchMoves.get(rand.nextInt(possibleBranchMoves.size()));
    
    // form a slice sampler for that edge
    WritableRealVar branchLen = new WritableRealVar() {
      double value = phylo.tree.getBranchLength(move.getLeft(), move.getRight());
      @Override
      public double doubleValue() {
        return value;
      }
      @Override
      public void set(double value) {
        this.value = value;
        if (value >= 0.0)
          phylo.updateBranchLength(move.getLeft(), move.getRight(), value); 
      }
    };
    LogScaleFactor factor = () -> {
      if (branchLen.doubleValue() <= 0.0)
        return Double.NEGATIVE_INFINITY;
      else return phylo.logDensity();
    };
    RealSliceSampler sampler = RealSliceSampler.build(branchLen, Arrays.asList(factor));
    sampler.execute(rand); 
  }

}
