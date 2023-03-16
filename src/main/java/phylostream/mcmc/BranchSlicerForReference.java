package phylostream.mcmc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import bayonet.distributions.Random;
import blang.core.LogScaleFactor;
import blang.core.WritableRealVar;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.RealSliceSampler;
import blang.mcmc.SampledVariable;
import blang.mcmc.Sampler;
import briefj.collections.UnorderedPair;
import conifer.TreeNode;
import conifer.UnrootedTree;
import static conifer.Utils.*;


public class BranchSlicerForReference implements Sampler {
  
  @SampledVariable UnrootedTree variable;
    
  @ConnectedFactor List<LogScaleFactor> numericFactors;

  @Override
  public void execute(Random rand) {
    // pick a branch near by to sample
    ArrayList<UnorderedPair<TreeNode,TreeNode>> allEdges = new ArrayList<>(variable.getTopology().edgeSet());
    final UnorderedPair<TreeNode, TreeNode> move = sample(allEdges, rand);
    
    // form a slice sampler for that edge
    WritableRealVar branchLen = new WritableRealVar() {
      double value = variable.getBranchLength(move.getFirst(), move.getSecond());
      @Override
      public double doubleValue() {
        return value;
      }
      @Override
      public void set(double value) {
        this.value = value;
        if (value >= 0.0)
          variable.updateBranchLength(move, value); 
      }
    };
    LogScaleFactor factor = () -> {
      if (branchLen.doubleValue() <= 0.0)
        return Double.NEGATIVE_INFINITY;
      else return ld();
    };
    RealSliceSampler sampler = RealSliceSampler.build(branchLen, Arrays.asList(factor));
    sampler.execute(rand); 
  }
  
  private double ld() {
    double sum = 0.0;
    for (LogScaleFactor f : numericFactors)
      sum += f.logDensity();
    return sum;
  }

}
