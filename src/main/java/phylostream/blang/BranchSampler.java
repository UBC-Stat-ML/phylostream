package phylostream.blang;

import bayonet.distributions.Random;
import blang.core.Constrained;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHSampler;
import blang.mcmc.SampledVariable;
import blang.mcmc.internals.Callback;

public class BranchSampler extends MHSampler {
  
  @SampledVariable BayesianPhylo tree;
  
  @ConnectedFactor Constrained ___;

  @Override
  public void propose(Random random, Callback callback) {
    System.out.println("Sampling branch length!");
  }

}
