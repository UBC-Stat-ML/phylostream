package phylostream.mcmc;

import java.util.List;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Random;
import bayonet.math.NumericalUtils;
import blang.core.Constrained;
import blang.mcmc.ConnectedFactor;
import blang.mcmc.MHSampler;
import blang.mcmc.SampledVariable;
import blang.mcmc.internals.Callback;
import conifer.TreeNode;
import phylostream.Utils;
import phylostream.blang.BayesianPhylo;

public class NNI extends MHSampler {
  
  @SampledVariable BayesianPhylo phylo;
  
  @ConnectedFactor Constrained ___;

  @Override
  public void propose(Random random, Callback callback) {
    if (numericFactors.size() != 2) throw new RuntimeException();
    
    // pick a branch near by to sample
    List<Pair<TreeNode, TreeNode>> possibleBranchMoves = Utils.possibleNNIs(phylo.tree);
    if (possibleBranchMoves.isEmpty()) return;
    
    final double initialLogDensity = phylo.logDensity();
    Pair<TreeNode, TreeNode> move = possibleBranchMoves.get(random.nextInt(possibleBranchMoves.size()));
    callback.setProposalLogRatio(0.0);
    Pair<TreeNode, TreeNode> undoMove = phylo.interchange(move.getLeft(), move.getRight());
    if (!callback.sampleAcceptance()) {
      phylo.interchange(undoMove.getLeft(), undoMove.getRight());
      NumericalUtils.checkIsClose(phylo.logDensity(), initialLogDensity);
    }
  }

}
