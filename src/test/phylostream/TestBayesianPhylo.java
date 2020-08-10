package phylostream;

import org.junit.Test;

import bayonet.distributions.Random;
import blang.distributions.Exponential;
import conifer.io.TreeObservations;
import phylostream.blang.BayesianPhylo;
import phylostream.io.Synthetic;
import phylostream.io.Synthetic.Realization;

public class TestBayesianPhylo {

  
  @Test
  public void testAdd() {
    Synthetic synt = new Synthetic();
    Realization generated = synt.next();
    TreeObservations dataset = generated.nextDataset();
    BayesianPhylo phylo = new BayesianPhylo(generated.trueModel, dataset, Exponential.distribution(() -> 10.0));
    phylo.annealingParameter._set(() -> 1.0);
    Random rand = new Random(1);
    phylo.sampleInitialCherry(rand);
    phylo.updateBranchLength(dataset.getObservedTreeNodes().get(0), dataset.getObservedTreeNodes().get(1), 1.1);
    phylo.addTipContinuousUniform(dataset.getObservedTreeNodes().get(2), rand); 
  }
}
