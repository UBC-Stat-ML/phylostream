package phylostream;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import bayonet.distributions.Random;
import blang.core.RealConstant;
import blang.distributions.Exponential;
import blang.runtime.internals.objectgraph.DeepCloner;
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
  
  @Test
  public void testClone() {
    
    Synthetic synt = new Synthetic();
    synt.nLeaves = 100;
    Realization generated = synt.next();
    TreeObservations dataset = generated.nextDataset();
    BayesianPhylo phylo = new BayesianPhylo(generated.trueModel, dataset, Exponential.distribution(new RealConstant(10.0)));
    phylo.annealingParameter._set(new RealConstant(1.0));
    Random rand = new Random(1);
    phylo.sampleInitialCherry(rand);
    phylo.updateBranchLength(dataset.getObservedTreeNodes().get(0), dataset.getObservedTreeNodes().get(1), 1.1);
    
    for (int i = 2; i < synt.nLeaves; i++)
      phylo.addTipContinuousUniform(dataset.getObservedTreeNodes().get(i), rand); 

    List<Object> keepThem = new ArrayList<>();
    for (int i = 0; i < 10; i++) {
      System.out.println(i);
      keepThem.add(DeepCloner.deepClone(phylo));
    }
  }
}
