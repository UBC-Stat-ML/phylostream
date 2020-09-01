package phylostream;


import org.junit.Test;

import blang.mcmc.internals.SamplerBuilderOptions;
import blang.validation.ExactInvarianceTest;
import blang.validation.Instance;
import conifer.SequenceAlignment;
import conifer.UnrootedTreeUtils;
import conifer.moves.SingleBranchScaling;
import conifer.moves.SingleNNI;
import phylostream.blang.ReferenceModel;
import phylostream.io.Synthetic;
import phylostream.io.Synthetic.Realization;
import phylostream.mcmc.BranchSlicerForReference;

public class TestReferenceModel {
  
  @SuppressWarnings("unchecked")
  @Test
  public void eit() {
    ExactInvarianceTest test = new ExactInvarianceTest();
    
    Synthetic synt = new Synthetic();
    synt.nSites = 20;
    synt.invariantSiteProbability = 0.0;
    synt.nPositiveCategories = 1;
    synt.nLeaves = 2;
    Realization realization = synt.next();
    SequenceAlignment alignment = (SequenceAlignment) realization.nextDataset();
    
    ReferenceModel model = new ReferenceModel.Builder().
      setData(synt).
      setObservations(alignment).
      setTree(realization.trueTree).
      setEvoModel(realization.trueModel).
        build(); 
    
    SamplerBuilderOptions samplerOptions = new SamplerBuilderOptions();
    samplerOptions.additional.add(BranchSlicerForReference.class);
    samplerOptions.excluded.add(SingleNNI.class);
    samplerOptions.excluded.add(SingleBranchScaling.class);
    Instance<ReferenceModel> instance = new Instance<ReferenceModel>(
        model, 
        samplerOptions , 
        (ReferenceModel m) -> UnrootedTreeUtils.totalTreeLength(m.getTree()));
        
    
    test.add(instance);
  
    test.check();
  }
}
