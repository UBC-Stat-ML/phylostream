package phylostream;

import org.junit.Test;

import blang.validation.ExactInvarianceTest;
import conifer.SequenceAlignment;
import conifer.UnrootedTreeUtils;
import phylostream.blang.ReferenceModel;
import phylostream.io.Synthetic;
import phylostream.io.Synthetic.Realization;

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
    
    test.add(model, 
        (ReferenceModel m) -> UnrootedTreeUtils.totalTreeLength(m.getTree()));
  
    test.check();
  }
}
