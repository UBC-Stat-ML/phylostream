package phylostream;

import org.junit.Assert;
import org.junit.Test;

import bayonet.distributions.Random;
import conifer.SequenceAlignment;

public class TestSynthetic {
  
  
  @Test
  public void generate() {
    Synthetic generator = new Synthetic();
    generator.nLeaves = 10000;
    generator.nSites = 200;
    Random rand = new Random(1);
    SequenceAlignment nextDataset = (SequenceAlignment) generator.next(rand).nextDataset(rand);
    Assert.assertEquals(nextDataset.getObservedTreeNodes().size(),  generator.nLeaves);
    Assert.assertEquals(nextDataset.nSites(),  generator.nSites);
  }
}
