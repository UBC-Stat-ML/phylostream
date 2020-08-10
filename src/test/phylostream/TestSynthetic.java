package phylostream;

import org.junit.Assert;
import org.junit.Test;

import conifer.SequenceAlignment;
import phylostream.io.Synthetic;

public class TestSynthetic {
  
  
  @Test
  public void generate() {
    Synthetic generator = new Synthetic();
    generator.nLeaves = 10000;
    generator.nSites = 200;
    SequenceAlignment nextDataset = (SequenceAlignment) generator.next().nextDataset();
    Assert.assertEquals(nextDataset.getObservedTreeNodes().size(),  generator.nLeaves);
    Assert.assertEquals(nextDataset.nSites(),  generator.nSites);
  }
}
