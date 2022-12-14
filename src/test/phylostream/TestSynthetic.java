package phylostream;

import org.junit.Assert;
import org.junit.Test;

import conifer.SequenceAlignment;
import conifer.TreeNode;
import conifer.UnrootedTree;
import phylostream.io.Synthetic;
import phylostream.io.Synthetic.Realization;

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
  
  
  @Test
  public void generateTree() {
    Synthetic generator = new Synthetic();
    generator.nLeaves = 50;
    generator.nSites = 200;    
    Realization realization = generator.next();         
    UnrootedTree trueTree = realization.trueTree;
    TreeNode treeNode = realization.trueRoot;
    System.out.println(trueTree.toNewick());    
    System.out.println(treeNode);
  }

}
