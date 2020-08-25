package phylostream.blang;

import blang.runtime.internals.DefaultPostProcessor;
import briefj.BriefIO;
import phylostream.io.Globals;

public class PhyloStreamPostProcessor extends DefaultPostProcessor {

  public static String GENERATING_TREE = "generatingTree.newick";
  
  @Override
  public void run() {
    super.run();
    if (Globals.generatingTree != null)
      BriefIO.stringToFile(results.getFileInResultFolder(GENERATING_TREE), Globals.generatingTree.toNewick());
  }
  
}
