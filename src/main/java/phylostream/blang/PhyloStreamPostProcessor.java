package phylostream.blang;

import java.io.File;
import java.util.Map;

import org.eclipse.xtext.xbase.lib.Pair;

import blang.inits.experiments.tabwriters.TabularWriter;
import blang.inits.experiments.tabwriters.factories.CSV;
import blang.runtime.Runner;
import blang.runtime.internals.DefaultPostProcessor;
import briefj.BriefIO;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;

public class PhyloStreamPostProcessor extends DefaultPostProcessor {

  @Override
  public void run() {
    File samplesFile = new File(blangExecutionDirectory.get(), Runner.SAMPLES_FOLDER);
    File treesCsv = CSV.csvFile(samplesFile, "tree");
    TabularWriter out = results.child("samples").getTabularWriter("branch");
    int i = 0;
    for (Map<String,String> line : BriefIO.readLines(treesCsv).indexCSV()) {
      String newick = line.get("value");
      UnrootedTree tree = UnrootedTreeUtils.fromNewickString(newick);
      double len = tree.getBranchLength(TreeNode.withLabel("synthetic-0"), TreeNode.withLabel("synthetic-1"));
      out.write(Pair.of("sample", i++), Pair.of("value", len));
    }
    out.close();
    
    super.run();
    
  }
  
}
