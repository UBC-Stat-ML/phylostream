package phylostream.blang;

import java.io.File;
import java.util.Map;

import org.eclipse.xtext.xbase.lib.Pair;

import blang.inits.experiments.tabwriters.TabularWriter;
import blang.inits.experiments.tabwriters.TidySerializer;
import blang.inits.experiments.tabwriters.factories.CSV;
import blang.runtime.Runner;
import blang.runtime.internals.DefaultPostProcessor;
import briefj.BriefIO;
import conifer.UnrootedTree;
import conifer.UnrootedTreeUtils;
import conifer.UnrootedTreeUtils.TreeMetric;
import phylostream.io.Globals;

public class PhyloStreamPostProcessor extends DefaultPostProcessor {

  public static String GENERATING_TREE = "generatingTree.newick";
  
  @Override
  public void run() {
    super.run();
    if (Globals.generatingTree != null) {
      BriefIO.stringToFile(results.getFileInResultFolder(GENERATING_TREE), Globals.generatingTree.toNewick());
      
      reportSampleMetrics();
    }
  }
  
  /**
   * When the data is synthetic, compute E[d(t*, T)|data] where 
   * t* is the tree used to generated the data and d(.,.) is based 
   * on 3 standard metrics. 
   * Also plot that result for convenience
   */
  void reportSampleMetrics() {
    File samples = results.getFileInResultFolder(Runner.SAMPLES_FOLDER);
    File trees = CSV.csvFile(samples, "tree");
    String metricsFilePrefix = "sampleMetrics";
    TabularWriter sampleMetrics = results.getTabularWriter(metricsFilePrefix);
    for (Map<String,String> line : BriefIO.readLines(trees).indexCSV()) {
      String newickStr = line.get(TidySerializer.VALUE);
      UnrootedTree reconstruction = UnrootedTreeUtils.fromNewickString(newickStr);
      Map<TreeMetric, Double> metrics = UnrootedTreeUtils.computeTreeMetrics(Globals.generatingTree, reconstruction);
      TabularWriter current = sampleMetrics.child(Runner.sampleColumn, Integer.parseInt(line.get(Runner.sampleColumn)));
      for (TreeMetric metric : metrics.keySet())
        current.write(Pair.of("metric", metric), Pair.of(TidySerializer.VALUE, metrics.get(metric)));
    }
    sampleMetrics.close();
    File metricsFile = CSV.csvFile(results.resultsFolder, metricsFilePrefix);
    Map<String, Class<?>> types = TidySerializer.types(metricsFile);
    createPlot(
        new DensityPlot(metricsFile, types, this),
        results.resultsFolder
      );
  }
  
}
