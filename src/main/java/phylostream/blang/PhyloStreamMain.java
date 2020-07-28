package phylostream.blang;

import blang.core.ModelBuilder;
import blang.inits.experiments.Experiment;
import blang.runtime.Runner;

public class PhyloStreamMain extends Runner {

  public PhyloStreamMain(ModelBuilder builder) {
    super(builder);
  }
  
  public static void main(String [] args) {
    Experiment.startAutoExit(args);
  }

}
