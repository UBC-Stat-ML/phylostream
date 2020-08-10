package phylostream.blang;

import java.util.List;

import bayonet.smc.ParticlePopulation;
import blang.System;
import blang.engines.AdaptiveJarzynski;
import blang.engines.internals.PosteriorInferenceEngine;
import blang.inits.GlobalArg;
import blang.inits.experiments.ExperimentResults;
import blang.runtime.SampledModel;
import blang.runtime.internals.objectgraph.GraphAnalysis;
import conifer.TreeNode;
import conifer.io.TreeObservations;
import phylostream.io.Globals;

public class PhyloStreamSMC extends AdaptiveJarzynski implements PosteriorInferenceEngine {

  @GlobalArg public ExperimentResults results = new ExperimentResults();
  
  SampledModel model;
  
  @Override
  public void setSampledModel(SampledModel model) 
  { 
    this.model = model;
  }
  
  @Override
  public void performInference() 
  {
    System.out.indentWithTiming("initialCherry");
    ParticlePopulation<SampledModel> approximation = getApproximation(model);
    System.out.popIndent();
    
    TreeObservations observations = Globals.alignment;
    List<TreeNode> observedTreeNodes = observations.getObservedTreeNodes();
    
    for (int i = 2; i < observedTreeNodes.size(); i++) {
      TreeNode toAdd = observedTreeNodes.get(i);
      System.out.indentWithTiming(toAdd.toString());
      for (SampledModel sample : approximation.particles) {
        sample.setExponent(0.0);
        PhyloModel model = (PhyloModel) (sample.model);
        model.getTree().addTipContinuousUniform(toAdd, random);
      }
      approximation = getApproximation(approximation, 1.0, null, parallelRandomStreams, true);
      System.out.popIndent();
    }
  }
  
  @Override
  public void check(GraphAnalysis analysis) 
  {
    return;
  }
}
