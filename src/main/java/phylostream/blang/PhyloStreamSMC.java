package phylostream.blang;

import java.util.List;

import org.eclipse.xtext.xbase.lib.Pair;

import bayonet.smc.ParticlePopulation;
import blang.System;
import blang.engines.AdaptiveJarzynski;
import blang.engines.internals.PosteriorInferenceEngine;
import blang.inits.GlobalArg;
import blang.inits.experiments.ExperimentResults;
import blang.io.BlangTidySerializer;
import blang.runtime.Runner;
import blang.runtime.SampledModel;
import blang.runtime.internals.objectgraph.GraphAnalysis;
import conifer.TreeNode;
import conifer.io.TreeObservations;
import phylostream.io.Globals;

import static blang.engines.internals.factories.SCM.*;

public class PhyloStreamSMC extends AdaptiveJarzynski implements PosteriorInferenceEngine {

  @GlobalArg public ExperimentResults results = new ExperimentResults();
  
  SampledModel model;
  
  @Override
  public void setSampledModel(SampledModel model) 
  { 
    this.model = model;
  }
  
  int tipId = -1;
  int globalIteration = -1;
  
  @SuppressWarnings("unchecked")
  @Override
  public void performInference() 
  {
    globalIteration = 0;
    System.out.indentWithTiming("initialCherry");
    ParticlePopulation<SampledModel> approximation = getApproximation(model);
    System.out.popIndent();
    
    TreeObservations observations = Globals.alignment;
    List<TreeNode> observedTreeNodes = observations.getObservedTreeNodes();
    
    for (tipId = 2; tipId < observedTreeNodes.size(); tipId++) {
      TreeNode toAdd = observedTreeNodes.get(tipId);
      System.out.indentWithTiming(toAdd.toString());
      for (SampledModel sample : approximation.particles) {
        sample.setExponent(0.0);
        PhyloModel model = (PhyloModel) (sample.model);
        model.getTree().addTipContinuousUniform(toAdd, random);
      }
      approximation = getApproximation(approximation, 1.0, null, parallelRandomStreams, true);
      System.out.popIndent();
    }
    
    approximation = approximation.resample(random, resamplingScheme);
    
    BlangTidySerializer tidySerializer = new BlangTidySerializer(results.child(Runner.SAMPLES_FOLDER)); 
    BlangTidySerializer densitySerializer = new BlangTidySerializer(results.child(Runner.SAMPLES_FOLDER)); 
    int particleIndex = 0;
    for (SampledModel model : approximation.particles)  
    {
      model.getSampleWriter(tidySerializer).write(Pair.of(Runner.sampleColumn, particleIndex)); 
      densitySerializer.serialize(model.logDensity(), "logDensity", Pair.of(Runner.sampleColumn, particleIndex));
      particleIndex++;
    }
  }
  
  public static final String
  
    tipIdColumn = "tipId";
  
  @Override
  protected void recordPropagationStatistics(int iteration, double temperature, double ess) {
    results.child(Runner.MONITORING_FOLDER).getTabularWriter(propagationFileName).printAndWrite(
        Pair.of(iterationColumn, globalIteration++),
        Pair.of(tipIdColumn, tipId),
        Pair.of(annealingParameterColumn, temperature),
        Pair.of(essColumn, ess)
    );
  }

  @Override
  protected void recordResamplingStatistics(int iteration, double temperature, double logNormalization) {
    results.child(Runner.MONITORING_FOLDER).getTabularWriter(resamplingFileName).printAndWrite(
        Pair.of(iterationColumn, globalIteration),
        Pair.of(tipIdColumn, tipId),
        Pair.of(annealingParameterColumn, temperature),
        Pair.of(logNormalizationColumn, logNormalization)
    );
  }
  
  @Override
  public void check(GraphAnalysis analysis) 
  {
    return;
  }
}
