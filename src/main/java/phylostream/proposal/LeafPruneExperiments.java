package phylostream.proposal;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

import org.eclipse.xtext.xbase.lib.Pair;

import bayonet.distributions.Random;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import blang.inits.GlobalArg;
import blang.inits.experiments.Experiment;
import blang.inits.experiments.ExperimentResults;
import blang.inits.experiments.tabwriters.TabularWriter;
import blang.runtime.Observations;
import conifer.SequenceAlignment;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;


public class LeafPruneExperiments extends Experiment {
	
	 
    @GlobalArg
	public ExperimentResults result;
		
	@Arg       @DefaultValue("data/data27-1949.nex.run1.newick")
	public 	String tree = "data/data27-1949.nex.run1.newick";
	
	@Arg       @DefaultValue("data/M336_27.fasta")
	public String data = "data/M336_27.fasta";
	
	@Arg       @DefaultValue("1")
	public int randSeed = 1;

	@Arg       @DefaultValue("10")
	public int nReplicates = 10;
	
	@Arg       @DefaultValue({"1", "2"})
	public List<Integer> neighborhoodRadius = Arrays.asList(1, 2);  
	

	@Override
	public void run() {
		
		Random rand = new Random(randSeed);
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		File file = new File(data);		
		Observations obs = new Observations();				
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);		
		LeafPrune leafPrune = new LeafPrune(urt, data);		
		
    Map<org.apache.commons.lang3.tuple.Pair<TreeNode, TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);
		
		TabularWriter allTV = result.getTabularWriter("totalVariation_all");
		TabularWriter edges = result.getTabularWriter("edgeNumberInNeighborhood");
		for(Integer k=0;k<neighborhoodRadius.size();k++)
		{
			//System.out.println(neighborhoodRadius.get(k));
      double[] tv0 = leafPrune.totalVariationSequenceNearestNeighbor(re, 20, neighborhoodRadius.get(k), 0.000001, 0.01);
			double acceptanceRate0 = leafPrune.getAcceptanceRate();
			double[] tvHalf = leafPrune.totalVariationSequenceNearestNeighbor(re, 20, neighborhoodRadius.get(k), 0.5, 0.01);
			double acceptanceRateHalf = leafPrune.getAcceptanceRate();
			double[] tv1 = leafPrune.totalVariationSequenceNearestNeighbor(re, 20, neighborhoodRadius.get(k), 1.0, 0.01);			
			TabularWriter tvFile = result.getTabularWriter("totalVariation_"+neighborhoodRadius.get(k));
			double acceptanceRate1 = leafPrune.getAcceptanceRate();
			for(int i=0; i<tvHalf.length; i++)
			{
			  List<Pair> toPrint = new ArrayList<>(Arrays.asList(
			      pair("neighborhoodRadius", neighborhoodRadius.get(k)),
            pair("log2n", i),
            pair("tv0", tv0[i]), 
            pair("tvHalf", tvHalf[i]),
            pair("tv1", tv1[i]),
            pair("neighborEdgeNumbers", leafPrune.getAverageEdgeNumberInNeighbor()),
            pair("acceptanceRate0", acceptanceRate0),
            pair("acceptanceRateHalf", acceptanceRateHalf),
            pair("acceptanceRate1", acceptanceRate1)));
				tvFile.write(toPrint.toArray(new Pair[0]));
				toPrint.add(pair("neighborRadius", k));
				allTV.write(toPrint.toArray(new Pair[0]));
			}	
			
	
			int [] edgeNumbers = leafPrune.getEdgeNumbers();

			for(int i=0; i<edgeNumbers.length; i++)
			{			
				edges.write(
						pair("neighborhoodRadius", neighborhoodRadius.get(k)),						
						pair("edgeIndex", i),
						pair("neighborEdgeNumbers", edgeNumbers[i])
						);
			}

		}
				
		UnrootedTree treeAfterPruning = leafPrune.urtAfterOneLeafRemoval();					
		try {
			result.getAutoClosedBufferedWriter("urt.newick").append(urt.toNewick());
			result.getAutoClosedBufferedWriter("prunedSubtree.newick").append(leafPrune.getPrunedSubtree());			
			result.getAutoClosedBufferedWriter("treeAfterPruning.newick").append(leafPrune.getTreeAfterPruning());
		} catch (IOException e) {
			e.printStackTrace();
		}
	    

		TabularWriter csv = result.getTabularWriter("treelikelihood");
		for(org.apache.commons.lang3.tuple.Pair<TreeNode, TreeNode> edge : re.keySet()) {			
			csv.write(
					pair("parent", edge.getLeft().toString()),
					pair("node", edge.getRight().toString()),
					pair("branch.length", treeAfterPruning.getBranchLength(edge.getLeft(), edge.getRight())), 
          pair("likelihood", re.get(edge).toString())
					);			
		}
		
		
	}

	public static void main(String[] args)
	{
	  //new LeafPruneExperiments().run();
		Experiment.start(args);
	}

	 private static <K,V> org.eclipse.xtext.xbase.lib.Pair<K,V> pair(K key ,V value) 
	  {
	    return org.eclipse.xtext.xbase.lib.Pair.of(key, value);
	  }
}

