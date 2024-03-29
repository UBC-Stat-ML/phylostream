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

	@Arg       @DefaultValue("10")
	public int K = 10;
	
	@Arg       @DefaultValue(" ")
	String rateMatrix = "/phylostream/COVID_GTR.txt";
	
	@Arg       @DefaultValue({"1", "2"})
	public List<Integer> neighborhoodRadius = Arrays.asList(1, 2);  
	
	
	
	@Override
	public void run() {
		
		Random rand = new Random(randSeed);
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		File file = new File(data);		
		Observations obs = new Observations();				
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);		
		LeafPrune leafPrune = new LeafPrune(urt, data, rateMatrix);		
		
        Map<org.apache.commons.lang3.tuple.Pair<TreeNode, TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);
		
		TabularWriter allTV = result.getTabularWriter("totalVariationDistance");
		TabularWriter edges = result.getTabularWriter("edgeNumberInNeighborhood");
		TabularWriter acc = result.getTabularWriter("acceptanceRate");
		
		
		double mixtureProportion = 1.0/urt.leaves().size(); 
		
		
		double[] powerList = new double[] {0, 0.5, 1};
		
		
		
		for(Integer k=0;k<neighborhoodRadius.size();k++)
		{		
			for(int l=0;l<powerList.length;l++)
			{
				double[] tv = leafPrune.totalVariationSequenceNearestNeighbor(re, K, neighborhoodRadius.get(k), powerList[l], mixtureProportion);
				double acceptanceRate = leafPrune.getAcceptanceRate();
				
				acc.write(
						pair("neighborRadius", neighborhoodRadius.get(k)),
						pair("Proposal", powerList[l]),						
						pair("acceptanceRate", acceptanceRate));
				
				for(int i=0; i<tv.length; i++)
				{
					List<Pair> toPrint = new ArrayList<>(Arrays.asList(
							pair("neighborRadius", neighborhoodRadius.get(k)),
							pair("Proposal", powerList[l]),
							pair("log2n", i),							
							pair("tv", tv[i])
							));
//					toPrint.add(pair("neighborRadius", k));
					allTV.write(toPrint.toArray(new Pair[0]));
				}	
			}
	
			int [] edgeNumbers = leafPrune.getEdgeNumbers();

			for(int i=0; i<edgeNumbers.length; i++)
			{			
				edges.write(
						pair("neighborRadius", neighborhoodRadius.get(k)),						
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

