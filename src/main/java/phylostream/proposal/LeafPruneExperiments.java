package phylostream.proposal;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import org.apache.commons.lang3.tuple.Pair;
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
//import org.eclipse.xtext.xbase.lib.Pair;


public class LeafPruneExperiments extends Experiment {
	
	 
    @GlobalArg
	public ExperimentResults result;
		
	@Arg       @DefaultValue("data/data27-1949.nex.run1.newick")
	public 	String tree;
	
	@Arg       @DefaultValue("data/M336_27.fasta")
	public String data;
	
	@Arg       @DefaultValue("1")
	public int randSeed;

	@Arg       @DefaultValue("10")
	public int nReplicates;
	
	@Arg       @DefaultValue("1 2")
	public List<Integer> neighborhoodRadius;  
	

	@Override
	public void run() {
		
		Random rand = new Random(randSeed);
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		File file = new File(data);		
		Observations obs = new Observations();				
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);		
		LeafPrune leafPrune = new LeafPrune(urt, data);		
		Map<Pair<TreeNode,TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);
		
		TabularWriter edges = result.getTabularWriter("edgeNumberInNeighborhood");
		for(Integer k=0;k<neighborhoodRadius.size();k++)
		{
			System.out.println(neighborhoodRadius.get(k));
			double[] tv0 = leafPrune.totalVariationSequenceNearestNeighbor(re, 10, neighborhoodRadius.get(k), 0.000001);
			double acceptanceRate0 = leafPrune.getAcceptanceRate();
			double[] tvHalf = leafPrune.totalVariationSequenceNearestNeighbor(re, 10, neighborhoodRadius.get(k), 0.5);
			double acceptanceRateHalf = leafPrune.getAcceptanceRate();
			double[] tv1 = leafPrune.totalVariationSequenceNearestNeighbor(re, 10, neighborhoodRadius.get(k), 1.0);			
			TabularWriter tvFile = result.getTabularWriter("totalVariation_"+neighborhoodRadius.get(k));
			double acceptanceRate1 = leafPrune.getAcceptanceRate();
			for(int i=0; i<tvHalf.length; i++)
			{			
				tvFile.write(
						org.eclipse.xtext.xbase.lib.Pair.of("neighborhoodRadius", neighborhoodRadius.get(k)),
						org.eclipse.xtext.xbase.lib.Pair.of("log2n", i),
						org.eclipse.xtext.xbase.lib.Pair.of("tv0", tv0[i]), 
						org.eclipse.xtext.xbase.lib.Pair.of("tvHalf", tvHalf[i]),
						org.eclipse.xtext.xbase.lib.Pair.of("tv1", tv1[i]),
						org.eclipse.xtext.xbase.lib.Pair.of("neighborEdgeNumbers", leafPrune.getAverageEdgeNumberInNeighbor()),
						org.eclipse.xtext.xbase.lib.Pair.of("acceptanceRate0", acceptanceRate0),
						org.eclipse.xtext.xbase.lib.Pair.of("acceptanceRateHalf", acceptanceRateHalf),
						org.eclipse.xtext.xbase.lib.Pair.of("acceptanceRate1", acceptanceRate1)
						
						);
			}	
			
	
			int [] edgeNumbers = leafPrune.getEdgeNumbers();

			for(int i=0; i<edgeNumbers.length; i++)
			{			
				edges.write(
						org.eclipse.xtext.xbase.lib.Pair.of("neighborhoodRadius", neighborhoodRadius.get(k)),						
						org.eclipse.xtext.xbase.lib.Pair.of("edgeIndex", i),
						org.eclipse.xtext.xbase.lib.Pair.of("neighborEdgeNumbers", edgeNumbers[i])
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
		for(Pair<TreeNode,TreeNode> edge : re.keySet()) {			
			csv.write(
					org.eclipse.xtext.xbase.lib.Pair.of("parent", edge.getLeft().toString()),
					org.eclipse.xtext.xbase.lib.Pair.of("node", edge.getRight().toString()),
					org.eclipse.xtext.xbase.lib.Pair.of("branch.length", treeAfterPruning.getBranchLength(edge.getLeft(), edge.getRight())), 
                    org.eclipse.xtext.xbase.lib.Pair.of("likelihood", re.get(edge).toString())
					);			
		}
		
		
	}

	public static void main(String[] args)
	{
		Experiment.start(args);
	}

}

