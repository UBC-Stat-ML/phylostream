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

	
	@Override
	public void run() {
		
		Random rand = new Random(randSeed);
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		File file = new File(data);		
		Observations obs = new Observations();				
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);		
		LeafPrune leafPrune = new LeafPrune(urt, data);		
		Map<Pair<TreeNode,TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);		
		double[] loglikelihoodsVec = leafPrune.loglikelihoodsVec(re);
		
//		boolean normalized = leafPrune.stationaryDist(loglikelihoodsVec);
//		System.out.println(normalized);
//		for(int i=0;i<loglikelihoodsVec.length;i++)
//		{
//			System.out.println(loglikelihoodsVec[i]);
//		}

		int[] neighborhoodRadius = new int[] {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,16,17,18,19,20,21,22,23,24,25};
		
		for(int k=0;k<neighborhoodRadius.length;k++)
		{
			double[] tv0 = leafPrune.totalVariationSequenceNearestNeighbor(re, 10, neighborhoodRadius[k], 0.000001);
			double[] tvHalf = leafPrune.totalVariationSequenceNearestNeighbor(re, 10, neighborhoodRadius[k], 0.5);			
			double[] tv1 = leafPrune.totalVariationSequenceNearestNeighbor(re, 10, neighborhoodRadius[k], 1.0);
			TabularWriter tvFile = result.getTabularWriter("totalVariation_"+neighborhoodRadius[k]);
			for(int i=0; i<tvHalf.length; i++)
			{			
				tvFile.write(
						org.eclipse.xtext.xbase.lib.Pair.of("neighborhoodRadius", neighborhoodRadius[k]),
						org.eclipse.xtext.xbase.lib.Pair.of("log2n", i),
						org.eclipse.xtext.xbase.lib.Pair.of("tv0", tv0[i]), 
						org.eclipse.xtext.xbase.lib.Pair.of("tvHalf", tvHalf[i]),
						org.eclipse.xtext.xbase.lib.Pair.of("tv1", tv1[i])
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

