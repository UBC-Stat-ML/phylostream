package phylostream.proposal;
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
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.TreeObservations;
import phylostream.io.Synthetic;
import phylostream.io.Synthetic.Realization;


public class LeafPruneSyntheticDataExperiments extends Experiment {
	
	 
    @GlobalArg
	public ExperimentResults result;
    
	@Arg       @DefaultValue("1")
	public int treeRandSeed = 1;

		
	@Arg       @DefaultValue("2")
	public int randSeed = 2;

	@Arg       @DefaultValue("10")
	public int nReplicates = 10;

	@Arg       @DefaultValue("10")
	public int K = 10;
	
	@Arg       @DefaultValue("/phylostream/COVID_GTR.txt")
	String rateMatrix = "/phylostream/COVID_GTR.txt";
	
	@Arg       @DefaultValue({"1", "2"})
	public List<Integer> neighborhoodRadius = Arrays.asList(1, 2);  
	
	
	@Arg       @DefaultValue("50")
	public int nLeaves = 50;
	

	@Arg       @DefaultValue("1000")
	public int nSites = 1000;
	
	@Arg       @DefaultValue("0.001")
	public double branchMeanLength = 0.001;  
	
	@Arg       @DefaultValue("0.5")	
	public double mixtureProportion =  0.5;
	
	@Arg       @DefaultValue("false")	
	public boolean globalMixture =  false;
	
	 @Arg       @DefaultValue("")
	 public 	String treeFile = "";

	
	@Override
	public void run() {
		Synthetic generator = new Synthetic();
		generator.rand =  new Random(treeRandSeed);
	    generator.nLeaves = nLeaves;
	    generator.nSites = 1000;	    
	    generator.branchMeanLength = branchMeanLength;
	    generator.errorProbability = 0;
	    generator.rateMatrix = rateMatrix;   
//	    generator.treeFile = "data/data27-1949.nex.run1.newick";	    
	    generator.treeFile = treeFile; //"/Users/liangliangwang/projects/phylostream/data/pruned.tree.nwk.txt";
	    

	    Realization realization = generator.next();         
	    UnrootedTree urt = realization.trueTree;
//	    TreeNode treeNode = realization.trueRoot;
	    	    
		try {
			result.getAutoClosedBufferedWriter("simulatedTree.newick").append(urt.toNewick());
		} catch (IOException e) {
			e.printStackTrace();
		}
	    	    		
		TreeObservations data =  generator.next().nextDataset();
	    System.out.println(data.nSites());
	    
		try {
			result.getAutoClosedBufferedWriter("data.txt").append(data.toString());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		Random rand = new Random(randSeed); 
		LeafPrune leafPrune = new LeafPrune(urt, data, rateMatrix);
		leafPrune.branchRate = 1.0/generator.branchMeanLength;
        Map<org.apache.commons.lang3.tuple.Pair<TreeNode, TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);

		
		UnrootedTree treeAfterPruning = leafPrune.urtAfterOneLeafRemoval();					
		try {
			result.getAutoClosedBufferedWriter("urt.newick").append(urt.toNewick());
			result.getAutoClosedBufferedWriter("prunedSubtree.newick").append(leafPrune.getPrunedSubtree());			
			result.getAutoClosedBufferedWriter("treeAfterPruning.newick").append(leafPrune.getTreeAfterPruning());
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		
		int treeRadius = leafPrune.treeRadius(treeAfterPruning); 	    
		 
		
		TabularWriter allTV = result.getTabularWriter("totalVariationDistance");
		TabularWriter edges = result.getTabularWriter("edgeNumberInNeighborhood");
		TabularWriter acc   = result.getTabularWriter("acceptanceRate");
		TabularWriter method = result.getTabularWriter("proposalMethod");
		TabularWriter treeradiusFile = result.getTabularWriter("treeRadius");
		treeradiusFile.write(
				pair("treeRadius", treeRadius)
				);
				
		double[] powerList = new double[] {0, 0.5, 1};
		
		for(Integer k=0;k<neighborhoodRadius.size();k++)
		{		
			for(int l=0;l<powerList.length;l++)
			{
				double[] tv = null;
				
				if(globalMixture && powerList[l]==0) {                                // Global uniform proposal
					System.out.println("global uniform:");
					tv = leafPrune.totalVariationGlobal(re, K, neighborhoodRadius.get(k), powerList[l], 1);
				}
				else
				{				
				if(globalMixture && powerList[l]>0)                                // Local proposal mixed with the global proposal                                        
					tv = leafPrune.totalVariationGlobal(re, K, neighborhoodRadius.get(k), powerList[l], mixtureProportion);
				else                 //  uniform local proposal, or local proposal mixed with the local uniform proposal  
					tv = leafPrune.totalVariationSequenceNearestNeighbor(re, K, neighborhoodRadius.get(k), powerList[l], mixtureProportion);
				}
				double acceptanceRate = leafPrune.getAcceptanceRate();
				
				acc.write(
						pair("neighborRadius", neighborhoodRadius.get(k)),
						pair("Proposal", powerList[l]),						
						pair("acceptanceRate", acceptanceRate));
				
				
				String methodName =  (globalMixture? "Global":"Local") + "; Proposal " + powerList[l] +"; "+ (powerList[l]==0?"":mixtureProportion);
				method.write(
						pair("Method", methodName),
						pair("MixtureType", (globalMixture? "Global":"Local")),
						pair("mixtureProportion", (powerList[l]==0?0:mixtureProportion)),									
						pair("Proposal", powerList[l])						
						);
								
				for(int i=0; i<tv.length; i++)
				{
					List<Pair> toPrint = new ArrayList<>(Arrays.asList(
							pair("neighborRadius", neighborhoodRadius.get(k)),							
							pair("mixtureProportion", (powerList[l]==0?0:mixtureProportion)),							
							pair("ProposalType", methodName),
							pair("LocalProposal", powerList[l]),
							pair("Global", (globalMixture? 1:0)),
							pair("log2n", i),							
							pair("tv", tv[i])
							));
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

