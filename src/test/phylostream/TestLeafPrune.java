package phylostream;

import java.io.File;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.junit.Test;

import bayonet.distributions.Random;
import blang.runtime.Observations;
import conifer.SequenceAlignment;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;
import phylostream.proposal.LeafPrune;

public class TestLeafPrune {
  

	@Test
	public void testLeafPrune() {
				
		String tree = "data/data27-1949.nex.run1.newick";
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		
		Random rand = new Random(1);
		
		File file = new File("data/M336_27.fasta");
		Observations obs = new Observations();
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);
								
		System.out.print(data);		
				
		String rateMatrix = "/phylostream/COVID_GTR.txt";
		
		LeafPrune leafPrune = new LeafPrune(urt, data, rateMatrix);		
		
		Map<Pair<TreeNode,TreeNode>, Double> result = leafPrune.attachmentPointsLikelihoods(rand, 10);
		
		for(Pair<TreeNode,TreeNode> edge : result.keySet()) {
			
			System.out.println(edge+": "+result.get(edge));
		}
		
		leafPrune.expNormalizeLikelihoods(result);
		for(Pair<TreeNode,TreeNode> edge : result.keySet()) {
			
			System.out.println(edge+": "+result.get(edge));
		}		
						
	}
	
	
	@Test
	public void findNearestNeighborsOfAnEdge() {
				
		String tree = "data/data27-1949.nex.run1.newick";
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		
		Random rand = new Random(1);
		
		File file = new File("data/M336_27.fasta");
		Observations obs = new Observations();
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);
								
		String rateMatrix = "/phylostream/COVID_GTR.txt";
		LeafPrune leafPrune = new LeafPrune(urt, data, rateMatrix);		
		int nReplicates =5;
		Map<Pair<TreeNode,TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);		

		
		for(Pair<TreeNode,TreeNode> key : re.keySet())
		{
			Map<Pair<TreeNode,TreeNode>, Integer> neighborMap = leafPrune.neighborhoods(urt, key, 5);
			System.out.println("The neighbors of "+key.getLeft() + " "+ key.getRight()+" are: ");
			for(Pair<TreeNode,TreeNode> neighbor: neighborMap.keySet())
				System.out.println(neighbor.getLeft()+" "+ neighbor.getRight()+": "+ neighborMap.get(neighbor));
		}

						
	}
	
	
	@Test
	public void neighborDistance() {
		String tree = "data/data27-1949.nex.run1.newick";
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));		
		Random rand = new Random(1);		
		File file = new File("data/M336_27.fasta");
		Observations obs = new Observations();
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);
		String rateMatrix = "/phylostream/COVID_GTR.txt";
		LeafPrune leafPrune = new LeafPrune(urt, data, rateMatrix);		
		int nReplicates =5;
		Map<Pair<TreeNode,TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);		
//		double[][] distance = leafPrune.neighborDistance(urt, re, 5);
//		for(int i=0; i<distance.length; i++) {
//			for(int j=0; j<distance.length; j++)
//				System.out.print(distance[i][j]+" ");
//		System.out.println();
//		}
		
		double[][] neighborLikelihood = leafPrune.neighborLikelihood(urt, re, 5);
//		for(int i=0; i<neighborLikelihood.length; i++) {
//			for(int j=0; j<neighborLikelihood.length; j++)
//				System.out.print(neighborLikelihood[i][j]+" ");
//		System.out.println();
//		}				

		double power=0.5;
		double[][]  weights = leafPrune.neighborWeight(neighborLikelihood, power);
		for(int i=0; i<weights.length; i++) {
			for(int j=0; j<weights.length; j++)
				System.out.print(weights[i][j]+" ");
		System.out.println();
		}				

		System.out.println();
		System.out.println();
		System.out.println();
		
		power=0.001;
		weights = leafPrune.neighborWeight(neighborLikelihood, power);
		for(int i=0; i<weights.length; i++) {
			for(int j=0; j<weights.length; j++)
				System.out.print(weights[i][j]+" ");
		System.out.println();
		}
		
		System.out.println();
		System.out.println();
		System.out.println();
		
		power=1;
		weights = leafPrune.neighborWeight(neighborLikelihood, power);
		for(int i=0; i<weights.length; i++) {
			for(int j=0; j<weights.length; j++)
				System.out.print(weights[i][j]+" ");
		System.out.println();
		}				

		
	}
	
	
	
	
}
