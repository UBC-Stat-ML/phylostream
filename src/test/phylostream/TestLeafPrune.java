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
				
		LeafPrune leafPrune = new LeafPrune(urt, data);		
		
		Map<Pair<TreeNode,TreeNode>, Double> result = leafPrune.attachmentPointsLikelihoods(rand, 10);
		
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
								
		
		LeafPrune leafPrune = new LeafPrune(urt, data);		
		int nReplicates =5;
		Map<Pair<TreeNode,TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, nReplicates);		

		
		for(Pair<TreeNode,TreeNode> key : re.keySet())
		{
			List<Pair<TreeNode,TreeNode>> neighborList = leafPrune.neighborhoods(urt, key);
			System.out.println("The neighbors of "+key.getLeft() + " "+ key.getRight()+" are: ");
			for(Pair<TreeNode,TreeNode> neighbor: neighborList)
				System.out.println(neighbor.getLeft()+" "+ neighbor.getRight());
		}

						
	}
}
