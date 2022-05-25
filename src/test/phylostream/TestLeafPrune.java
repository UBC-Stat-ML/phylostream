package phylostream;

import java.io.File;
import org.junit.Test;

import bayonet.distributions.Random;
import blang.runtime.Observations;
import conifer.SequenceAlignment;
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
		leafPrune.propose(rand);
	}
	
}
