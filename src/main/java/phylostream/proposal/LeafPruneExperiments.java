package phylostream.proposal;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;

import bayonet.distributions.Random;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import blang.runtime.Observations;
import conifer.SequenceAlignment;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;


public class LeafPruneExperiments implements Runnable {
		
	@Arg       @DefaultValue("data/data27-1949.nex.run1.newick")
	public 	String tree = "data/data27-1949.nex.run1.newick";
	
	@Arg       @DefaultValue("data/M336_27.fasta")
	public String data = "data/M336_27.fasta";
	
	
	@Arg       @DefaultValue("results/likelihood.txt")
	public String out = "results/likelihood.txt";		

	@Override
	public void run() {

		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));

		Random rand = new Random(1);
		File file = new File(data);
		
		Observations obs = new Observations();
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);

		System.out.print(data);		

		LeafPrune leafPrune = new LeafPrune(urt, data);		

		Map<Pair<TreeNode,TreeNode>, Double> result = leafPrune.attachmentPointsLikelihoods(rand, 10);
		
		File outfile = new File(out);
		FileWriter fileWriter = null;
		try {
			fileWriter = new FileWriter(outfile);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


		for(Pair<TreeNode,TreeNode> edge : result.keySet()) {
			StringBuilder line = new StringBuilder();
			line.append(edge.getLeft());
			line.append('	');
			line.append(edge.getRight());
			line.append('	');
			line.append(result.get(edge));
			line.append("\n");
			try {
				fileWriter.write(line.toString());
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		try {
			fileWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		try {
			leafPrune.writePrunedSubtree("results/prunedTree.newick");
			leafPrune.writeTreeAfterPruning("results/treeAfterPruning.newick");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}

	public static void main(String[] args)
	{
		 LeafPruneExperiments leafPruneExperiments = new LeafPruneExperiments();
		 leafPruneExperiments.run();
	}

}

