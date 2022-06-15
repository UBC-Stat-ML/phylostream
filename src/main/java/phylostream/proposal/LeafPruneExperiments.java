package phylostream.proposal;
import java.io.File;
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
	
	
	@Override
	public void run() {
		
		Random rand = new Random(randSeed);
		UnrootedTree urt = UnrootedTree.fromNewick(new File(tree));
		File file = new File(data);		
		Observations obs = new Observations();
		TreeObservations data = SequenceAlignment.loadObservedData(file, PhylogeneticObservationFactory.nucleotidesFactory(), obs);
//		System.out.print(data);		
		LeafPrune leafPrune = new LeafPrune(urt, data);		
		Map<Pair<TreeNode,TreeNode>, Double> re = leafPrune.attachmentPointsLikelihoods(rand, 5);
				
		TabularWriter csv = result.getTabularWriter("treelikelihood");
		for(Pair<TreeNode,TreeNode> edge : re.keySet()) {			
			csv.write(
					org.eclipse.xtext.xbase.lib.Pair.of("node1", edge.getLeft().toString()),
					org.eclipse.xtext.xbase.lib.Pair.of("node2", edge.getRight().toString()),
					org.eclipse.xtext.xbase.lib.Pair.of("likelihood", re.get(edge).toString())
					);			
		}
		
		
	}

	public static void main(String[] args)
	{
		Experiment.start(args);
	}

}

