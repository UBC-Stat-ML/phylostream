package phylostream.blang;

import java.io.File;
import blang.runtime.internals.DefaultPostProcessor;
import briefj.BriefIO;
import phylostream.io.Globals;

public class PhyloStreamPostProcessor extends DefaultPostProcessor {

	public static String GENERATING_TREE = "generatingTree.newick";
	public static String neighborPath = "/Users/liangliangwang/Dropbox/software/phylip-3.695/exe//neighbor";
	public static String trees = "results/latest/samples/tree.csv";
	@Override
	public void run() {
		super.run();
		if (Globals.generatingTree != null)
			BriefIO.stringToFile(results.getFileInResultFolder(GENERATING_TREE), Globals.generatingTree.toNewick());				
		File _get_1 = this.blangExecutionDirectory.get();
		System.out.println(_get_1);
		ConsensusTree consensus = new ConsensusTree(trees, _get_1.toString()); 		 
		consensus.getConsensus();		
	}
	
	public static void main(String[] args){
		(new PhyloStreamPostProcessor()).run();
	}
		
}
