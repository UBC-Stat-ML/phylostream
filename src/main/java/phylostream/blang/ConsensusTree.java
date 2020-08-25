package phylostream.blang;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Map.Entry;
import binc.Command;
import briefj.BriefIO;
import briefj.BriefIO.ReadLineIterable;
import briefj.collections.Counter;
import briefj.collections.UnorderedPair;
import conifer.TreeNode;	
import conifer.UnrootedTree;

public class ConsensusTree {
	public static String neighborPath = "/Users/liangliangwang/Dropbox/software/phylip-3.695/exe//neighbor";
	private final Counter<UnorderedPair<TreeNode,TreeNode>> _meanDistances = new Counter<UnorderedPair<TreeNode,TreeNode>>();
	private double norm = 0.0;
	private List<TreeNode> taxa = null;
	private String trees  = "results/latest/samples/tree.csv";
	private String dir = "results/latest"; 
	public ConsensusTree(String trees){
		this.trees = trees;
	}

	public ConsensusTree(String trees, String dir){
		this.trees = trees;
		this.setDir(dir);	
	}

	public ConsensusTree(){
		this.trees = "results/latest/samples/tree.csv";		
	}

	public static void main(String[] args) {

		ConsensusTree consensus = new ConsensusTree(); 
		UnrootedTree urt = consensus.getConsensus();
		System.out.println(urt);
	}

	public UnrootedTree getConsensus()
	{   
		ReadLineIterable lineIterator = BriefIO.readLines(trees);		
		Iterator<String> iterator = lineIterator.iterator();
		while(iterator.hasNext()){
			String oneString = iterator.next().toString();
			String[] subStrings =  oneString.split("\"");
			if(subStrings.length>1) 
			{			
				process(UnrootedTree.fromNewickString(subStrings[1]), 1.0);
			}			
		}
		UnrootedTree njTree=inferTree(getUnrootedCladesPosterior());    
		return njTree;
	}

	public void process(UnrootedTree t, double w)
	{
		norm += w;
		_process(t, w);
	}

	private void _process(UnrootedTree ut, double w)
	{
		if (taxa == null)
			taxa = ut.leaves();

		Counter<UnorderedPair<TreeNode,TreeNode>> c = ut.allTotalBranchLengthDistances();

		synchronized (this)
		{
			for (Entry<UnorderedPair<TreeNode,TreeNode>,Double> keyPair : c.entries.entrySet())
			{
				double count = keyPair.getValue();
				UnorderedPair<TreeNode, TreeNode> key =keyPair.getKey();
				_meanDistances.incrementCount(key, w * count);
			}
		}
	}

	public Counter<UnorderedPair<TreeNode, TreeNode>> getUnrootedCladesPosterior()
	{
		Counter<UnorderedPair<TreeNode, TreeNode>> result = new Counter<UnorderedPair<TreeNode, TreeNode>>();
		for (UnorderedPair<TreeNode, TreeNode> key : _meanDistances.keySet())
			result.setCount(key, _meanDistances.getCount(key) / norm);
		return result;
	}
	
	public  UnrootedTree inferTree(Counter<UnorderedPair<TreeNode,TreeNode>> distances)
	{
		File directory = new File(this.dir, "concensusTree");		
		System.out.println("making a directory "+directory);
		boolean bool = directory.mkdir();
		if(bool){
			System.out.println("Directory concensusTree created successfully.");
		}else{
			System.out.println("Directory concensusTree cannot be created successfully!");
		}

		if (distances.keySet().size() == 1) // this situation makes NJ crash
		{
			File fileConsensus = new File(directory, "consensus.tree");
			UnrootedTree tree = inferTreeFromPair(distances);
			BriefIO.write(fileConsensus, tree.toNewick());
			return tree;
		}

		Map<TreeNode, TreeNode> conversion = new LinkedHashMap<TreeNode, TreeNode>();
		int i = 0;
		for (UnorderedPair<TreeNode, TreeNode> key : distances.keySet())
		{
			i = add(key.getFirst(), conversion, i);
			i = add(key.getSecond(), conversion, i);
		}
		Counter<UnorderedPair<TreeNode, TreeNode>> distances2 = new Counter<UnorderedPair<TreeNode, TreeNode>>();
		for (UnorderedPair<TreeNode, TreeNode> key : distances.keySet())
		{
			TreeNode c1 = conversion.get(key.getFirst());
			TreeNode c2 = conversion.get(key.getSecond());
			distances2.setCount(new UnorderedPair<TreeNode, TreeNode>(c1,c2), distances.getCount(key));
		}
		String str = phylipDistanceMatrix(distances2);
		File distFile = new File(directory, "infile");
		BriefIO.write(distFile, str);		
		File scriptFile = new File(directory, "script.bat");
		String scriptString = "#!bin/csh/\n"+"cd "+directory +"\n"+neighborPath+" <<EOF\n"+"infile\n"+"y\n"+"EOF\n"+"cd -";
		BriefIO.write(scriptFile, scriptString);
		Command cmd = Command.cmd("/bin/sh");
		Command.call(cmd.appendArgs(scriptFile.getAbsolutePath()));
		File fileBeforeConversion = new File(directory, "outtree");
		String newickStr = BriefIO.fileToString(fileBeforeConversion);	
		for (TreeNode originalName : conversion.keySet())
			newickStr = newickStr.replaceAll(conversion.get(originalName).toString(), originalName.toString());		
		File fileAfterConversion = new File(directory, "consensus.tree");
		BriefIO.write(fileAfterConversion, newickStr);
		return UnrootedTree.fromNewick(fileAfterConversion);		
	}


	public static String phylipDistanceMatrix(
			Counter<UnorderedPair<TreeNode, TreeNode>> pairwiseDistances)
	{
		StringBuilder result = new StringBuilder();
		Set<TreeNode> taxaSet = new LinkedHashSet();
		for (UnorderedPair<TreeNode, TreeNode> key : pairwiseDistances.keySet())
		{
			taxaSet.add(key.getFirst()); 
			taxaSet.add(key.getSecond());
		}
		List<TreeNode> taxa = new ArrayList<TreeNode>(taxaSet);
		result.append(taxa.size() + "\n");
		for (TreeNode treeNode : taxa)
		{
			result.append(cleanForPhylip(treeNode.toString()));
			for (TreeNode treeNode2 : taxa)
				//				result.append("  " + EasyFormat.fmt(pairwiseDistances.getCount(new UnorderedPair<TreeNode,TreeNode>(lang,lang2))));
				result.append("  " + pairwiseDistances.getCount(new UnorderedPair<TreeNode,TreeNode>(treeNode,treeNode2)));
			result.append('\n');
		}
		return result.toString();
	}



	public static String cleanForPhylip(String s)
	{
		return fillWithSpaces(cleanedTaxaName(s,true,10),10);
	}

	public static String fillWithSpaces(String s, int n)
	{
		if (s.length() > n) 
			return s.substring(0,n);
		int delta = n - s.length();
		for (int i = 0; i < delta; i++)
			s += " ";
		return s;
	}

	// clip taxon names to 10 character for the stupid phylip input format
	// Warning: this may create duplicates!
	public static String cleanedLangName(String s)
	{ return cleanedLangName(s, true); }
	// useful for removing parentheses and accents in names and optionally clipping to 10
	public static String cleanedLangName(String s, boolean clip)
	{
		return cleanedTaxaName(s, clip, 10);
	}
	public static String cleanedTaxaName(String s, boolean clip, int clipL)
	{
		// keep only letters
		String result = "";
		for (char c : s.toCharArray())
			if (("" + c).matches("[A-Za-z0-9]"))
				result += c;
		if (result.length() == 0)
			throw new RuntimeException("orig:" + s);
		// up to length 10 (for phylip)
		if (clip)
			return result.substring(0,Math.min(clipL,result.length())); 
		else return result;
	}


	
	public static UnrootedTree inferTreeFromPair(Counter<UnorderedPair<TreeNode, TreeNode>> distances)
	{	     
		UnorderedPair<TreeNode,TreeNode> key = distances.keySet().size() == 0 ? null : distances.keySet().iterator().next(); 
		double half = distances.getCount(key) /2; 
		String newickStr = "(" + key.getFirst() + ":" + half + "," + key.getSecond() + ":" + half + ");";
		return UnrootedTree.fromNewick(new File(newickStr)); 
	}

	private static int add(TreeNode taxon, Map<TreeNode, TreeNode> conversion, int i)
	{
		if (conversion.keySet().contains(taxon)) return i;
		conversion.put(taxon, TreeNode.withLabel("I"+i) );	
		return i+1;
	}


	public String getDir() {
		return dir;
	}

	public void setDir(String dir) {
		this.dir = dir;
	}

}
