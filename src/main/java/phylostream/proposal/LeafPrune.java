package phylostream.proposal;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.Graphs;
import com.google.common.collect.Lists;
import bayonet.distributions.Random;
import bayonet.graphs.GraphUtils;
import bayonet.marginal.FactorGraph;
import bayonet.marginal.UnaryFactor;
import bayonet.marginal.algo.SumProduct;
import bayonet.math.NumericalUtils;
import blang.core.RealConstant;
import blang.core.RealDistribution;
import blang.distributions.Exponential;
import briefj.collections.UnorderedPair;
import conifer.EvolutionaryModel;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.Utils;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.SimpleRateMatrix;
import conifer.io.TreeObservations;
import conifer.models.DiscreteGammaMixture;
import conifer.models.EvolutionaryModelUtils;
import conifer.models.LikelihoodComputationContext;
import conifer.models.MultiCategorySubstitutionModel;
import phylostream.io.Synthetic;


public class LeafPrune 
{

	private UnrootedTree urtree;  
	private TreeObservations data;
	private List<SumProduct<TreeNode>> likelihoods; 
	private EvolutionaryModel evolutionaryModel;
	private UnrootedTree prunedSubtree = null; 
	private UnrootedTree urtAfterOneLeafRemoval = null;   // tree after pruning a leaf

	public LeafPrune(UnrootedTree urt, TreeObservations data) {
		this.urtree = urt;
		this.data = data;	  
		TreeNode arbitraryRoot = Synthetic.arbitraryRoot(urt);
		String rateMatrix = "/phylostream/COVID_GTR.txt";	    
		double [][] loadedRateMatrix = SimpleRateMatrix.fromResource(rateMatrix).getRateMatrix();	    
		CTMCParameters ctmc = new SimpleRateMatrix(loadedRateMatrix, null);	    
		double invariantSiteProbability = 0.01;	    	    
		double shapeParameter = 1;
		int nPositiveCategories = 4;
		int nSites = data.nSites();	    	    
		DiscreteGammaMixture rateMatrixMixture = new DiscreteGammaMixture(new RealConstant(invariantSiteProbability), new RealConstant(shapeParameter), ctmc, nPositiveCategories);
		this.evolutionaryModel = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(rateMatrixMixture, nSites);		
		List<FactorGraph<TreeNode>> graphs = EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, urt, arbitraryRoot, data);				
		this.likelihoods = new ArrayList<SumProduct<TreeNode>>(graphs.size());
		for (int category = 0; category < graphs.size(); category++) 
			this.likelihoods.add(new SumProduct<TreeNode>(graphs.get(category)));	 	  	  
	}



	public Pair<TreeNode, TreeNode>  edgeAttachedToLeaf(Random random, UnrootedTree urtree) {
		// nothing interesting to do if tree is only a single branch
		List<TreeNode> internalNodes = GraphUtils.internalNodes(urtree.getTopology());
		if (internalNodes.isEmpty())
			throw new RuntimeException("There are no internal nodes.");

		// get the leaves
		List<TreeNode> leavesNodes = GraphUtils.leaves(urtree.getTopology()); 

		// pick a leaf node at random	 	    
		TreeNode removedLeaf = Utils.sample(leavesNodes, random);

		// one neighbor will from the edge to disconnect, and another one, the new root
		List<UnorderedPair<TreeNode,TreeNode>> neighbors = Lists.newArrayList(urtree.getTopology().edgesOf(removedLeaf));


		if (neighbors.size() != 1)
			throw new RuntimeException("Currently only supporting leaf (only 1 edge connecting to it).");

		TreeNode  removedInternalNode  = GraphUtils.pickOther(neighbors.get(0),removedLeaf);

		return Pair.of(removedInternalNode, removedLeaf);
	}


	
	
	
	public String getPrunedSubtree() {
		return prunedSubtree.toNewick(); 
	}
	

	public String getTreeAfterPruning() {
		return urtAfterOneLeafRemoval.toNewick(); 
	}

		

	public Map<Pair<TreeNode,TreeNode>, Double>  attachmentPointsLikelihoods(Random random, int N) {

		Map<Pair<TreeNode,TreeNode>, Double> result =  new HashMap<Pair<TreeNode, TreeNode>,Double>();	
		 
		UnrootedTree urt = new UnrootedTree(urtree);

		Pair<TreeNode, TreeNode> randomEdge = edgeAttachedToLeaf(random, urt);
		TreeNode removedInternalNode = randomEdge.getLeft(), removedLeaf=randomEdge.getRight();

		// disconnect the leaf from the tree
		prunedSubtree = pruneLeaf(urt, removedLeaf, removedInternalNode); 
//		System.out.println(prunedSubtree);
		
		urtAfterOneLeafRemoval = new UnrootedTree(urt);

		List<List<UnaryFactor<TreeNode>>> prunedSubtreeMarginals = Lists.newArrayList();

		double  rate = 10; 
		RealDistribution brDist = Exponential.distribution(new RealConstant(rate));

		double br[] = new double[N]; 

		for (int i=0;i<N;i++)
		{
			br[i] = brDist.sample(random);
			//			System.out.println(br);
			prunedSubtree.updateBranchLength(prunedSubtree.getTopology().getEdge(removedLeaf, removedInternalNode), br[i]);	
			prunedSubtreeMarginals.add(EvolutionaryModelUtils.getRootMarginalsFromFactorGraphs(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, prunedSubtree, removedInternalNode, data), removedInternalNode));

		}


		Map<Pair<TreeNode, TreeNode>, List<TreeNode>>  attachmentPointsMap = addMultipleAuxiliaryInternalNodes(urt, removedInternalNode, N);

		// run the sum product on the main tree
		List<SumProduct<TreeNode>> mainTreeSumProducts = EvolutionaryModelUtils.getSumProductsFromFactorGraphs(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, urt, removedInternalNode, data, false), removedInternalNode);

		for (Pair<TreeNode, TreeNode> edge:attachmentPointsMap.keySet()) {

			List<TreeNode> attachmentPoints = attachmentPointsMap.get(edge);
			
			double edgeLoglikelihood = Double.NEGATIVE_INFINITY;  

			for (int i = 0; i < attachmentPoints.size(); i++)
			{

				//	      System.out.println("i is "+i);
				TreeNode attachmentPoint = attachmentPoints.get(i);
				//	      System.out.println(attachmentPoint);
				List<List<UnaryFactor<TreeNode>>> currentFullUnaries = Lists.newArrayList();
				for (int j=0; j<N; j++)  currentFullUnaries.add(Lists.newArrayList());

				for (int f = 0; f < mainTreeSumProducts.size(); f++)
				{
					SumProduct<TreeNode> currentMainTreeSumProduct = mainTreeSumProducts.get(f);
					//	        System.out.println("currentMainTreeSumProduct log normal "+currentMainTreeSumProduct.logNormalization());

					UnaryFactor<TreeNode> currentMainTreeMarginal = currentMainTreeSumProduct.computeMarginal(attachmentPoint);

					//	        System.out.println(currentMainTreeMarginal.logNormalization());

					for (int j=0; j<N; j++)
					{
						UnaryFactor<TreeNode> prunedSubtreeMarginal = prunedSubtreeMarginals.get(j).get(f);
						currentFullUnaries.get(j).add(currentMainTreeSumProduct.getFactorGraph().factorOperations().pointwiseProduct(Arrays.asList(prunedSubtreeMarginal, currentMainTreeMarginal)));
					}
				}

				double branchLoglikelihood = Double.NEGATIVE_INFINITY; 
				for (int j=0; j<N; j++) {
					LikelihoodComputationContext context = new LikelihoodComputationContext(currentFullUnaries.get(j));		      
					branchLoglikelihood = NumericalUtils.logAdd(branchLoglikelihood, evolutionaryModel.computeLogLikelihood(context)- brDist.logDensity(br[j]));
				}
				edgeLoglikelihood = NumericalUtils.logAdd(edgeLoglikelihood, branchLoglikelihood-Math.log(N));
			}
			
			result.put(edge, edgeLoglikelihood-Math.log(attachmentPoints.size()));
		}

		return result; 	    
	}


	/**
	 * This splits the tree into two parts relative to the edge e=(removedLeaf, detachedNode).
	 * This will return a subtree containing e and the leaf. 
	 * The unrootedTree urtree will be modified to remove all the edges and nodes in the returned
	 * tree, except for the node detachedNode.  
	 */
	public UnrootedTree pruneLeaf(UnrootedTree urtree, final TreeNode removedLeaf, final TreeNode detachedNode) {
		UnrootedTree result = new UnrootedTree();
		result.addNode(removedLeaf);
		result.addNode(detachedNode);
		result.addEdge(removedLeaf, detachedNode, urtree.getBranchLength(detachedNode, removedLeaf));

		urtree.removeEdge(detachedNode,removedLeaf);
		urtree.getTopology().removeVertex(removedLeaf);

		List<TreeNode> neighbors = Graphs.neighborListOf(urtree.getTopology(), detachedNode);
		if (neighbors.size() != 2)
			throw new RuntimeException();
		return result;
	}



	/**
	 * Iterate the edge (oriented with the provided root) and add a dummy internal node on
	 * each edge, except for edges connected to current.
	 * The nodes are placed at a uniform fractions from the bottom node of each edge.
	 */
	public List<TreeNode> addAuxiliaryInternalNodes(UnrootedTree urt, Random rand, TreeNode current) {

		List<TreeNode> result = Lists.<TreeNode>newArrayList();    
//		System.out.println(current);
		List<Pair<TreeNode, TreeNode>> _rootedEdges = urt.getRootedEdges(current);  

		for (final Pair<TreeNode, TreeNode> edge : _rootedEdges) {
			if ((edge.getLeft().equals(current) || edge.getRight().equals(current))) {
				result.add(current);
			} else {
				double ratio = rand.nextDouble();
				double originalBL = (urt.getBranchLength(edge.getLeft(), edge.getRight())).doubleValue();
				double bottomBL = (ratio * originalBL);
				double top_BL = ((1.0 - ratio) * originalBL);
				urt.removeEdge(edge.getLeft(), edge.getRight());
				TreeNode dummyNode = TreeNode.nextUnlabelled();
				urt.addNode(dummyNode);
				urt.addEdge(edge.getLeft(), dummyNode, top_BL);
				urt.addEdge(dummyNode, edge.getRight(), bottomBL);
				result.add(dummyNode);
			}
		}
		return result;
	}


	
	

	public Map<Pair<TreeNode, TreeNode>, List<TreeNode>>  addMultipleAuxiliaryInternalNodes(UnrootedTree urt,  TreeNode current, int N) {

		Map<Pair<TreeNode, TreeNode>, List<TreeNode>> result = new HashMap<Pair<TreeNode, TreeNode>, List<TreeNode>>();	
		//    List<TreeNode> result = Lists.<TreeNode>newArrayList();    
//		System.out.println(current);
		List<Pair<TreeNode, TreeNode>> _rootedEdges = urt.getRootedEdges(current);  

		for (final Pair<TreeNode, TreeNode> edge : _rootedEdges) {
			result.put(edge,  Lists.<TreeNode>newArrayList());
//			if ((edge.getLeft().equals(current) || edge.getRight().equals(current))) {
//				result.get(edge).add(current);
//			} else {
				double originalBL = (urt.getBranchLength(edge.getLeft(), edge.getRight())).doubleValue();
				double BL = originalBL/N;
				TreeNode dummyNode1 =  edge.getLeft(), dummyNode2=null;	        
				for (int i=0;i<N;i++) {
					dummyNode2 = TreeNode.nextUnlabelled();
					urt.addNode(dummyNode2);    	        
					result.get(edge).add(dummyNode2);
					urt.addEdge(dummyNode1, dummyNode2, BL);
					dummyNode1 =dummyNode2;     	            		  
				} 

				urt.addEdge(dummyNode1, edge.getRight(), BL);
				urt.removeEdge(edge.getLeft(), edge.getRight());
			}
//		}
		return result;
	}






	private static double sum(Pair<Double,Double> pair)
	{
		return pair.getLeft() + pair.getRight();
	}

	private static Pair<Double,Double> neighborBranchLengths(TreeNode treeNode, UnrootedTree tree)
	{
		List<TreeNode> attachmentNeighbors = Graphs.neighborListOf(tree.getTopology(), treeNode);
		if (attachmentNeighbors.size() != 2)
			throw new RuntimeException();
		double
		b1 = tree.getBranchLength(attachmentNeighbors.get(0), treeNode), 
		b2 = tree.getBranchLength(attachmentNeighbors.get(1), treeNode);
		return Pair.of(b1, b2);
	}





}
