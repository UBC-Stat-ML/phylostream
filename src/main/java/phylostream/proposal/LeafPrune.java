package phylostream.proposal;

//import bayonet.distributions.Multinomial;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

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
import briefj.Indexer;
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
	private int[]  edgeNumbers; 
	private int edgeNumberInNeighbor; 
	private double acceptanceRate=0; 

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


	public UnrootedTree urtAfterOneLeafRemoval() {
		return urtAfterOneLeafRemoval; 
	}

	
	
public boolean  expNormalizeLikelihoods(Map<Pair<TreeNode,TreeNode>, Double> likelihoods) {		
	    double max = Double.NEGATIVE_INFINITY;
	    for(Pair<TreeNode,TreeNode> key : likelihoods.keySet()) 	    
	    	max = Math.max(max, likelihoods.get(key));
	    double sum = 0;
	    for(Pair<TreeNode,TreeNode> key : likelihoods.keySet()) 
	    {
	    	double temp = Math.exp(likelihoods.get(key)-max);
	    	likelihoods.put(key, temp);
	    	sum += temp;	    	
	    }
	    if(sum == 0) return false;
	    for(Pair<TreeNode,TreeNode> key : likelihoods.keySet())
	    	likelihoods.put(key, likelihoods.get(key)/sum);
	    return true;
	}

	
	public Map<Pair<TreeNode,TreeNode>, Double>  attachmentPointsLikelihoods(Random random, int N) {

		Map<Pair<TreeNode,TreeNode>, Double> result =  new HashMap<Pair<TreeNode, TreeNode>,Double>();			 
		UnrootedTree urt = new UnrootedTree(urtree);
		Pair<TreeNode, TreeNode> randomEdge = edgeAttachedToLeaf(random, urt);
		TreeNode removedInternalNode = randomEdge.getLeft(), removedLeaf=randomEdge.getRight();

		// disconnect the leaf from the tree
		prunedSubtree = pruneLeaf(urt, removedLeaf, removedInternalNode); 
//		System.out.println(prunedSubtree);		
		urtAfterOneLeafRemoval = urt;  
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
		List<Pair<TreeNode, TreeNode>> rootedEdges = urtAfterOneLeafRemoval.getRootedEdges(removedInternalNode);
		for (Pair<TreeNode, TreeNode> edge:rootedEdges) {			
			UnrootedTree urt1 = new UnrootedTree(urtAfterOneLeafRemoval);
			List<TreeNode> attachmentPoints = addMultipleAuxiliaryInternalNodesToOneEdge(urt1, edge, N);  
			// run the sum product on the main tree
			List<SumProduct<TreeNode>> mainTreeSumProducts = EvolutionaryModelUtils.getSumProductsFromFactorGraphs(EvolutionaryModelUtils.buildFactorGraphs(evolutionaryModel, urt1, removedInternalNode, data, false), removedInternalNode);
			double edgeLoglikelihood = Double.NEGATIVE_INFINITY;  
			for (int i = 0; i < attachmentPoints.size(); i++)
			{
				TreeNode attachmentPoint = attachmentPoints.get(i);
				List<List<UnaryFactor<TreeNode>>> currentFullUnaries = Lists.newArrayList();
				for (int j=0; j<N; j++)  currentFullUnaries.add(Lists.newArrayList());
				for (int f = 0; f < mainTreeSumProducts.size(); f++)
				{
					SumProduct<TreeNode> currentMainTreeSumProduct = mainTreeSumProducts.get(f);
					UnaryFactor<TreeNode> currentMainTreeMarginal = currentMainTreeSumProduct.computeMarginal(attachmentPoint);
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


	
public double[][] multiply(double[][] A, double[][] B) {
		double[][] C = new double[A.length][B[0].length];

		for(int i=0; i<C.length; i++){
			for(int k=0; k<A[0].length; k++){
				if(A[i][k]!=0){
					for(int j=0; j<C[0].length; j++){
						C[i][j] += A[i][k]*B[k][j];
					}
				}
			}
		}

		return C;
	}	
	
	
//public double[][] multiplyTransitionProbMat(double[][] A, int n){
//	if(n==1) return A; 
//	double[][] result = A; 
//	for(int i=0;i<n;i++)
//	result = multiply(result, A);
//	return(result);
//}


public static int log2(int N)
{
    // calculate log2 N indirectly
    // using log() method
    int result = (int)(Math.log(N) / Math.log(2));
    return result;
}


public double[][] multiplyTransitionProbMat(double[][] A, int K){
	
	if(K==1) return A; 
	double[][] result=A;  			
	for(int k=0; k<K;k++)
	{
		result = multiply(result, result);
	}	
	return(result);
}




public double totalVariation(double[] stationaryDist, double[][] A) {
	if(stationaryDist.length!=A.length) throw new RuntimeException("Length of A doesn't match the length of the stationary distribution!");
	double tv = 0; 
	for(int i=0; i<A.length; i++)
	{	
		double sum=0;
		for(int j=0; j<stationaryDist.length; j++) 
			sum = sum + 0.5*Math.abs(stationaryDist[j] - A[i][j]);
		if(sum>tv) tv=sum;
	}
	return tv;
}


public static boolean expNormalize(double[] probs) {
    // Input: log probabilities (unnormalized too)
    // Output: normalized probabilities
    // probs actually contains log probabilities; so we can add an arbitrary constant to make
    // the largest log prob 0 to prevent overflow problems
    double max = Double.NEGATIVE_INFINITY;
    for(int i = 0; i < probs.length; i++)
      max = Math.max(max, probs[i]);
    for(int i = 0; i < probs.length; i++)
      probs[i] = Math.exp(probs[i]-max);
    return normalize(probs);
  }


public static boolean normalize(double[] data) {
    double sum = 0;
    for(double x : data) sum += x;
    if(sum == 0) return false;
    for(int i = 0; i < data.length; i++) data[i] /= sum;
    return true;
  }

public boolean stationaryDist(double[] loglikelihoodsVec)
{
	return(expNormalize(loglikelihoodsVec));
}


public double[] totalVariationSequence(double[] likelihoodsVec, int n) {	
	double[] result = new double[n];
	if(stationaryDist(likelihoodsVec)) {	  //	likelihoodsVec expNormalized  
	double[][] A= transitionProb(likelihoodsVec);  	
	for(int i=0;i<n;i++) 	
		result[i]=totalVariation(likelihoodsVec,A);
		A = multiply(A, A);
	}
	return(result);
}




public double[] totalVariationSequenceNearestNeighbor(Map<Pair<TreeNode,TreeNode>, Double> likelihoods, int K, int neighborSize, double power) {
//	int K = log2(n);
	double[] result = new double[K];
	double[] likelihoodsVec=loglikelihoodsVec(likelihoods);
	if(stationaryDist(likelihoodsVec)) {		   //	likelihoodsVec expNormalized
	double[][] A= transitionProb(urtAfterOneLeafRemoval, likelihoods, neighborSize, power);  // Note likelihoods are in log scale   
	for(int i=0;i<K;i++) 	
	result[i]=totalVariation(likelihoodsVec,multiplyTransitionProbMat(A, (i+1)));
	}
	return(result);
}


public double[] totalVariationSequenceNearestNeighbor(Map<Pair<TreeNode,TreeNode>, Double> likelihoods, int K, int neighborSize, double power, double mixtureProportion) {
//	int K = log2(n);
	double[] result = new double[K];
	double[] likelihoodsVec=loglikelihoodsVec(likelihoods);
	if(stationaryDist(likelihoodsVec)) {		   //	likelihoodsVec expNormalized
	double[][] A= transitionProbMixtureMove(urtAfterOneLeafRemoval, likelihoods, neighborSize, power, mixtureProportion);  // Note likelihoods are in log scale   
	for(int i=0;i<K;i++) 	
	result[i]=totalVariation(likelihoodsVec,multiplyTransitionProbMat(A, (i+1)));
	}
	return(result);
}




public double[] loglikelihoodsVec(Map<Pair<TreeNode,TreeNode>, Double> likelihoods) {
	int nEdge = likelihoods.keySet().size();
	int k= 0;
	double[] loglikelihoodsVec = new double[nEdge];
	for(Pair<TreeNode, TreeNode> edge:likelihoods.keySet()) {
		loglikelihoodsVec[k] = likelihoods.get(edge); 
		k=k+1;
	}
	return(loglikelihoodsVec);
}
	
	


public double[][] neighborDistance(UnrootedTree urt, Map<Pair<TreeNode,TreeNode>, Double> likelihoods, int neighborSize) {	
	Set<Pair<TreeNode,TreeNode>> edgeSet = likelihoods.keySet();
	Indexer<Pair<TreeNode,TreeNode>> index = new Indexer<Pair<TreeNode,TreeNode>>(edgeSet);		
	int nEdge=index.size();	
//	System.out.println("Edge number is "+ nEdge);
	double distance[][] = new double[nEdge][nEdge];
	for(int i = 0; i< nEdge; i++)
		for(int j=0; j<nEdge; j++)
			distance[i][j] = Double.MAX_VALUE; 
	
	for(int i = 0; i< nEdge; i++)
	{		
		Pair<TreeNode, TreeNode> edge =  index.i2o(i); 
		Map<Pair<TreeNode,TreeNode>, Integer>  neighbors = neighborhoods(urt, edge, neighborSize);
		for(int j=0; j< nEdge; j++)
		{
			if(i!=j) {
			List<Pair<TreeNode, TreeNode>> unorderedEdge = Lists.newArrayList();
			Pair<TreeNode, TreeNode> edge2 = index.i2o(j); 
			unorderedEdge.add(edge2);
			unorderedEdge.add(Pair.of(edge2.getRight(), edge2.getLeft()));
			if(!Collections.disjoint(unorderedEdge, neighbors.keySet())) { 	// edge2 is in the neighbor of edge
				
				if(neighbors.keySet().contains(edge2))
					distance[i][j] = neighbors.get(edge2); 
				else
					distance[i][j] = neighbors.get(unorderedEdge.get(1));
			}
			}
		}
	}
	return(distance); 
}



// likelihoods are in log scale
public double[][] neighborLikelihood(UnrootedTree urt, Map<Pair<TreeNode,TreeNode>, Double> likelihoods, int neighborSize) {	
	Set<Pair<TreeNode,TreeNode>> edgeSet = likelihoods.keySet();
	Indexer<Pair<TreeNode,TreeNode>> index = new Indexer<Pair<TreeNode,TreeNode>>(edgeSet);		
	int nEdge=index.size();	
	System.out.println("Edge number is "+ nEdge);
	double neighborLikelihood[][] = new double[nEdge][nEdge];
		
	for(int i = 0; i< nEdge; i++)
	{		
		Pair<TreeNode, TreeNode> edge =  index.i2o(i); 
		Map<Pair<TreeNode,TreeNode>, Integer>  neighbors = neighborhoods(urt, edge, neighborSize);		
		
		for(int j=0; j< nEdge; j++)
		{
			if(i!=j) {
			List<Pair<TreeNode, TreeNode>> unorderedEdge = Lists.newArrayList();
			Pair<TreeNode, TreeNode> edge2 = index.i2o(j); 
			unorderedEdge.add(edge2);
			unorderedEdge.add(Pair.of(edge2.getRight(), edge2.getLeft()));
			if(!Collections.disjoint(unorderedEdge, neighbors.keySet())) { 	// edge2 is in the neighbor of edge				
				if(likelihoods.keySet().contains(edge2))
					neighborLikelihood[i][j] =  likelihoods.get(edge2); 
				else {
					neighborLikelihood[i][j] = likelihoods.get(unorderedEdge.get(1));
				}
			}
		else
		{
			neighborLikelihood[i][j] = Double.NEGATIVE_INFINITY; 
		}
	}
			else
				neighborLikelihood[i][j] = likelihoods.get(edge);
		}
	}
	return(neighborLikelihood); 
}




// neighborLikelihood in log scale
public double[][] neighborWeight(double[][] neighborLikelihood, double power) {
	int sz = neighborLikelihood.length;

	double neighborWeight[][] = new double[sz][sz];	
	for(int i=0; i<sz; i++)
		for(int j=0; j<sz; j++)
			if(neighborLikelihood[i][j]== Double.NEGATIVE_INFINITY) 
			neighborWeight[i][j] = Double.NEGATIVE_INFINITY;
	
	for(int i=0; i<sz; i++)
	{
		double max = Double.NEGATIVE_INFINITY;
		for(int j=0; j<sz; j++)
			if(neighborLikelihood[i][j]!= Double.NEGATIVE_INFINITY)
				max = Math.max(max, power*neighborLikelihood[i][j]);
	    double sum = 0;
	    for(int j=0; j<sz; j++) 
	    {
	    	if(neighborLikelihood[i][j]!= Double.NEGATIVE_INFINITY)
	    	{
	    	double temp = Math.exp(power*neighborLikelihood[i][j]-max);
	    	neighborWeight[i][j] = temp; 	    	
	    	sum += temp;
	    	}
	    }
	    if(sum == 0) throw new RuntimeException();
	    for(int j=0; j<sz; j++)
	    	if(neighborLikelihood[i][j]!= Double.NEGATIVE_INFINITY)
	    	  neighborWeight[i][j] = neighborWeight[i][j]/sum;		
	}		
	return(neighborWeight);
}



public double[][] transitionProb(double[] likelihoodsVec) {		
	int nEdge=likelihoodsVec.length;	
	double transitionProb[][] = new double[nEdge][nEdge];
	for(int i=0;i<nEdge;i++) {
		double sum = 1.0/nEdge;
		for(int j=0;j<nEdge;j++)
			if(i!=j) {
			double mhRatio = Math.min(likelihoodsVec[j]/likelihoodsVec[i], 1);	
			transitionProb[i][j] = (1.0/nEdge)*mhRatio; // assuming uniform prior
			sum +=  (1.0/nEdge)*(1-mhRatio);  
			}
		transitionProb[i][i] = sum;    
	}
	
//	for(int i=0;i<nEdge;i++) { 
//		double sum=0;
//		for(int j=0;j<nEdge;j++) { 
//			System.out.print(transitionProb[i][j]+" ");
//			sum = sum+transitionProb[i][j];
//		}
//		System.out.println("add up to "+sum);		
//	}
	return(transitionProb); 
}



//likelihoods are in log scale
public double[][] transitionProbMixtureMove(UnrootedTree urt, Map<Pair<TreeNode,TreeNode>, Double> likelihoods, int neighborSize, double power, double mixtureProportion) {		
	double[][] neighborLikelihood = neighborLikelihood(urt, likelihoods, neighborSize);
	int nEdge = neighborLikelihood.length;	
	double[][]  weights = neighborWeight(neighborLikelihood, power);
	
	double transitionProb[][] = new double[nEdge][nEdge];
	int edgeNumber = 0; 
	edgeNumbers = new int[nEdge];
	
	acceptanceRate = 0; 
	for(int i = 0; i< nEdge; i++)
	{
		int edges = 0;
		double sum = mixtureProportion*1.0/nEdge+ (1-mixtureProportion)*weights[i][i];				
		for(int j=0; j< nEdge; j++)
		{			
			if(i!=j) {
//			if(neighborLikelihood[i][j]!= Double.NEGATIVE_INFINITY) {	
			if(weights[i][j]!= Double.NEGATIVE_INFINITY) {
				edges += 1;  					
			double likelihoodRatio = Math.exp(neighborLikelihood[i][j] - neighborLikelihood[i][i]); 
			double mhRatioUniform = Math.min(likelihoodRatio, 1); 
			double mhRatioLocalInformative = Math.min(likelihoodRatio*weights[j][i]/weights[i][j], 1); 			
			transitionProb[i][j] = mixtureProportion*(1.0/nEdge)*mhRatioUniform + (1-mixtureProportion)*weights[i][j] *mhRatioLocalInformative;
			acceptanceRate += transitionProb[i][j]; 
			sum +=  mixtureProportion*(1.0/nEdge)*(1-mhRatioUniform) + (1-mixtureProportion)*weights[i][j]*(1-mhRatioLocalInformative);  
			}
			else 
			{
				double likelihoodRatio = Math.exp(neighborLikelihood[i][j] - neighborLikelihood[i][i]); 
				double mhRatioUniform = Math.min(likelihoodRatio, 1); 
				transitionProb[i][j] = mixtureProportion*(1.0/nEdge)*mhRatioUniform;
				acceptanceRate += transitionProb[i][j]; 
				sum +=  mixtureProportion*(1.0/nEdge)*(1-mhRatioUniform);  				
			}
		}
		}
			transitionProb[i][i] = sum;
			edgeNumbers[i] = edges+1;
			edgeNumber += edgeNumbers[i]; 
			System.out.println(edges); 
	}
	
//	for(int i=0;i<nEdge;i++) { 
//	double sum=0;
//	for(int j=0;j<nEdge;j++) { 
//		System.out.print(transitionProb[i][j]+" ");
//		sum = sum+transitionProb[i][j];
//	}
//	System.out.println("add up to "+sum);		
//}
	edgeNumberInNeighbor = edgeNumber/nEdge;
	acceptanceRate /= nEdge;   //Uniformly pick one edge. 
	return(transitionProb); 
}




//likelihoods are in log scale
public double[][] transitionProb(UnrootedTree urt, Map<Pair<TreeNode,TreeNode>, Double> likelihoods, int neighborSize, double power) {		
	double[][] neighborLikelihood = neighborLikelihood(urt, likelihoods, neighborSize);
	int nEdge = neighborLikelihood.length;	
	double[][]  weights = neighborWeight(neighborLikelihood, power);
	
	double transitionProb[][] = new double[nEdge][nEdge];
	int edgeNumber = 0; 
	edgeNumbers = new int[nEdge];
	
	acceptanceRate = 0; 
	for(int i = 0; i< nEdge; i++)
	{
		int edges = 0;
		double sum = weights[i][i];				
		for(int j=0; j< nEdge; j++)
		{			
			if(i!=j) {
//			if(neighborLikelihood[i][j]!= Double.NEGATIVE_INFINITY) {	
			if(weights[i][j]!=  Double.NEGATIVE_INFINITY) {
				edges += 1;  	
				
			double mhRatio = Math.min(Math.exp(neighborLikelihood[i][j] - neighborLikelihood[i][i])*weights[j][i]/weights[i][j], 1);
			
			if(mhRatio==0.0) System.out.println("WARNING: MH is 0!!");
			transitionProb[i][j] = weights[i][j] *mhRatio;
			acceptanceRate += transitionProb[i][j]; 
			sum +=  weights[i][j]*(1-mhRatio);  
			}
		}
		}
			transitionProb[i][i] = sum;
			edgeNumbers[i] = edges+1;
			edgeNumber += edgeNumbers[i]; 
			System.out.println(edges); 
	}
	edgeNumberInNeighbor = edgeNumber/nEdge;
	acceptanceRate /= nEdge;   //Uniformly pick one edge. 
	return(transitionProb); 
}






public int getAverageEdgeNumberInNeighbor()
{
	return (edgeNumberInNeighbor); 
}


public int[] getEdgeNumbers()
{
	return(edgeNumbers); 
}

public double getAcceptanceRate()
{
	return(acceptanceRate);	
}



public List<Pair<TreeNode,TreeNode>>  neighborhoods(UnrootedTree urt, Pair<TreeNode, TreeNode> edge) {
	TreeNode left = edge.getLeft();
	TreeNode right = edge.getRight();	
	List<TreeNode>  leftNodeNeighbors  = Graphs.neighborListOf(urt.getTopology(), left);
	List<TreeNode>  rightNodeNeighbors = Graphs.neighborListOf(urt.getTopology(), right);
	leftNodeNeighbors.remove(right);
	rightNodeNeighbors.remove(left);
	List<Pair<TreeNode,TreeNode>> result = Lists.newArrayList(); 	
	for (TreeNode node : leftNodeNeighbors)
		result.add(Pair.of(left,node));   		
	for (TreeNode node : rightNodeNeighbors)
		result.add(Pair.of(right,node));  
    return result;
}




public Map<Pair<TreeNode,TreeNode>, Integer>  neighborhoods(UnrootedTree urt, Pair<TreeNode, TreeNode> edge, int radius) {
	Map<Pair<TreeNode, TreeNode>, Integer> result = new HashMap<Pair<TreeNode, TreeNode>, Integer>();
	result.putAll(descendant(urt, edge, radius)); 
	result.putAll(descendant(urt, Pair.of(edge.getRight(), edge.getLeft()), radius)); 	
    return result;
}



public Map<Pair<TreeNode,TreeNode>, Integer>  descendant(UnrootedTree urt, Pair<TreeNode, TreeNode> edge, int radius)
{
	Map<Pair<TreeNode,TreeNode>, Integer> result = new HashMap<Pair<TreeNode,TreeNode>, Integer>();  		
	List<Pair<TreeNode,TreeNode>> edges = Lists.newArrayList();
	edges.add(edge);
	for(int i=0; i<radius; i++)
	{
		List<Pair<TreeNode,TreeNode>> children = Lists.newArrayList();
		for(Pair<TreeNode,TreeNode> edge0: edges)
			children.addAll(children(urt, edge0));		
		for(int j=0; j<children.size(); j++)
			result.put(children.get(j), i); 
		edges.removeAll(edges); 
		edges.addAll(children);			
	}
   return result;
}




public List<Pair<TreeNode,TreeNode>>  children(UnrootedTree urt, Pair<TreeNode, TreeNode> edge)
{
	TreeNode left = edge.getLeft();
	TreeNode right = edge.getRight();	
	List<TreeNode>  rightNodeNeighbors = Graphs.neighborListOf(urt.getTopology(), right);
	rightNodeNeighbors.remove(left);	
	List<Pair<TreeNode,TreeNode>> result = Lists.newArrayList(); 	
	for (TreeNode node : rightNodeNeighbors)
		result.add(Pair.of(right,node));  
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





	public List<TreeNode>  addMultipleAuxiliaryInternalNodesToOneEdge(UnrootedTree urt, Pair<TreeNode, TreeNode> edge,  int N) {
		List<TreeNode> result = Lists.<TreeNode>newArrayList();	
		double originalBL = (urt.getBranchLength(edge.getLeft(), edge.getRight())).doubleValue();
		double BL = originalBL/N;
		TreeNode dummyNode1 =  edge.getLeft(), dummyNode2=null;	        
		for (int i=0;i<N;i++) {
			dummyNode2 = TreeNode.nextUnlabelled();
			urt.addNode(dummyNode2);    	        
			result.add(dummyNode2);
			urt.addEdge(dummyNode1, dummyNode2, BL);
			dummyNode1 =dummyNode2;     	            		  
		} 
		urt.addEdge(dummyNode1, edge.getRight(), BL);
		urt.removeEdge(edge.getLeft(), edge.getRight());
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


}
