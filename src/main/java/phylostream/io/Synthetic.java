package phylostream.io;


import bayonet.distributions.Random;
import blang.core.RealConstant;
import blang.core.RealDistribution;
import blang.distributions.Exponential;
import blang.inits.Arg;
import blang.inits.DefaultValue;
import conifer.EvolutionaryModel;
import conifer.SequenceAlignment;
import conifer.TopologyUtils;
import conifer.TreeNode;
import conifer.UnrootedTree;
import conifer.ctmc.CTMCParameters;
import conifer.ctmc.NoisyEmissionModel;
import conifer.ctmc.RateMatrixToEmissionModel;
import conifer.ctmc.SimpleRateMatrix;
import conifer.factors.NonClockTreePriorUtils;
import conifer.io.PhylogeneticObservationFactory;
import conifer.io.TreeObservations;
import conifer.models.DiscreteGammaMixture;
import conifer.models.MultiCategorySubstitutionModel;

/**
 * Generates trees and datasets. 
 * 
 * See TestSynthetic for usage.
 */
public class Synthetic implements Dataset {
    
  @Arg @DefaultValue("1")
  public Random rand = new Random(1);
  
  @Arg              @DefaultValue("0.01")
  public double branchMeanLength = 0.01;
  
  @Arg  @DefaultValue("10")
  public int nLeaves = 10;
  
  @Arg                      @DefaultValue("0.01")
  public double invariantSiteProbability = 0.01;
  
  @Arg            @DefaultValue("1")
  public double shapeParameter = 1;
  
  @Arg              @DefaultValue("4")
  public int nPositiveCategories = 4;
  
  @Arg              @DefaultValue("0.005")
  public double errorProbability = 0.005;
  
  @Arg         @DefaultValue("/phylostream/COVID_GTR.txt")
  public String rateMatrix = "/phylostream/COVID_GTR.txt";
  
  @Arg @DefaultValue("1000")
  public int nSites = 1000;
  
  public Realization next() {
    RealDistribution branchDistribution = Exponential.distribution(new RealConstant(1.0 / branchMeanLength)); 
    UnrootedTree tree = NonClockTreePriorUtils.sample(rand, branchDistribution, TopologyUtils.syntheticTaxaList(nLeaves));
    TreeNode arbitraryRoot = arbitraryRoot(tree);
    
    double [][] loadedRateMatrix = SimpleRateMatrix.fromResource(rateMatrix).getRateMatrix();
    int [] latent2observed = new int[] {0, 1, 2, 3};
    RateMatrixToEmissionModel emissionModel = errorProbability == 0.0 ? null : new NoisyEmissionModel(latent2observed, 4, new RealConstant(errorProbability));
    CTMCParameters ctmc = new SimpleRateMatrix(loadedRateMatrix, emissionModel);
    DiscreteGammaMixture rateMatrixMixture = new DiscreteGammaMixture(new RealConstant(invariantSiteProbability), new RealConstant(shapeParameter), ctmc, nPositiveCategories);
    MultiCategorySubstitutionModel<DiscreteGammaMixture> model = new MultiCategorySubstitutionModel<DiscreteGammaMixture>(rateMatrixMixture, nSites);
    
    SequenceAlignment data = new SequenceAlignment(PhylogeneticObservationFactory.nucleotidesFactory(), nSites);
    
    return new Realization(arbitraryRoot, tree, model, data);
  }
  
  /**
   * Test in testAddTip requires the root not to be a leaf, so try to return that if possible. 
   */
  public static TreeNode arbitraryRoot(UnrootedTree tree) {
    if (tree.leaves().size() <= 2) return TopologyUtils.arbitraryNode(tree);
    for (TreeNode node : tree.getTopology().vertexSet())
      if (node.isLabelled())
        return node;
    throw new RuntimeException();
  }
  
  public class Realization {
    public final TreeNode trueRoot;
    public final UnrootedTree trueTree;
    public final EvolutionaryModel trueModel;
    final TreeObservations data;
    
    Realization(TreeNode trueRoot, UnrootedTree trueTree, EvolutionaryModel trueModel, TreeObservations data) {
      super();
      this.trueRoot = trueRoot;
      this.trueTree = trueTree;
      this.trueModel = trueModel;
      this.data = data;
    }

    public TreeObservations nextDataset() {
      trueModel.generateObservationsInPlace(rand, data, trueTree, trueRoot);
      return data;
    }
  }
}
