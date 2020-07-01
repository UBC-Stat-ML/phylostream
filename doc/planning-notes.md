## Logistics

- Project name
- Git group

## Possible initial simplifying assumptions

- Fixed evolutionary parameters
- No rejuvenation
- No real data
- Adding a single taxon on a fixed $n$ leaves tree
- 2 main moves, simplify to one?

Question: how many adaptive anneal steps to add one taxon? How does this number grow in $n$?


## Useful code to borrow


### Bayonet's ``SumProduct``

- Terminology: [``FactorGraph``](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/main/java/bayonet/marginal/FactorGraph.java#L16): a factor graph. The one we will be using is [``DiscreteFactorGraph``](https://github.com/UBC-Stat-ML/bayonet/blob/master/src/main/java/bayonet/marginal/DiscreteFactorGraph.java) which is a tree shaped factor graph with discrete state space and embed in a plate (the plate will be ranging over sites)


- [``SumProduct``](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/main/java/bayonet/marginal/algo/SumProduct.java): Takes a ``FactorGraph``, and provides various exact inference algorithms on that factor graph. 

- Why this is useful to borrow: [implements the necessary scalings black magic needed to get reasonable speed](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/main/java/bayonet/marginal/DiscreteFactorGraph.java#L679), [takes care of the annoying subcases of arbitrary aritcy](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/main/java/bayonet/marginal/DiscreteFactorGraph.java#L408), while presenting a clean level of abstraction to define our more advance quick updates for sum product algorithms, [``SumProduct``](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/main/java/bayonet/marginal/algo/SumProduct.java)

- [Pedagogical example of how to setup a factor graph](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/test/java/bayonet/factors/SumProductTests.java#L32) and [how to use it](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/test/java/bayonet/factors/SumProductTests.java#L20). Note: other code already exists that builds factor graphs for phylogenetic tree, see conifer below.

- Initialization of the real-time algorithm: compute the incoming message to the left and right of an edge: [SumProduct.computeSubtreeMarginal](https://github.com/UBC-Stat-ML/bayonet/blob/32a40492f66574fec8906001bcdc8b20c9af01f6/src/main/java/bayonet/marginal/algo/SumProduct.java#L91)

- A few more methods will need to be added to ``SumProduct`` to support parsimonious updates when local updates are made to the factor graph. This is delicate but a simple way to test is to have a test based on an inefficient subclass which check each incremental update with a fresh recompute. This test will be imperative. 

- Currently unclear: FactorGraph was designed with immutability in mind. Not sure yet what would be the right way to do expose APIs in SumProduct and FactorGraph for doing small updates on the factor graphs. API should be on SumProduct while propagating changes to FactorGraph. Would be easier if messages would have been stored with an object strategy? Yes, perhaps, by using a trick where message structure is based on immutable pointer to its parent. Might want to bypass SumProduct and use only DiscreteFactorGraph directly.

### Conifer's ``EvolutionaryModel``

- Terminology: 
    - [``UnrootedTree``](https://github.com/UBC-Stat-ML/conifer/blob/master/src/main/java/conifer/UnrootedTree.xtend) = discrete tree + branch lengths
    - [``EvolutionaryModel``](https://github.com/UBC-Stat-ML/conifer/blob/master/src/main/java/conifer/EvolutionaryModel.java): a recipe to produce a ``FactorGraph`` from an ``UnrootedTree``, observations and potentially other model parameters. 
    - The one we would be using is the standard multi-category mixture, implemented in [``MultiCategorySubstitutionModel``](https://github.com/UBC-Stat-ML/conifer/blob/master/src/main/java/conifer/models/MultiCategorySubstitutionModel.java)

- To create a fresh factor graph, use  [EvolutionaryModelUtils.buildFactorGraphs](https://github.com/UBC-Stat-ML/conifer/blob/99036b150865371064280b3cba05756c0e0b56b4/src/main/java/conifer/models/EvolutionaryModelUtils.java#L31). 

- We will need to pluck a leaf and move it around. This should be possible to achieve by specializing the code in the [Subtree Prune Regraft (SPR) move](https://github.com/UBC-Stat-ML/conifer/blob/master/src/main/java/conifer/moves/SPRMove.java)


### SMC code

- 2 versions available
    - Blang's SMC engine
    - Annealed SMC repo

- Which to use?


## Possible next step

New class: RealTimeSumProduct?

contain:

- factor graph
- latestTip

queries:

- computeLikelihood

edit operations:

- rebuild(EvolutionaryModel model): build factor graph, clear the caches
- addTip(newLatestTip)
- editEdge
- move(int i) i is in 1,2,3,4


Debug strategy: test which computes everything twice, once with cache and once from scratch and check they match