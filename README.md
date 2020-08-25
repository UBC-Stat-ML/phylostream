Summary [![Build Status](https://travis-ci.org/UBC-Stat-ML/phylostream.png?branch=master)](https://travis-ci.org/UBC-Stat-ML/phylostream) 
-------

Development setup (eclipse)
-----------------

Requires Open SDK 8.

- Type ``./setup-eclipse.sh`` from the root of the repository
- From eclipse:
  - ``Import`` in ``File`` menu
  - ``Import existing projects into workspace``
  - Select the root of this repository


Development setup (CLI)
-----------------

- Type ``./setup-cli.sh`` from the root of the repository


Running phylostream
-----------------

To run ``PhyloModel.bl``, from cli:

```
phylostream --engine phylostream.blang.PhyloStreamSMC --engine.nThreads Max    --engine.nParticles 100 --engine.temperatureSchedule.threshold 0.8 --model.data Synthetic --model.data.nLeaves 100 --model.data.nSites 100
```

From eclipse: 

- Open ``PhyloModel.bl``
- Right-click on the editor showing the model and select ``Open generated file``
- Right-click on the editor showing the generated source and select ``Run as Java application``
- You can now edit the command line arguments in ``Run > Run Configurations... > Arguments``