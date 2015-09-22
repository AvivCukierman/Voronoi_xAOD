# Voronoi_xAOD
This is a tool for calculating and applying cluster-based pileup mitigation and subtraction algorithms. It also includes tools for writing out certain quantities to an output tree for performance testing.

## Dependencies
This package makes use of [UChicago](https://github.com/UCATLAS)'s [xAODAnaHelpers](https://github.com/UCATLAS/xAODAnaHelpers) package and [Giordon Stark](https://github.com/kratsg)'s [xAODJetReclustering](https://github.com/kratsg/xAODJetReclustering) package.

## Installing
To install,
```bash
mkdir myRootCore && cd $_
rcSetup Base,2.3.23
git clone https://github.com/kratsg/xAODJetReclustering
git clone https://github.com/UCATLAS/xAODAnaHelpers
rc find_packages
rc compile
```

## Running
Two sample run scripts (`run_lsf_test.sh`, for running on a file stored on the grid) and (`run_test.sh`, for running on a file stored locally) are included in the scripts folder. The `--driver` option allows you to run on batch, by specifying `--driver lsf`. I'll document this better later but just ask me for help if you need it for now.

## How it works
### VoronoiWeights
This tool calculates the Voronoi subtraction algorithm (in `MakeVoronoiClusters`) and applies a weight to clusters as a `SG::AuxElement::Decorator< float >`. A different algorithm would simply calculate a different weight than `MakeVoronoiClusters` but the weight could be applied in the same way. The assumption is that if the weight is 0, then that cluster should be ignored in later jet clustering. If the weight is not 0, then that is the value used as the new pT in jet clustering.
### VoronoiJets
This tool clusters jets using the weights applied in `VoronoiWeights`. The weights are used appropriately as described in the previous section. For now the output jet container is called "AntiKt4VoronoiJets" but I guess I can make that configurable in the run script in a future revision.
### JetMatching
This tool calculates two things. First, it matches jets in a jet container to jets in a truth jet container. Then, it attaches a decoration to each truth jet as the minimum distance to any other truth jet (pT > 5 GeV). Both of these are only for testing, as we look reconstructed jets matched to isolated truth jets for performance studies.
### WriteTree
This tool writes out a tree containing some quantities necessary for performance testing, such as reconstructed pT, matched truth jet pT, etc. Again, this tool is only for testing.
