# Voronoi_xAOD
This is a tool for calculating and applying cluster-based pileup mitigation and subtraction algorithms. It also includes tools for writing out certain quantities to an output tree for performance testing.

## Dependencies
This package makes use of [UChicago](https://github.com/UCATLAS)'s [xAODAnaHelpers](https://github.com/UCATLAS/xAODAnaHelpers) package and [Giordon Stark](https://github.com/kratsg)'s [xAODJetReclustering](https://github.com/kratsg/xAODJetReclustering) package.

## Installing
To install,
```bash
mkdir myRootCore && cd $_
setupATLAS
rcSetup Base,2.4.33
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/AvivCukierman/VoronoiWeightTool
git clone https://github.com/AvivCukierman/Voronoi_xAOD
mv Voronoi_xAOD MyAnalysis
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJet/tags
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJetContrib/tags
rc find_packages
rc compile
```

## Running
Two sample run scripts (`run_lsf_test.sh`, for running on a file stored on the grid) and (`run_test.sh`, for running on a file stored locally) are included in the scripts folder. The `--driver` option allows you to run on batch, by specifying `--driver lsf`. I'll document this better later but just ask me for help if you need it for now.

## How it works
### MyxAODAnalysis
This analysis applied the [VoronoiWeightTool](https://github.com/AvivCukierman/VoronoiWeightTool) tool for applying Voronoi weights to input clusters. It then clusters the resulting clusters into jets using xAODJetReclustering.
### JetMatching
This tool calculates two things. First, it matches jets in a jet container to jets in a truth jet container. Then, it attaches a decoration to each truth jet as the minimum distance to any other truth jet (pT > 5 GeV) and to each reco jet as the minimum distance to any other reco jet (pT > 2 GeV). There is also a "isPU" label attached to reco jets that are more than 0.6 from any truth jet above 5 GeV. These metrics are for performance studies, as we look at reconstructed jets matched to isolated truth jets for performance studies, and also calculate the fake rate.
### WriteTree
This tool writes out a tree containing some quantities necessary for performance testing, such as reconstructed pT, matched truth jet pT, etc. Again, this tool is for performance studies.
