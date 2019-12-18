# On the Trail of the Most Massive Galaxy Clusters in the Universe

This is all of the source code for the papers associated with the Planck cluster follow up project. Papers include:

- [High Confidence Optical Confirmation Among The High Signal-to-Noise Planck Cluster Candidates](https://arxiv.org/abs/1809.06378)

Here is the basic structure of the repo. Also please have a look at the git branches for features not immediately apparent in the master branch.


- analysis_op -- all of the analysis associated with the optical cluster search
- analysis_ir -- all of the analysis associated with the IR cluster search
- catalogs -- original PSZ catalogs
- CLUSTERpipe -- pipeline for cluster finding
- data -- all of the original imaging associated with the project
- legacy -- old code I didn't want to throw away
- MOSAICpipe -- the main pipeline to process all of the MOSAIC1/2 and NEWFIRM data
- observing -- scripts to help with the observing
- papers -- papers from other people who have done similar work
- plots -- script to make (some of) the plots for the paper
- results -- results from everyone's cluster finding
- scripts -- misc. utility scripts, a catchall.
- sims -- simulations to understand the depth of the imaging
- snippets -- snippets of code that I found useful.

## Authors

- Steven Boada (Rutgers University) -- Principal Maintainer
- Jack Hughes (@jphastro, Rutgers University)
- Peter Doze (@pddoze, Rutgers University)
- Felipe Menanteau (@menanteau, University of Illinois)

## License

Copyright 2017 Steven Boada (Rutgers University).
**All rights reserved.**

## Citation
If you would like to reference our [paper](https://ui.adsabs.harvard.edu/abs/2019ApJ...871..188B),
please use the following citation, produced by
[NASA ADS](https://ui.adsabs.harvard.edu/abs/2019ApJ...871..188B):
```
@ARTICLE{2019ApJ...871..188B,
       author = {{Boada}, Steven and {Hughes}, John P. and {Menanteau}, Felipe and
         {Doze}, Peter and {Barrientos}, L. Felipe and {Infante}, L.},
        title = "{High Confidence Optical Confirmations among the High Signal-to-noise Planck Cluster Candidates}",
      journal = {\apj},
     keywords = {cosmology: observations, galaxies: clusters: general},
         year = "2019",
        month = "Feb",
       volume = {871},
       number = {2},
          eid = {188},
        pages = {188},
          doi = {10.3847/1538-4357/aaf3a0},
       adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJ...871..188B},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```

*Any use of the content of this project or repository requires citation and acknowledgment.*
