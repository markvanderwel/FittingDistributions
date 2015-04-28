# FittingDistributions
R code for inverse CAIN analysis across eastern US

Created 4/3/15

File descriptions:

##### GitHub files
* .gitignore
* README.md

##### R Project file
* Fitting distributions.Rproj 

##### R Scripts
* cain gof - use this to create figures used in ESA 2014 talk (incl animation).R 	
  * generates figures from esa talk
* compare modelled diam distributions to fia.R 	
  * generates figure comparing diameter distributions
* estimating harvest rates 2.R 
  * estimates (and plots) harvest rates for each grid cell/plot class
* Graphs.R
  * generates graphs
* kriging testing 2 - use this to generate trend surface or krigged priors from ncdf results.R 
  * creates (and plots) smooth trend surface for rate parameters, which are used in hierarchical model 
* plot rates against climate.R 
  * compares prior/posterior distributions for rate paramters against temp, precip
* plotcompclustering.R 
  * clusters plots within grid cells based on species composition
* Slopes comparison within and among grid cells.R
  * compares slopes of demographic rates vs. predicted BA among and within grid cells

##### Data files
* FittingOutput.nc 
  * main model output file
* SimOutWithClimData.nc
  * output from original model fit only to demographic data
* HarvestPlotsList.csv 
  * estimated harvest rate for each fia plot
* plotclasses.csv 
  * assignments of fia plots to plot classes within grid cells
