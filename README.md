# FittingDistributions
R code for inverse CAIN analysis across eastern US

Created 4/3/15

File descriptions:

##GITHUB FILES
* .gitignore
* README.md

##R PROJECT FILE
* Fitting distributions.Rproj 

##R SCRIPTS
* cain gof - use this to create figures used in ESA 2014 talk (incl animation).R 	
  * generates figures from esa talk
* compare modelled diam distributions to fia.R 	
  * generates figure comparing diameter distributions
* estimating harvest rates 2.R 
  * estimates (and plots) harvest rates for each grid cell/plot class
* kriging testing 2 - use this to generate trend surface or krigged priors from ncdf results.R 
  * creates (and plots) smooth trend surface for rate parameters, which are used in hierarchical model 
* plot rates against climate.R 
  * compares prior/posterior distributions for rate paramters against temp, precip
* plotcompclustering.R 
  * clusters plots within grid cells based on species composition

##DATA FILES
* FittingOutput.nc 
  * main model output file
* HarvestPlotsList.csv 
  * estimated harvest rate for each fia plot
* plotclasses.csv 
  * assignments of fia plots to plot classes within grid cells
