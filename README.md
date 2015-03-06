# FittingDistributions
R code for inverse CAIN analysis across eastern US

Created 4/3/15

File descriptions:

.gitignore
README.md

Fitting distributions.Rproj 

cain gof - use this to create figures used in ESA 2014 talk (incl animation).R 	GENERATES FIGURES FROM ESA TALK
compare modelled diam distributions to fia.R 	GENERATES FIGURE COMPARING DIAMETER DISTRIBUTIONS
estimating harvest rates 2.R ESTIMATES (AND PLOTS) HARVEST RATES FOR EACH GRID CELL/PLOT CLASS
kriging testing 2 - use this to generate trend surface or krigged priors from ncdf results.R CREATES (AND PLOTS) SMOOTH TREND SURFACE FOR RATE PARAMETERS, WHICH ARE USED TO AS A HIERARCHICAL PARAMETERS
plot rates against climate.R COMPARES PRIOR/POSTERIOR DISTRUBTIONS FOR RATE PARAMETERS AGAINST TEMP, PRECIP
plotcompclustering.R CLUSTERS PLOTS WITHIN GRID CELLS BASED ON SPECIES COMPOSITION

FittingOutput.nc MAIN MODEL OUTPUT FILE
HarvestPlotsList.csv ESTIMATED HARVEST RATE FOR EACH FIA PLOT
plotclasses.csv ASSIGNMENTS OF FIA PLOTS TO PLOT CLASSES WITHIN GRID CELLS
