# README #

This repository contain the following folders used on RdelaCruz et al (2017) paper: Stochastic, hybrid, hybrid bis, coarse grained and mean field. 

We introduce a modification on these programs. Now, we study the behavour of two population (host and invader). Host population is in equilibrium and invader has a faster cell-cycle and few cells. 

### Stochastic ###

* Stochastic method with age structure
* Version 2017
* 1-D and 1-Pop


### Mean Field (mf) ###

* Mean field model with age structure
* Version 2017
* 1-D and 1-Pop

### Hybrid and Hybridbis ###

* Determinist population if it is over a thershold and if the number of cells is above this thershold the populatios is stochastic
* Version 2017
* 1-D and 1-Pop

### Coarse grained ###

* Faster code which aproximate mean field. Cg integrates age term in mean field.
* Version 2017
* 1-D and 1-Pop

## Two population ##

### Coarse grained ###

* Include two different population. We read six different initial files:
	* INITIAL_DATUM_host
	* INITIAL_DATUM_invader
	* Initial Oxigen distribution
	* Host parameter
	* Invader parameter
	* General parameter (dt, dh, ...)

