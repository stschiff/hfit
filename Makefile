build/hfit: hfit.d brent.d powell.d modelScores.d popGenFunc.d mcmc.d
	#dmd -odtmp -ofbuild/hfit -g -gc -gs -L-L/opt/local/lib -L-lgsl -L-lgslcblas -L--export-dynamic hfit.d brent.d powell.d modelScores.d popGenFunc.d
	dmd -odtmp -ofbuild/hfit -g -gc -gs -L-lgsl -L-lgslcblas hfit.d brent.d powell.d modelScores.d popGenFunc.d mcmc.d
