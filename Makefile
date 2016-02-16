build/hfit: hfit.d brent.d powell.d modelScores.d popGenFunc.d mcmc.d
	dmd -odtmp -ofbuild/hfit -g -gc -gs -L-lgsl -L-lgslcblas hfit.d brent.d powell.d modelScores.d popGenFunc.d mcmc.d

clean:
	rm -r build