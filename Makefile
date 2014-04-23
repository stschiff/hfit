build/hfit: hfit.d brent.d powell.d modelScores.d popGenFunc.d
	dmd -odtmp -ofbuild/hfit -g -debug -L-L/opt/local/lib -L-lgsl -L-lgslcblas hfit.d brent.d powell.d modelScores.d popGenFunc.d