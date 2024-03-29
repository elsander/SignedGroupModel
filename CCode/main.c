#include "Common.h"
#include "SearchAlgs.h"
#include "partition.h"

//Run multiple-chain Markov Chain Monte Carlo to partition nodes
//of a network into groups
//
//Author: Elizabeth Sander (esander@uchicago.edu)
//Version: 1.0.0
//
//Requirements: gcc compiler, GSL and GSL CBLAS libraries for C
//
//Takes in the number of species, an adjacency matrix, a random
//seed, the number of steps for each chain to take, and the number of
//chains.
//A chain swap is attempted every 20 steps. This is hardcoded in the
//global variables in SearchAlgs.c
//
//
//Example function call from terminal:
// ./FindGroupGibbs N ./path/to/file seed steps nchains
//
//Prints a vector of length N to file. This vector contains the best
//partition found by the algorithm. Each integer i in the vector
//represents the group identity of species i.
//The output file will be named "<./path/to/file>-Marginal<marginal loglikelihood>.txt"
//Therefore, you should copy your adjacency matrix to the folder where
//you want your output file.

int main (int argc, char *argv[]){
	int N = atoi(argv[1]); // number of species
	char * FileName = argv[2]; // file storing the signed adjacency mat
	int seed = atoi(argv[3]); // random seed
	int Steps = atoi(argv[4]); // number of steps
	int nChains = atoi(argv[5]); // for MC3, this is the number of
							   // chains

	// Read the adjacency matrix
	gsl_matrix_short * Adj = gsl_matrix_short_calloc(N,N);
	FILE * F;
	F=fopen(FileName, "rb");
	gsl_matrix_short_fscanf(F, Adj);
	fclose(F);

	// Set the random seed
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
	gsl_rng_set (r, seed);
	
	// Best Solution
	double MarginalLikelihood = 0.0;
	// Stuff for output
	char OutFileName[1000];
	gsl_vector_short * BestSolution = gsl_vector_short_calloc(N);

	int maxLog = 0; //used for lookup table later
	int i,j;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			if(gsl_matrix_short_get(Adj, i, j) != 0){
				maxLog++;
			}
		}
	}

	//set up log gamma function lookup table
	//max lgamma is 3 + N*N
	gsl_vector * lgammaLookup = gsl_vector_calloc(3+(N*N));
	for(i=0; i<(3+(N*N)); i++){
		gsl_vector_set(lgammaLookup, i, gsl_sf_lngamma(i+1));
	}
	gsl_vector * logLookup = gsl_vector_calloc(2+maxLog);
	for(i=0; i<(2+maxLog); i++){
		gsl_vector_set(logLookup, i, log(i+1));
	}

	// Run the search
	fprintf(stderr, "starting MC3 with %d steps\n", Steps);
	MarginalLikelihood = MC3(N, Adj, Steps, nChains,
							 BestSolution, r, lgammaLookup, logLookup);
	fprintf(stderr, "MC3 complete. Best solution marginal %.4f\n", MarginalLikelihood);

	// Save the results
	sprintf(OutFileName,"%s-Marginal%f", FileName, MarginalLikelihood);
	F=fopen(OutFileName,"w");
	// print best solution
	PrintVectorShort(F, BestSolution);
	// close file
	fclose(F);
	// Free memory
	gsl_matrix_short_free(Adj);
	gsl_vector_free(lgammaLookup);
	gsl_vector_free(logLookup);
	gsl_vector_short_free(BestSolution);
	gsl_rng_free(r);
	return 0;	
}
