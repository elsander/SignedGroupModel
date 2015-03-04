extern int PrintMatrixShort(FILE * F, gsl_matrix_short * A);
extern int PrintVectorFloat(FILE * F, gsl_vector * V);
extern int PrintVectorShort(FILE * F, gsl_vector_short * V);
extern int Partition_Initialize(gsl_vector_short * V,
								int N,
								gsl_rng * r);
extern int RGF(int N, gsl_vector_short * Chain, gsl_vector_short * RGFswap);
extern double Partition_Marginal(gsl_vector_short * V,
								 gsl_vector_short * VCopy,
								 gsl_vector_short * RGFswap,
								 gsl_matrix_short * Adj,
								 int N,
								 gsl_vector * lgammaLookup,
								 gsl_vector * logLookup);
extern int TrySwap(int N,
				   gsl_matrix_short * Adj,
				   gsl_vector_short * HotChain, gsl_vector_short * ColdChain,
				   gsl_vector_short * ChainCopy,
				   gsl_vector_short * RGFswap,
				   double HotTemp, double ColdTemp,
				   gsl_rng * r,
				   gsl_vector * lgammaLookup,
				   gsl_vector * logLookup);
