# SignedGroupModel

Elizabeth Sander and Stefano Allesina

## Overview

This package contains code to search for and analyze group model 
structure of signed adjacency networks.

## Search Algorithm

All code to search for network groupings can be found in the `CCode` folder. 
To compile the code, run `make` from terminal while in the folder. See
`CCode/example.sh` for an explanation of how to run the algorithm.

We use a Metropolis-Coupled Markov Chain Monte Carlo (MCMCMC) search 
algorithm with a Gibbs sampler to search for the optimal grouping for 
a given network. Bayes factors are used for model selection.

## Collecting Taxonomic Data

Code to create taxonomic records for a species list from the ITIS database
is in the `PythonCode` folder. Detailed instructions can be found in 
`PythonCode/README-GetTaxonomy.txt`.

## Partition Analysis

Several functions to compare and analyze network groupings can be found in 
the `RCode` folder. `MutualInformation.R` contains functions to calculate 
the  similarity between partitions. `JackknifeMIs.R` calculates the MI
between two partitions relative to the maximum possible MI. 
`RandomizationPvals.R` calculates the p-value associated with a MI
value (here, the p-value is the probability of getting an equal or
higher MI when the group identities in the partitions are randomized).


License: MIT