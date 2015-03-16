#!/bin/bash
# Here is an example line to run the FindGroups code
# 110 = number of species
# ../Data/TatooshAll.txt = path to adjacency matrix
# 587184 = random seed
# 5000 = number of MCMC steps (should probably be higher for a network this size, maybe 100000 and with multiple runs)
# 20 = number of MCMC chains
./FindGroups 110 ../Data/TatooshAll.txt 587184 5000 20
