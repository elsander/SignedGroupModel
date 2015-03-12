## Author: Elizabeth Sander (esander@uchicago.edu)
## Version: 1.0.0
##
## Calculate similarity between partition and jackknifed partition
## relative to the max possible similarity
## Used in Sander et al. (in review) to determine sensitivity of
## groupings to inclusion/exclusion of individual species
## Note that this code assumes that the two partitions are the same length!

source('MutualInformation.R')

Jackknife <- function(part1, part2){
    ##calculate ratio between MI and max MI possible (=lower entropy of the two)

    MIall <- MutInf(part1, part2, ent = TRUE)
    minEnt <- min(MIall[1:2])
    return(MIall[3]/minEnt)
}
