source('MutualInformation.R')

Jackknife <- function(part1, part2){
    ##calculate ratio between MI and max MI possible (=lower entropy of the two)

    MIall <- MutInf(part1, part2, ent = TRUE)
    minEnt <- min(MIall[1:2])
    return(MIall[3]/minEnt)
}
