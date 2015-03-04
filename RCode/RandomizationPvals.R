source('MutualInformation.R')

RandomizationPvals <- function(part1, part2, reps = 1000000){
    ##Performs reps randomizations of two partitions to
    ##generate a null distribution of Mutual Information values
    ##Then calculates a p-value based on this null.
    ##part1 and part2 are vectors containing two groupings for
    ##the same set of species.
    
    S <- length(part1) ##should be the same as in part2
    if(length(part1) != length(part2)) stop('Partition lengths must be identical')

    oldSim <- MutInf(part1, part2)

    nMoreSimilar <- 0

    ##Using lapply for speed here
    ##limiting list size to 10000 for memory reasons
    ##the number of reps will be rounded up to the
    ##nearest 10000
    if(reps %% 10000) print('# of reps will be rounded up to the nearest 10,000')
    samppart1 <- rep(list(part1), 10000)
    samppart2 <- rep(list(part2), 10000)
    for(i in 1:ceiling(reps/10000)){
        print(i*10000)
        ##this has the partitions sample themselves every replicate
        ##but since it's just reshuffling the original partition over
        ##and over, the only bias this will introduce will be a result
        ##of the randomizer used
        samppart1 <- lapply(samppart1, sample, size = S,
                            replace = FALSE)
        samppart2 <- lapply(samppart2, sample, size = S,
                            replace = FALSE)

        similarities <- unlist(mapply(MutInf,
                                      A = samppart1,
                                      B = samppart2))

        nMoreSimilar <- nMoreSimilar + sum(similarities >= oldSim)
    }

    SimPval <- nMoreSimilar/(ceiling(reps/10000)*10000)

    return(SimPval)
}
