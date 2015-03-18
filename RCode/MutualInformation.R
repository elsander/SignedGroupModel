## Author: Elizabeth Sander (esander@uchicago.edu)
## Version: 1.0.0
##
## MutInf calculates the mutual information (in nats)
## between two partitions. If ent = TRUE, the function
## will also return the entropies (in nats) of each partition.
##
## allMIs is a utility wrapper for MutInf. It takes in
## a data.frame of partitions, with each partition stored
## in a column. It then returns all pairwise MIs in a matrix.
## Note that MI is a symmetric metric (i.e., MI(A,B) = MI(B,A)),
## so the matrix contains values only in the upper triangular.
## MI(A,A) = H(A), where H(.) is the entropy in nats.

MutInf <- function(A, B, ent = FALSE){
    ##make sure that the partitions are in vector form
    A <- as.vector(as.matrix(A))
    B <- as.vector(as.matrix(B))
    if(length(A) != length(B)) stop('partitions must have same number of nodes')
    
    S <- length(A)
    ##number of modules in A and B
    Amod <- length(unique(A))
    Bmod <- length(unique(B))
    
    uniqA <- unique(A)
    uniqB <- unique(B)
    niA <- unlist(lapply(uniqA, function(x)return(length(A[A==x]))))
    njB <- unlist(lapply(uniqB, function(x)return(length(B[B==x]))))

    nijAB <- data.frame(Amod = rep(unique(A), each = Bmod),
                        Bmod = rep(unique(B), times = Amod))

    MI <- 0
    for(i in 1:(Amod*Bmod)){
        ABoverlap <- length(intersect(which(A == nijAB$Amod[i]),
                                      which(B == nijAB$Bmod[i])))
        if(ABoverlap > 0){
            MI <- MI + ABoverlap/S*log(ABoverlap*S/
                                       (niA[which(uniqA == nijAB$Amod[i])]
                                        *njB[which(uniqB == nijAB$Bmod[i])]))
        }
    }
    
    ##choose whether to return MI only or entropies as well
    if(ent){
        Apart <- unlist(lapply(niA, S=S, function(x,S){
            return(x/S*log(x/S))}))
        entA <- -1*sum(Apart)
        Bpart <- unlist(lapply(njB, S=S, function(x,S){
            return(x/S*log(x/S))}))
        entB <- -1*sum(Bpart)
        return(c(entA, entB, MI))
    } else {
        return(MI)
    }
}

allMIs <- function(parts){
    pnames <- names(parts)
    nparts <- ncol(parts)
    MIs <- matrix(0, nparts, nparts)
    for(i in 1:(nparts-1)){
        for(j in (i+1):nparts){
            MIs[i,j] <- MutInf(parts[,i], parts[,j])
        }
    }

    if(!is.null(pnames)) row.names(MIs) <- colnames(MIs) <- pnames
    return(MIs)
}
