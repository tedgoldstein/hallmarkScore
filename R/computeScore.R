library(dplyr)
library(pracma)

options(shiny.maxRequestSize=50*1024^2) 

    "Evading_growth_suppressors",
    "Sustaining_proliferative_signaling",
    "Reprogramming_energy_metabolism",
    "Resisting_cell_death",
    "Genome_instability",
    "Sustained_angiogenesis",
    "Tissue_invasion_and_metastasis",
    "Tumor_promoting_inflammation",
    "Replicative_immortality",
    "Evading_immune_destruction")

rank.normalize <- function(x, FUN=qnorm, ties.method = "average", na.action) {
    if (missing(na.action)) {
        na.action <- get(getOption("na.action"))
    }
    if(! is.function(na.action)) {
        stop("'na.action' must be a function")
    }
    x <- na.action(x)
    ret = FUN(rank(x, ties.method = ties.method)/(length(x)+1))
    ret
}

#' computeSignatureScore Function
#'
#' This function computes the oncology models fidelity score based on the hallmarks of cancer
#' @param X gene dataframe expression dataset either in z-score or log2(n+1) normalized 
#' @param cancer string that shows the names of the cancer
#' @keywords score
#' @export
#' @examples
#' computeSignatureScore(df, "")
computeSignatureScore = function(X, cancer) {
    signaturesForTissue <- Filter(function(ss) ss$cancer == cancer, Signatures$signatures)

    possible = as.character(row.names(X))
    X = apply(X, 2, function(x) scale(rank.normalize(x), scale=TRUE, center=TRUE))

    row.names(X) <- possible
    X <- data.frame(X)
    scores = data.frame()
    
    n = length(signaturesForTissue)
    signature <- NULL
 
    for (i in 1:n) {
        # Increment the progress bar, and update the detail text.
        # incProgress(1/n, detail = paste("Doing part", i, "of", n))

        signature    <- signaturesForTissue[[i]];
        hallmark <- signature$hallmark;
        
        should  <- names(signature$w)
        genes    <- as.character(intersect(should, possible))

        # printf("should=%d possible=%d actual=%d\n", length(should),length(possible),length(genes));

        score = data.frame();
        posScale <- signature$posScale;
        negScale <- signature$negScale;
        w = signature$w[genes]
    
        XX <- t(X[genes,])
        #cat(XX);
      
    
    
        raw = -XX %*% w + signature$b;
        #heat= XX * w + signature$b;
    for (j in 1:length(raw)) {
            value = raw[j];
            if (value < 0) {
                score[1,j] = round(500  - (negScale * raw[j]));
            } else {
                score[1,j] = round( (posScale * raw[j]) + 500);
            }
            score[1,j] = max(score[1,j], 1) # never less than one
        }
        scores = rbind(scores, score);
    }


    scores = t(scores)
    rownames(scores) = colnames(X);
    colnames(scores) = unlist(lapply(signaturesForTissue,function(sig) sig$hallmark))
    scores = scores[,1:10]

    Hallmark = apply(scores, 1, function(x)  round(exp(mean(log(x)))))

    return (scores)
}

