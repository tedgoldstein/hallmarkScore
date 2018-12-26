library(dplyr)
library(pracma)

hallmark_columns = c(
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

convertGeneNamesToGene_Id = function(d) {

  #convert gene Symbol to geneID
  genename <- as.character(rownames(d))
  GeneID <- as.character(with(Mapgene, geneID[match(genename,gene)]))
  d <- cbind(GeneID,d)
  #remove non-value Gene ID
  d <- subset(d, GeneID != "NA")

  if ( any(d[2:nrow(d),2:ncol(d)] > 1000) ) {
    #it appears to be raw counts
    d = d %>% group_by(GeneID) %>% summarise_all(funs(sum))
    e = "Aggregating duplicate rows by summing counts"
  } else if ( all(d[2:nrow(d),2:ncol(d)] >= 0 & d[2:nrow(d),2:ncol(d)] <= 20) ) {
    #it appears to be normalized log in some way
    d = d %>% group_by(GeneID) %>% summarise_all(funs(geometric_mean))
    e = "Aggregating duplicate rows by geometric mean averaging"

  } else {
    d = d %>% group_by(GeneID) %>% summarise_all(funs(mean))
    e = "Aggregating duplicate rows by averaging"
  }
  din = as.data.frame(d)
  rownames(din) <- as.character(din[,1])
  din[,1] <- NULL

  gid = as.character(d$GeneID)
  HGene <- as.character(with(Hgenes, gene[match(gid,geneID)]))
  MGene <- as.character(with(Mgenes, gene[match(gid,geneID)]))
  dout <- cbind(MGene,HGene,d)
  dout = as.data.frame(dout)
  rownames(dout) <- as.character(dout[,3])
  dout[,1:3] <- NULL
  return(dout)
}

#' cancers Function
#'
#' This function returns a list of the currently supported cancers for the computeSignatureScore
#' @examples
#' cancers()
cancers = function() return(names(Signatures$index))

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
    X = convertGeneNamesToGene_Id(X)

    signaturesForTissue <-  Signatures$signatures[Signatures$index[[cancer]]]

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


Y = computeSignatureScore(X, "Brain")


