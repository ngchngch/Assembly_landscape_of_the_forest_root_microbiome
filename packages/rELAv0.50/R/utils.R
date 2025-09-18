#'applymap: supporting function for Formatting
applymap <- function(df, fun){
  return(t(apply(df, 1, function(row){sapply(row, fun)})))
}

#'applyvec: supporting function for Formatting
applyvec <- function(df, fun){
  if(length(rownames(df)) > 1){
    output <- t(data.frame(lapply(t(df), fun))) %>% 'rownames<-'(rownames(df))
  }else(
    output <- data.frame(lapply(df, fun)) %>% 'colnames<-'(colnames(df))
  )
  return(output)
}

#'T: supporting function for Formatting
T <- function(df){return(as.data.frame(t(df)))}

#  Arguments: baseabtable, basemetadata: DataFrame
# normalize: 0 or 1
# parameters: list
# Returns: ocmatrix, abmatrix, fctmatrix: DataFrame
# samplelabel, taxonlabel, factorlabel: list

#'Formatting: generate data matrices from rawdata
#' @useDynLib rELA, .registration=TRUE
#' @export
Formatting <- function(baseabtable, basemetadata=NULL, normalize, parameters, grouping=0, grouping_th=0, reporting=TRUE){
  ath <- parameters[1]     # Relative frequency threshold; if less than this value, set to 0
  minocth <- parameters[2] # Lower limit of mean occurence; species with mean occurence lower than this value are excluded.
  maxocth <- parameters[3] # Upper limit of mean occurence; species with mean occurence higher than this value are excluded.
  
  if(length(basemetadata) != 0){
    # with basemetadata
    # Get label information
    specieslabel = colnames(baseabtable)
    absamplelabel = rownames(baseabtable)
    mdsamplelabel = rownames(basemetadata)
    sharedsamplelabel = sort(intersect(absamplelabel, mdsamplelabel))

    # Generating occurence matrix
    patab = match(sharedsamplelabel, absamplelabel)
    patmd = match(sharedsamplelabel, mdsamplelabel)
    abnum = baseabtable[patab,]
    if (normalize == 1){ # Nomalization within samples (optional)
      abnum <- t(apply(abnum, 1, function(x){x / (sum(x) + 10^-16)}))
    }
    aboc = applymap(abnum, function(x){if(x >= ath){1}else{0}}) # 0/1 transformation with ath
    aboc_mean = apply(aboc, 2, mean)
    occrit = names(aboc_mean[(aboc_mean > minocth) & (aboc_mean < maxocth)]) # Calculation of mean occurence//determination of the target to be removed

    # Object generation
    ocmatrix = aboc[, occrit] # Occurence matrix
    abmatrix = abnum[, occrit] # Abundance matrix
    fctmatrix = as.matrix(as.data.frame(basemetadata[patmd, ])) # Environmental factors
    rownames(fctmatrix) <- rownames(ocmatrix) 
    colnames(fctmatrix) <- colnames(basemetadata)
    samplelabel = sharedsamplelabel # Sample labels
    specieslabel = occrit # Species labels
    factorlabel = colnames(basemetadata) # Factor labels
  }else{
    # without basemetadata
    # Get label information
    specieslabel = colnames(baseabtable)
    absamplelabel = rownames(baseabtable)
    sharedsamplelabel = sort(absamplelabel)

    # Generating occurence matrix
    patab = match(sharedsamplelabel, absamplelabel)
    abnum = baseabtable[patab,]
    if (normalize == 1){ # Nomalization within samples (optional)
      abnum <- t(apply(abnum, 1, function(x){x / (sum(x) + 10^-16)}))
    }
    aboc = applymap(abnum, function(x){if(x >= ath){1}else{0}}) # 0/1 transformation with ath
    aboc_mean = apply(aboc, 2, mean)
    occrit = names(aboc_mean[(aboc_mean > minocth) & (aboc_mean < maxocth)]) # Calculation of mean occurence//determination of the target to be removed

    # Object generation
    ocmatrix = aboc[, occrit] # Occurence matrix
    abmatrix = abnum[, occrit] # Abundance matrix
    fctmatrix = NULL # Environmental factors
    samplelabel = sharedsamplelabel # Sample labels
    specieslabel = occrit # Species labels
    factorlabel = NULL # Factor labels
  }

  # Grouping
  if(grouping == 1){
    ocvstr <- apply(ocmatrix, 2, function(y){return(paste(y, collapse = ""))})
    oo <- order(ocvstr)
    ov <- ocvstr[oo]
    j <- 0
    cc <- 0
    binlist <- list()
    aaa <- foreach (i=ov) %do% {
        cc <- cc + 1
        if(stringdist(i, j, method='hamming')/nchar(i) <= grouping_th){a <- rev(binlist)[[1]]
                aa <- append(a, cc)
                binlist <- append(binlist[-length(binlist)], list(aa))
                }else{binlist <- append(binlist, list(cc))}
        j=i
        }
    ltx <- length(specieslabel)
    specieslabel <- foreach (i=binlist, .combine='c') %do% {paste(specieslabel[oo][i], collapse=", ")}
    ocmatrix <- foreach(i=binlist, .combine='cbind') %do% {ocmatrix[,oo][,i[1]]}
    abmatrix <- foreach(i=binlist, .combine='cbind') %do% {if(length(i)>1){apply(abmatrix[,oo][,i],1,sum)}else{abmatrix[,oo][,i]}}
    colnames(ocmatrix) <- specieslabel
    colnames(abmatrix) <- specieslabel
  }

  if(reporting){
    cat('Processed', as.character(length(samplelabel)), 'samples.\n')
    cat('Relative abundance threshold =', as.character(ath), '\n')
    cat('Occurrence threshold (lower) =', as.character(minocth), '\n')
    cat('Occurrence threshold (upper) =', as.character(maxocth), '\n')
    if(grouping == 1){cat(as.character(ltx-length(specieslabel)), ' groups were found.\n')}
    cat('Selected ', as.character(length(specieslabel)), ' out of ',
                as.character(ncol(baseabtable)), 'species.\n')
  }
  return(list(ocmatrix, abmatrix, fctmatrix, samplelabel, specieslabel, factorlabel))}


#'dec2bin: supporting function for CInegerDigits
dec2bin <- function(num, digit=0){
  if(num <= 0 && digit <= 0){
    return(NULL)
  }else{
    return(append(Recall(num%/%2,digit-1), num%%2))
  }
}

#'CInegerDigits: map an integer to a binary vector
#' @useDynLib rELA, .registration=TRUE
#' @export
CIntegerDigits <- function(ssid, n, index=NULL){
  if(is.null(index)){index <- 1:n}
  return(data.frame(dec2bin(ssid, n), row.names=index))
  }


#'str2i: map a binary vector to an integer (for large integers)
#' @export
str2i <- function(x) {
    y <- as.numeric(strsplit(x, "")[[1]])
    sum(y * 2^rev((seq_along(y)-1)))
}

#'bin2id: map a binary vector to a 36base characters
#' @export
bin2id <- function(bin){
  digs <- "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_/"
  bi <- unlist(bin)
  modi <- length(bi)%%6
  if(modi != 0){b0 <- rep(0, 6-modi)}else{b0<-c()}
  bii <- c(b0, bi)
  lb <- length(bii)
  aa <- foreach(i=seq(1, length(bii), 6), .combine="c") %do% {a <- rev(bii)[i:(i+5)]
                                s <- str2i(paste(as.character(rev(a)), collapse=""))
                                substr(digs, s+1, s+1)}
  return(paste(rev(aa), collapse=""))}

#'id2bin: map a 36base characters to a binary vector
#' @export
id2bin <-  function(id, s){
    digs <- "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz_/"
    aa <- foreach (i=seq(1:nchar(id)), .combine="rbind") %do% {
        a <- substr(id,i,i)
        ia <- (unlist(gregexpr(a, digs))[1]) - 1
        CIntegerDigits(ia, 6)}
    uaa <- unlist(aa)
    aaa <- as.vector(uaa[(length(uaa)-s+1):length(uaa)])
    return(aaa)}

