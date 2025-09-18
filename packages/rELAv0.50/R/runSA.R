#'Stochastic Approximation
#'@description Estimating species interaction, environemnt preference.
#'
#'@param ocmat : array of the vectors of presence/absence status for each sample. Row is sample, column is species.
#'@param enmat : array of the vectors of environmental factors for each sample. Normalization required.
#'@param qth : threshold parameter update size as a criterion for stopping the optimization.
#'@param rep : number of parallel processes for the optimization.
#'@param threads : number of cores (threads) used for calculation. 
#'@param reporting : enable/disable the reporting function.
#'
#' @useDynLib rELA, .registration=TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom foreach foreach
#' @export
runSA <- function(ocmat, enmat=NULL, qth=10^-5, rep=128, threads=1, reporting=TRUE){

    cluster = makeCluster(threads)
    registerDoParallel(cluster)
    on.exit(stopCluster(cluster))

    if(is.null(enmat)){

    	## ============================================== ##
        ## -- without explicit variables

        if(reporting){cat('Start parameter fitting:\n')}
        s <- proc.time()[3]
        fittingMat <- foreach(i = 1:rep, .packages="rELA") %dopar% {

	      if(reporting){cat( sprintf('Trial %s...',i) )}
            SAres <- SA(ocData = ocmat, qTh = qth)
            return(SAres)
        }

        if(reporting){cat(sprintf('SA: elapsed time %.2f sec\n\n', proc.time()[3]-s))}
        ## ============================================== ##
		fittingRes <- matrix(0, ncol=ncol(ocmat)+1, nrow=ncol(ocmat),
        		     dimnames=list(colnames(ocmat), c('h', paste('J',colnames(ocmat),sep='.')) ))

    	for(i in 1:rep)fittingRes <- fittingRes + fittingMat[[i]]

	    }else{

        ## ============================================== ##
        ## -- with explicit variables

        if(reporting){cat('Start parameter fitting\n')}
        s <- proc.time()[3]
        fittingMat <- foreach(i = 1:rep, .packages="rELA") %dopar% {

            if(reporting){cat( sprintf('Trial %s...',i) )}
            SAres <- fullSA(ocData = ocmat, envData=as.matrix(enmat), qTh = qth)
            return(SAres)
        }

        if(reporting){cat(sprintf('SA: elapsed time %.2f sec\n\n', proc.time()[3]-s))}
        ## ============================================== ##
	if( ncol(as.matrix(enmat))==1 ){
	    gname <- paste('g', colnames(enmat), sep='.')
	    ncols <- ncol(ocmat)+2
	}else{
	    gname <- paste('g',1:(1+(ncol(enmat)-1)), sep='.')
	    ncols <- ncol(ocmat)+2+(ncol(enmat)-1)
	}

        fittingRes <- matrix(0, ncol= ncols, nrow=ncol(ocmat),
        		     dimnames=list(colnames(ocmat), c('h', gname,
					   paste('J',colnames(ocmat),sep='.')) ))

    	for(i in 1:rep)fittingRes <- fittingRes + fittingMat[[i]]
    }

    return(list(fittingRes/rep, enmat))
}


#'sa2paramss: extract model parameters from runSA output
#' @useDynLib rELA, .registration=TRUE
sa2params <- function(sa, env=NULL){
  sadim <- dim(sa[[1]])
  na <- rownames(sa[[1]])
  if(sadim[2] > sadim[1]+1){
      if(length(env)==0){env <- rep(0, sadim[2]-sadim[1]-1)} ## !!!
        he <- sa[[1]][, 1] %>% 'names<-'(na)
        je <- sa[[1]][, ((sadim[2]-sadim[1])+1):sadim[2]] %>% 'dimnames<-'(list(na,na))
        ge <- as.matrix(sa[[1]][,2:(sadim[2]-sadim[1])])  %>% 'dimnames<-' (list(na,colnames(sa[[2]])))
        hge <- (he + ge %*% env)[,1]
  }else{
        he <- sa[[1]][, 1] %>% 'names<-'(na)
        je <- sa[[1]][, -1] %>% 'dimnames<-'(list(na,na))
        ge <- NULL
        hge <- he
  }
  return(list(he, je, ge, hge))}
