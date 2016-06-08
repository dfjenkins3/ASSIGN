#' Optimize GFRN gene lists lengths
#'
#' This function runs ASSIGN pathway prediction on gene list lengths from 5 to
#' 500 to find the optimum gene list length for the GFRN pathways by correlating
#' the ASSIGN predictions to a matrix of correlation data that you provide. This
#' function takes a long time to run because you are running ASSIGN many times
#' on many pathways, so I recommend parallelizing by pathway or running the
#' ASSIGN predictions first (long and parallelizable) and then running the
#' correlation step (quick) separately.
#'
#' @param indata The list of data frames from ComBat.step2
#' @param correlation A matrix of data to correlate ASSIGN predictions to.
#' The number of rows should be the same and in the same order as indata
#' @param correlationList A list that shows which columns of correlation should
#' be used for each pathway. See below for more details
#' @param run specifys the pathways to predict. The default list will
#' cause all eight pathways to be run in serial. Specify a pathway ("akt",
#' "bad", "egfr", etc.) or list of pathways to run those pathways only.
#' @param run_ASSIGN_only a logical value indicating if you want to run the
#' ASSIGN predictions only. Use this to parallelize ASSIGN runs across a compute
#' cluster or across compute threads 
#' @param correlation_only a logical value indicating if you want to run the
#' correlation step only. The function will find the ASSIGN runs in the cwd and
#' optimize them based on the correlation data matrix.
#' @param keep_optimized_only a logical value indicating if you want to keep
#' all of the ASSIGN run results, or only the runs that provided the optimum
#' ASSIGN correlations
#' @param pathway_lengths The gene list lengths that should be run. The default
#' is the 20 pathway lengths that were used in the paper, but this list can
#' be customized to which pathway lengths you are willing to accept
#' @param iter The number of iterations in the MCMC.
#' @param burn_in The number of burn-in iterations. These iterations are
#' discarded when computing the posterior means of the model parameters.
#' 
#' @return ASSIGN runs and correlation data are output to the current working
#' directory. This function returns the optimized gene lists that you can use
#' with runassignGFRN to try these lists on other data.
#'
#' @export optimizeGFRN
#'
optimizeGFRN <- function(indata, correlation, correlationList,
                         run=c("akt","bad","egfr","her2", "igf1r", "krasgv",
                               "krasqh", "raf"), run_ASSIGN_only=FALSE,
                         correlation_only=FALSE, keep_optimized_only=FALSE,
                         pathway_lengths=c(seq(5,20,5), seq(25,275,25),
                                           seq(300,500,50)), iter=100000,
                         burn_in=50000) {
  if(!(correlation_only)){
    for (curr_path in run){
      for (curr_len in pathway_lengths){
        geneList <- list()
        geneList[[curr_path]] <- c(gfrn_geneList[[paste(curr_path,"up",sep="_")]][1:floor(curr_len/2)],
                                   gfrn_geneList[[paste(curr_path,"down",sep="_")]][1:ceiling(curr_len/2)])
        runassignGFRN(indata, run=curr_path, optimized_geneList=geneList,
                      iter=iter, burn_in=burn_in)   
      }
    }
  }
}
