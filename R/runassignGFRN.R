#' Run optimized single pathway ASSIGN
#'
#' This function runs eight ASSIGN runs based on the pathway optimizations
#' from the paper. You can run all eight pathways in serial, or call this
#' function and specify the run parameter to run a specific pathway.
#' Some ASSIGN parameters can be customized using this function. The default
#' values were used in the analysis for the paper.
#'
#' @param indata The list of data frames from ComBat.step2
#' @param run specifys the pathways to predict. The default list will
#' cause all eight pathways to be run in serial. Specify a pathway ("akt",
#' "bad", "egfr", etc.) or list of pathways to run those pathways only.
#' @param optimized_geneList a list of custom optimized gene lists for the gfrn
#' pathways either created manually or output by optimizeGFRN
#' @param use_seed Set the seed before running ASSIGN. This will make the result
#' consistant between runs. The default is 1234. Set use_seed as FALSE to not
#' set a seed.
#' @param sigma_sZero Each element of the signature matrix (S) is modeled by a
#' spike-and-slab mixuture distribution. Sigma_sZero is the variance of the
#' spike normal distribution. The default is 0.05.
#' @param sigma_sNonZero Each element of the signature matrix (S) is modeled by
#' a spike-and-slab mixuture distribution. Sigma_sNonZero is the variance of the
#' slab normal distribution. The default is 0.5.
#' @param S_zeroPrior Logicals. If TRUE, the prior distritribution of signature
#' follows a normal distribution with mean zero. The default is FALSE.
#' @param iter The number of iterations in the MCMC. The default is 100000.
#' @param burn_in The number of burn-in iterations. These iterations are
#' discarded when computing the posterior means of the model parameters. The
#' default is 50000.
#'
#' @return Data is output to the current working directory in a results
#' directory.
#' 
#' @examples
#' \dontrun{
#' testData <- read.table("https://dl.dropboxusercontent.com/u/62447/ASSIGN/icbp_Rsubread_tpmlog.txt",
#'                        sep='\t', row.names=1, header=1)
#' combat.data <- ComBat.step2(testData, pcaPlots = TRUE)
#' runassignGFRN(combat.data)
#' }
#'
#' @export runassignGFRN
#'
runassignGFRN <- function(indata, run=c("akt","bad","egfr","her2","igf1r",
                                        "krasgv","krasqh","raf"),
                          optimized_geneList=NULL, use_seed=1234,
                          sigma_sZero=0.05, sigma_sNonZero=0.5,
                          S_zeroPrior=FALSE, iter=100000, burn_in=50000) {

  #list of anchor genes
  anchorGeneList <- list(akt="AKT1", bad="BAD", egfr="EGFR", her2="ERBB2",
                         igf1r="IGF1R", krasgv="KRAS", krasqh="KRAS",
                         raf="RAF1")
  
  #list of corresponding controls for each pathway
  gfpList <- list(akt="gfp", bad="gfp", egfr="egfr_gfp", her2="gfp",
                  igf1r="gfp", krasgv="kras_gfp", krasqh="kras_gfp", raf="gfp")
  
  if(is.null(optimized_geneList)){
    utils::data('gfrn_geneList', package='ASSIGN', envir=environment()) 
    gfrn_geneList <- get("gfrn_geneList", envir=environment())
    optimized_geneList=list(akt=c(gfrn_geneList$akt_up[1:10],
                                  gfrn_geneList$akt_down[1:10]),
                            bad=c(gfrn_geneList$bad_up[1:125],
                                  gfrn_geneList$bad_down[1:125]),
                            egfr=c(gfrn_geneList$egfr_up[1:25],
                                   gfrn_geneList$egfr_down[1:25]),
                            her2=c(gfrn_geneList$her2_up[1:5],
                                   gfrn_geneList$her2_down[1:5]),
                            igf1r=c(gfrn_geneList$igf1r_up[1:50],
                                    gfrn_geneList$igf1r_down[1:50]),
                            krasgv=c(gfrn_geneList$krasgv_up[1:87],
                                     gfrn_geneList$krasgv_down[1:88]),
                            krasqh=c(gfrn_geneList$krasqh_up[1:150],
                                     gfrn_geneList$krasqh_down[1:150]),
                            raf=c(gfrn_geneList$raf_up[1:175],
                                  gfrn_geneList$raf_down[1:175]))
  }
  
  for (curr_path in run){
    trainingLabel <- list()
    trainingLabel[['control']] <- list()
    trainingLabel[['control']][[curr_path]] <- 1:
      ncol(indata[[gfpList[[curr_path]]]])
    trainingLabel[[curr_path]] <- (ncol(indata[[gfpList[[curr_path]]]])+1):
      (ncol(indata[[gfpList[[curr_path]]]])+ncol(indata[[curr_path]]))
    
    if(!(anchorGeneList[curr_path] %in% rownames(indata[['test']]))){
      warning(anchorGeneList[curr_path], " not in input data. No anchor gene ",
              "will be used.")
      anchorGeneList[curr_path] <- list(NULL)
    }
    
    if(use_seed){
      set.seed(use_seed)
    }
    
    assign.wrapper(trainingData=cbind(indata[[gfpList[[curr_path]]]],
                                      indata[[curr_path]]),
                   testData=indata[['test']],
                   anchorGenes=anchorGeneList[curr_path],
                   trainingLabel=trainingLabel,
                   geneList=optimized_geneList[curr_path],
                   n_sigGene=NULL,
                   adaptive_B=TRUE,
                   adaptive_S=TRUE,
                   mixture_beta=FALSE,
                   S_zeroPrior=S_zeroPrior,
                   outputDir=paste(curr_path,"_",
                                   length(optimized_geneList[[curr_path]]),
                                   "_gene_list", sep=""),
                   sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                   iter=iter, burn_in=burn_in)
  }
}
