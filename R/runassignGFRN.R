#' Run optimized single pathway ASSIGN
#'
#' This function runs eight ASSIGN runs based on the pathway optimizations
#' from the paper. You can run all eight pathways in serial, or call this
#' function and specify the run parameter to run a specific pathway.
#' Some ASSIGN parameters can be customized using this function. The default
#' values were used in the analysis for the paper.
#'
#' @param indata The list of data frames from ComBat.step2
#' @param run specifys the pathways to predict. The default NULL value will
#' cause all eight pathways to be run in serial. Specify a pathway ("akt",
#' "bad", "egfr", etc.) to run that pathway only.
#' @param use_seed a logical value indicating if you want to run the analysis
#' using a set seed. This will make the result consistant between runs. The
#' default is TRUE.
#' @param sigma_sZero Each element of the signature matrix (S) is modeled by a
#' spike-and-slab mixuture distribution. Sigma_sZero is the variance of the
#' spike normal distribution. The default is 0.05.
#' @param sigma_sNonZero Each element of the signature matrix (S) is modeled by
#' a spike-and-slab mixuture distribution. Sigma_sNonZero is the variance of the
#' slab normal distribution. The default is 0.5.
#' @param S_zeroPrior Logicals. If TRUE, the prior distritribution of signature
#' follows a normal distribution with mean zero. The default is FALSE.
#' @param iter The number of iterations in the MCMC. The default is 100000.
#' @param burn_in The number of burn-in iterations. These iterations are discarded
#' when computing the posterior means of the model parameters. The default is 50000.
#'
#' @return Data is output to the current working directory in a results directory.
#'
#' @export runassignGFRN
#'
runassignGFRN <- function(indata, run=NULL, use_seed=TRUE, sigma_sZero=0.05,
                             sigma_sNonZero=0.5, S_zeroPrior=FALSE, iter=100000,
                             burn_in=50000) {
  if(is.null(run) || run=="akt"){
    if(use_seed){
      set.seed(1234)
    }
    #AKT
    ASSIGN::assign.wrapper(trainingData=cbind(indata$gfp,indata$akt),
                           testData=indata$test,
                           anchorGenes=list(akt=c("AKT1")),
                           trainingLabel=list(control=list(akt=1:12),akt=13:18),
                           geneList=list(akt=geneLists$akt),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("akt_",length(geneLists$akt),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="bad"){
    if(use_seed){
      set.seed(1234)
    }
    #BAD
    ASSIGN::assign.wrapper(trainingData=cbind(indata$gfp,indata$bad),
                           testData=indata$test,
                           anchorGenes=list(bad=c("BAD")),
                           trainingLabel=list(control=list(bad=1:12),bad=13:18),
                           geneList=list(bad=geneLists$bad),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("bad_",length(geneLists$bad),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="egfr"){
    if(use_seed){
      set.seed(1234)
    }
    #EGFR
    ASSIGN::assign.wrapper(trainingData=cbind(indata$egfr_gfp,indata$egfr),
                           testData=indata$test,
                           anchorGenes=list(egfr=c("EGFR")),
                           trainingLabel=list(control=list(egfr=1:6),egfr=7:12),
                           geneList=list(egfr=geneLists$egfr),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("egfr_",length(geneLists$egfr),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="her2"){
    if(use_seed){
      set.seed(1234)
    }
    #HER2
    ASSIGN::assign.wrapper(trainingData=cbind(indata$gfp,indata$her2),
                           testData=indata$test,
                           anchorGenes=list(her2=c("ERBB2")),
                           trainingLabel=list(control=list(her2=1:12),her2=13:17),
                           geneList=list(her2=geneLists$her2),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("her2_",length(geneLists$her2),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="igf1r"){
    if(use_seed){
      set.seed(1234)
    }
    #IGF1R
    ASSIGN::assign.wrapper(trainingData=cbind(indata$gfp,indata$igf1r),
                           testData=indata$test,
                           anchorGenes=list(igf1r=c("IGF1R")),
                           trainingLabel=list(control=list(igf1r=1:12),igf1r=13:18),
                           geneList=list(igf1r=geneLists$igf1r),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("igf1r_",length(geneLists$igf1r),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="krasgv"){
    if(use_seed){
      set.seed(1234)
    }
    #KRASGV
    ASSIGN::assign.wrapper(trainingData=cbind(indata$kras_gfp,indata$krasgv),
                           testData=indata$test,
                           anchorGenes=list(krasgv=c("KRAS")),
                           trainingLabel=list(control=list(krasgv=1:9),krasgv=10:18),
                           geneList=list(krasgv=geneLists$krasgv),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("krasgv_",length(geneLists$krasgv),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="krasqh"){
    if(use_seed){
      set.seed(1234)
    }
    #KRASQH
    ASSIGN::assign.wrapper(trainingData=cbind(indata$kras_gfp,indata$krasqh),
                           testData=indata$test,
                           anchorGenes=list(krasqh=c("KRAS")),
                           trainingLabel=list(control=list(krasqh=1:9),krasqh=10:18),
                           geneList=list(krasqh=geneLists$krasqh),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("krasqh_",length(geneLists$krasqh),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
  if(is.null(run) || run=="raf"){
    if(use_seed){
      set.seed(1234)
    }
    #RAF
    ASSIGN::assign.wrapper(trainingData=cbind(indata$gfp,indata$raf),
                           testData=indata$test,
                           anchorGenes=list(raf=c("RAF1")),
                           trainingLabel=list(control=list(raf=1:12),raf=13:18),
                           geneList=list(raf=geneLists$raf),
                           n_sigGene=NULL,
                           adaptive_B=TRUE,
                           adaptive_S=TRUE,
                           mixture_beta=FALSE,
                           S_zeroPrior=S_zeroPrior,
                           outputDir=paste("raf_",length(geneLists$raf),"_gene_list", sep=""),
                           sigma_sZero=sigma_sZero, sigma_sNonZero=sigma_sNonZero,
                           iter=iter, burn_in=burn_in)
  }
}
