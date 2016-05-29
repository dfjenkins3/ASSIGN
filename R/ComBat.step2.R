#' Perform the second step of ComBat
#'
#' The first ComBat step (on the signatures only) has already been
#' performed. This step performs batch correction on the test data,
#' using reference batch ComBat, to prepare the test data for ASSIGN
#' analysis.
#'
#' @param testData The input test data to batch correct
#' @param pcaPlots a logical value indicating whether or not the function
#' should create PCA plots. The default is FALSE.
#'
#' @return A list of data.frames is returned, including control (GFP) and signature
#' data, as well as the batch corrected test data. This data can go directly into
#' the runassign.single and runassign.multi functions, or subsetted to go directly
#' into ASSIGN.
#'
#' @export ComBat.step2
ComBat.step2 <- function(testData, pcaPlots=FALSE) {
  dat <- merge_drop(combat_train,testData)
  sub <- c(6,6,12,6,6,5,6,6,9,9,9,9,ncol(testData))
  bat <- c(rep(1,ncol(combat_train)),rep(2,ncol(testData)))
  if(pcaPlots){
    grDevices::pdf("pca_refcombat_twostep.pdf")
    pcaplot(dat,sub)
  }
  combat_expr1 <- sva::ComBat(dat=dat, batch=bat, mod=NULL, ref.batch=1)
  if(pcaPlots){
    pcaplot(combat_expr1,sub)
    invisible(grDevices::dev.off())
  }
  c_gfp <- subset(combat_expr1, select=GFP.1:GFP.12)
  c_akt <- subset(combat_expr1, select=AKT.1:AKT.6)
  c_bad <- subset(combat_expr1, select=BAD.1:BAD.6)
  c_her2 <- subset(combat_expr1, select=HER2.1:HER2.6)
  c_igf1r <- subset(combat_expr1, select=IGF1R.1:IGF1R.6)
  c_raf <- subset(combat_expr1, select=RAF.1:RAF.6)
  c_egfr_gfp <- combat_expr1[,1:6]
  c_egfr <- combat_expr1[,7:12]
  c_kras_gfp <- subset(combat_expr1,select=GFP.31:GFP.39)
  c_kraswt <- subset(combat_expr1,select=KRASWT.1:KRASWT.9)
  c_krasqh <- subset(combat_expr1,select=KRASQH.1:KRASQH.9)
  c_krasgv <- subset(combat_expr1,select=KRASGV.1:KRASGV.9)
  c_test <- combat_expr1[,(ncol(combat_train)+1):ncol(combat_expr1)]
  results <- list(gfp=c_gfp, akt=c_akt, bad=c_bad, her2=c_her2, igf1r=c_igf1r,
                  raf=c_raf, egfr_gfp=c_egfr_gfp, egfr=c_egfr,
                  kras_gfp=c_kras_gfp, kraswt=c_kraswt, krasqh=c_krasqh,
                  krasgv=c_krasgv, test=c_test)
  return(results)
}