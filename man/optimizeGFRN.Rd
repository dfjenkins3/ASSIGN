% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimizeGFRN.R
\name{optimizeGFRN}
\alias{optimizeGFRN}
\title{Optimize GFRN gene lists lengths}
\usage{
optimizeGFRN(indata, correlation, correlationList, run = c("akt", "bad",
  "egfr", "her2", "igf1r", "krasgv", "krasqh", "raf"),
  run_ASSIGN_only = FALSE, correlation_only = FALSE,
  keep_optimized_only = FALSE, pathway_lengths = c(seq(5, 20, 5), seq(25,
  275, 25), seq(300, 500, 50)), iter = 1e+05, burn_in = 50000)
}
\arguments{
\item{indata}{The list of data frames from ComBat.step2}

\item{correlation}{A matrix of data to correlate ASSIGN predictions to.
The number of rows should be the same and in the same order as indata}

\item{correlationList}{A list that shows which columns of correlation should
be used for each pathway. See below for more details}

\item{run}{specifys the pathways to predict. The default list will
cause all eight pathways to be run in serial. Specify a pathway ("akt",
"bad", "egfr", etc.) or list of pathways to run those pathways only.}

\item{run_ASSIGN_only}{a logical value indicating if you want to run the
ASSIGN predictions only. Use this to parallelize ASSIGN runs across a compute
cluster or across compute threads}

\item{correlation_only}{a logical value indicating if you want to run the
correlation step only. The function will find the ASSIGN runs in the cwd and
optimize them based on the correlation data matrix.}

\item{keep_optimized_only}{a logical value indicating if you want to keep
all of the ASSIGN run results, or only the runs that provided the optimum
ASSIGN correlations. This will delete all directories in the current working
directory that match the pattern "_gene_list". The default is FALSE}

\item{pathway_lengths}{The gene list lengths that should be run. The default
is the 20 pathway lengths that were used in the paper, but this list can
be customized to which pathway lengths you are willing to accept}

\item{iter}{The number of iterations in the MCMC.}

\item{burn_in}{The number of burn-in iterations. These iterations are
discarded when computing the posterior means of the model parameters.}
}
\value{
ASSIGN runs are output to the current workingdirectory. This function
returns the correlation data and the optimized gene lists that you can use
with runassignGFRN to try these lists on other data.
}
\description{
This function runs ASSIGN pathway prediction on gene list lengths from 5 to
500 to find the optimum gene list length for the GFRN pathways by correlating
the ASSIGN predictions to a matrix of correlation data that you provide. This
function takes a long time to run because you are running ASSIGN many times
on many pathways, so I recommend parallelizing by pathway or running the
ASSIGN predictions first (long and parallelizable) and then running the
correlation step (quick) separately.
}
\examples{
\dontrun{
testData <- read.table("https://dl.dropboxusercontent.com/u/62447/ASSIGN/icbp_Rsubread_tpmlog.txt",
                       sep='\\t', row.names=1, header=1)
corData <- read.table("https://dl.dropboxusercontent.com/u/62447/ASSIGN/proteomics.txt",
                      sep='\\t', row.names=1, header=1)
corData$negAkt <- -1 * corData$Akt
corData$negPDK1 <- -1 * corData$PDK1
corData$negPDK1p241 <- -1 * corData$PDK1p241

corList <- list(akt=c("Akt","PDK1","PDK1p241"),
                bad=c("negAkt","negPDK1","negPDK1p241"),
                egfr=c("EGFR","EGFRp1068"),
                her2=c("HER2","HER2p1248"),
                igf1r=c("IGFR1","PDK1","PDK1p241"),
                krasgv=c("EGFR","EGFRp1068"),
                krasqh=c("EGFR","EGFRp1068"),
                raf=c("MEK1","PKCalphap657","PKCalpha"))

combat.data <- ComBat.step2(testData, pcaPlots = TRUE)

optimization_results <- optimizeGFRN(combat.data, corData, corList)
}

}

