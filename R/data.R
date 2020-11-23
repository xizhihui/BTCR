#' vdjdb slim dataset
#'
#' A data frame from vdjdb.slim.txt
#'
#' @format A data frame with 16 columns
#' \describe{
#'     \item{gene}{TCR chain: alpha or beta}
#'     \item{cdr3}{TCR complementarity determining region 3 (CDR3) amino acid sequence}
#'     \item{species}{TCR parent species}
#'     \item{antigen.epitope}{Amino acid sequence of the epitope}
#'     \item{antigen.gene}{Representative parent gene of the epitope}
#'     \item{antigen.species}{Representative parent species of the epitope}
#'     \item{complex.id}{TCR alpha and beta chain records having the same complex identifier belong to the same T-cell clone}
#'     \item{v.segm}{TCR Variable segment allele}
#'     \item{j.segm}{TCR Joining segment allele}
#'     \item{v.end}{End of v.segm}
#'     \item{j.start}{Start of j.segm}
#'     \item{mhc.a}{First MHC chain allele}
#'     \item{mhc.b}{Second MHC chain allele (defaults to Beta2Microglobulin for MHC class I)}
#'     \item{mhc.class}{MHC class MHC class (I or II)}
#'     \item{reference.id}{Pubmed reference / URL / or submitter details in case unpublished}
#'     \item{vdjdb.score}{VDJdb confidence score, the higher is the score the more confidence we have in the antigen specificity annotation of a given TCR clonotype/clone. Zero score indicates that there are insufficient method details to draw any conclusion}
#' }
#' @source \url{https://github.com/antigenomics/vdjdb-db}
#'
"vdjdb"
