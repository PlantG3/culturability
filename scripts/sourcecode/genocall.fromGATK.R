###############################################################################
### function to call genotypes
###############################################################################
genocall <- function(x, homo=2, homo.perc=0.9, heter=1, heter.perc=0.1,
                     ref.code=1, alt.code=2, heter.code=3, overall.perc=0.9, 
                     missing.code="-", nocall.code=0, help=F) {
if (help) {
  cat("x include two columns of 1.REF 2.ALT"); cat("\n")
  cat("homo: minimum reads for homozygous call"); cat("\n")
  cat("homo.perc: minimum read percentage for homozygous call"); cat("\n")
  cat("heter: minimum reads for heterozygous call"); cat("\n")
  cat("heter.perc: minimum read percentages for heterozygous call"); cat("\n")
  cat("ref.code: code of ref homozygous genotype"); cat("\n")
  cat("alt.code: code of alt homozygous genotype"); cat("\n")
  cat("heter.code: code of heterozygous genotype"); cat("\n")
  cat("missing.code: code of missing reads support (0 read)"); cat("\n")
  cat("no.code: code of no calls but with reads support"); cat("\n")
}
  names(x) <- c("REF", "ALT")
  refc <- as.numeric(as.character(x["REF"]))
  altc <- as.numeric(as.character(x["ALT"]))
  total <- refc + altc
  hetero.criteria <- (refc >= heter & refc/total >= heter.perc
                      & altc >= heter & altc/total >= heter.perc
                      & (refc+altc)/total >= overall.perc)
  ### homo criteria:
  homo.ref.criteria <- (refc >= homo & refc/total >= homo.perc & (refc+altc)/total >= overall.perc)
  homo.alt.criteria <- (altc >= homo & altc/total >= homo.perc & (refc+altc)/total >= overall.perc)
  geno <- nocall.code
  if (refc == 0 & altc == 0) { geno <- missing.code }
  if (hetero.criteria) { geno <- heter.code }
  if (homo.ref.criteria) { geno <- ref.code }
  if (homo.alt.criteria) { geno <- alt.code }
  return(geno)
}
