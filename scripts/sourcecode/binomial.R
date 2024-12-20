#### function to do binomial test ####
binomial.od <- function(data, chrcol=1, poscol=2, cat1pos, cat2pos, geno, lambda = 1) {
  library(dispmod)
  chr <- data[chrcol]
  pos <- data[poscol]
  cat1counts <- data[cat1pos]
  cat1counts <- as.numeric(cat1counts)
  cat2counts <- data[cat2pos]
  cat2counts <- as.numeric(cat2counts)
  geno <- as.character(geno)
  # check the levels of geno
  geno.level <- unique(geno)
  if (length(geno.level) == 2) {
    geno <- as.factor(geno)
    full <- glm(cbind(cat1counts, cat2counts) ~ geno,
            family = binomial(link=logit)) ## Binomial test
    reduced <- glm(cbind(cat1counts, cat2counts) ~ 1,
			family = binomial(link=logit)) ## Binomial test
	logm2b <- full[[1]][2] # log odd ratio
    log2m2b <- log2(exp(logm2b))
    ### method: suggest by Haiyan Wang
    full.od <-try(glm.binomial.disp(full, verbose = F), silent=TRUE)
	if (is(full.od, "try-error")) {
      cis.odp.control.pvalue <- NA
    } else {
	    res <- anova(full.od, test = 'Chisq')
		cis.odp.control.pvalue <- 1 - pchisq(res[2,2] / lambda, res[2,1])
	}
  } else {
	cis.odp.control.pvalue <- NA
    log2m2b <- NA
  }
  # output variables:
  return(c(data, log2m2b, cis.odp.control.pvalue))
}

