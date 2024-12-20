starRCmerger <- function(datapath = ".", selection = NULL, suffix = "ReadsPerGene.out.tab") {
### script to merge all counts data from STAR output
### Sanzhen Liu
### 5/4/2017
### datapath: directory for read count files
### suffix: suffix of file names
  count.files <- dir(path = datapath, pattern = suffix)
  if (!is.null(selection)) {
  	count.files <- grep(selection, count.files, value = T)
  }
  ### merge all counts
  allcounts <- NULL
  for (cf in count.files) {
    counts <- read.delim(paste0(datapath, "/", cf), header = F, stringsAsFactors = F, skip = 4)
    base <- gsub(suffix, "", cf)
    counts <- counts[, 1:2]
    colnames(counts) <- c("Gene", base)
    
    ### merge data
    if (is.null(allcounts)) {
      allcounts <- counts
    } else {
      allcounts <- merge(allcounts, counts, by = "Gene")
    }
  }
  
  ### output
  allcounts
}
