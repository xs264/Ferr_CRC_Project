
## TCGA is small, we use curl -o Colon_cancer_TCGA_AllSamples.csv https://granddb.s3.amazonaws.com/cancer/colon_cancer/networks/lioness/Colon_cancer_TCGA_AllSamples.csv
## to download all files

## gse is big, we download it one by one
library(httr)
gse = read.csv("data/GSE39582.csv")
geo_id = gse$geo_accession
for (i in 1:length(geo_id)) {
  print(i)
  id = geo_id[i]
  
  url <- paste0("https://granddb.s3.amazonaws.com/cancer/colon_cancer/networks/lioness/Colon_cancer_sample_", id, ".csv")
  destfile <- paste0("data/GEO_colon/lioness_", id, ".csv")
  response <- HEAD(url)
  
  if (status_code(response) == 200) {
    download.file(url, destfile = destfile, mode = "wb")
    #next
  } else {
    message(paste("File does not exist:", id))
  }
}
