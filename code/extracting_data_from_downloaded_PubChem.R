library(ChemmineR)
library(jsonlite)
library(R.utils)
library(stringi)
library(tidyverse)

setwd("C:/Users/annel/OneDrive - Kruvelab/Kruvelab/R codes for variouse things/Computational resources/PubChem")

get_PubChem_finger <- function(PubChem_as_SDF) {
  tryCatch(
    blockmatrix <- datablock2ma(datablocklist=datablock(PubChem_as_SDF)), # Converts data block to matrix 
    error = function(e) print("index not in the range")
  )
  tryCatch(
    numchar <- splitNumChar(blockmatrix=blockmatrix), # Splits to numeric and character matrix 
    error = function(e) print("index not in the range")
  )
  tryCatch(cbind(SDFID = sdfid(PubChem_as_SDF),
                 numchar[[1]][,1],
                 numchar[[2]][,12],
                 numchar[[2]][,1]),
           error = function(e) print("index not in the range"))
}

data_path <- "C:/Users/annel/OneDrive - Kruvelab/Kruvelab/R codes for variouse things/Computational resources/PubChem"
files <- dir(data_path, pattern = "*.sdf") # get file names

for (filename in files) {
  name_output <- paste("fingerprint_", filename, ".csv", col="", sep="")
  sdfStream(input = filename, 
            output = name_output, 
            append = FALSE, 
            fct = tryCatch(get_PubChem_finger, 
                            error = function(e) print("index not in the range"),
            Nlines = 1000))
}


