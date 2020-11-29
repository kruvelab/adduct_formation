library(jsonlite)
library(R.utils)
library(stringi)

fingerprint_base64 = "AAADcQBAOAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACgAACAAAAAAAgAAACAAAAgAIAACQCAIAAAAAAAAAAABAAAABAAAAAAAAAAAAAAABAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA=="



fingerprint_hex_to_dataframe <- function(fingerprint_base64) {
  fingerprint_hex <- base64_dec(fingerprint_base64)
  fingerprint_with_pre <- R.utils::intToBin(fingerprint_hex)
  fingerprint <- fingerprint_with_pre[5:115]
  fingerprint <- stri_join(fingerprint, collapse = "")
  fingerprint <- strsplit(fingerprint[1], "")
  fingerprint <- as_tibble(t(as.integer(fingerprint[[1]])))
  return(fingerprint)
}
