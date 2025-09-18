# get_new_seqs

library(tidyverse)
library(rentrez)


## (optional but speeds you up 3-fold)
key <- Sys.getenv("ENTREZ_KEY", unset = NA)

# get seqs
taxid <- 138950  # Enterovirus C :contentReference[oaicite:1]{index=1}
query  <- paste0("txid", taxid, "[Organism:exp]")

search <- entrez_search(db = "nuccore",
                        term = query,
                        retmax = 0,      # first call just to get record count
                        use_history = TRUE)
search$count        
wh <- search$web_history

gb_txt = entrez_fetch("nuccore", web_history = wh,
                      rettype = "gb", retmode = "text",
                      retstart = k, retmax = chunk_size)