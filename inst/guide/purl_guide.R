# regenerate inst/guide/JABBA3_user.R from the guide chapters (purl documentation = 2)
# run from the package root
old <- readLines("inst/guide/JABBA3_user.R"); hdr <- old[seq_len(grep("^## ----setup", old)[1] - 1)]
main <- readLines("inst/guide/JABBA3_user_guide.Rmd")
ch <- regmatches(main, regexpr('chapters/[^"]+\\.Rmd', main)); ch <- file.path("inst/guide", ch)
out <- c("## ----setup--------------------------------------------------------------------", "library(JABBA3)", "library(ggplot2)", "")
for(f in ch){ o <- tempfile(fileext = ".R"); knitr::purl(f, o, documentation = 2, quiet = TRUE); out <- c(out, readLines(o), "") }
writeLines(c(hdr, out), "inst/guide/JABBA3_user.R")
