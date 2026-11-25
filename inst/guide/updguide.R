# Render the user guide, regenerate the R script and publish to gh-pages. Run from the package root.
# Needs tmbstan/rstan for the NUTS chapters.
rmarkdown::render("inst/guide/JABBA3_user_guide.Rmd")
source("inst/guide/purl_guide.R")                       # inst/guide/JABBA3_user.R from the chapters
pages <- "../jabba3-pages"                               # gh-pages worktree
system(sprintf('git -C "%s" pull --ff-only origin gh-pages', pages))
file.copy("inst/guide/JABBA3_user_guide.html", file.path(pages, "index.html"), overwrite = TRUE)
file.copy("inst/guide/JABBA3_user_guide.html", file.path(pages, "docs/index.html"), overwrite = TRUE)
system(sprintf('git -C "%s" add index.html docs/index.html', pages))
system(sprintf('git -C "%s" commit -m "Update user guide"', pages))
system(sprintf('git -C "%s" push origin gh-pages', pages))
