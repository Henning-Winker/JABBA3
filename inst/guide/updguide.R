
rmarkdown::render("inst/guide/JABBA3_user_guide.Rmd")
pages <- "../jabba3-pages"

file.copy("inst/guide/JABBA3_user_guide.html", file.path(pages, "index.html"), overwrite = TRUE)

file.copy("inst/guide/JABBA3_user_guide.html", file.path(pages, "docs/index.html"), overwrite = TRUE)

# commit and push the gh-pages worktree
system(sprintf('git -C "%s" add index.html docs/index.html', pages))
system(sprintf('git -C "%s" push origin gh-pages', pages))

system(sprintf('git -C "%s" commit -m "Update user guide"', pages))


system('git -C "../jabba3-pages" commit -m "Update user guide"')
system('git -C "../jabba3-pages" push origin gh-pages')
