z <- readRDS("trial_a50_all_1_30.rds"); x <- z$res
x <- x[x$year == 2025, ]; x$ok <- x$tier %in% c("ok", "warning"); x$e <- log(x$est/x$tru_eb)
mods <- c("annual", "annual_prot", "annual_lag50", "quarter_lag50", "month_lag50")
grp <- c("pil.27.8abd" = "early", "pra.27.3a4a" = "early", "her.27.25-2932" = "early", "hke.27.8c9a" = "intermediate", "ple.27.420" = "intermediate",
         "hom.27.2a4a5b6a7a-ce-k8" = "intermediate", "bss.27.8ab" = "intermediate", "cod.27.5a" = "late", "anf.27.3a46" = "late")
x$group <- grp[x$stock]
f <- function(d){ o <- d[d$ok, ]; if(!nrow(o)) return("-")
  sprintf("%+.2f/%.2f/%.2f%s", median(o$e), sqrt(mean(o$e^2)), mean(o$tru_eb >= o$lci & o$tru_eb <= o$uci, na.rm = TRUE),
          ifelse(mean(d$ok) < 1, sprintf(" (%.0f%%)", 100*mean(d$ok)), "")) }
short <- c("pil.27.8abd" = "sardine", "pra.27.3a4a" = "shrimp", "her.27.25-2932" = "herring", "hke.27.8c9a" = "hake", "ple.27.420" = "plaice",
           "hom.27.2a4a5b6a7a-ce-k8" = "horse mackerel", "bss.27.8ab" = "sea bass", "cod.27.5a" = "cod", "anf.27.3a46" = "NS anglerfish")
cat("Terminal (2025) error vs true EB/EBmsy: median / RMSE / 95% coverage (convergence if < 100%)\n")
tab <- t(sapply(names(short), function(s) sapply(mods, function(m) f(x[x$stock == s & x$model == m, ]))))
rownames(tab) <- short; print(noquote(tab))
cat("\nBy group\n"); tg <- t(sapply(c("early", "intermediate", "late"), function(g) sapply(mods, function(m) f(x[x$group == g & x$model == m, ]))))
print(noquote(rbind(tg, all = sapply(mods, function(m) f(x[x$model == m, ])))))
