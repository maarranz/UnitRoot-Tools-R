Rcpp::sourceCpp("dfgls_kernel.cpp")   # run from data-raw/ after installing the package
.libPaths(c("/tmp/rlib", .libPaths())); ns <- asNamespace("unitroottools")
set.seed(777); R <- 50000
for (tr in c("c", "ct")) for (T in c(45, 90, 175, 264, 350)) {
  ks <- c(0L, 3L, 7L)
  st <- dfgls_null_stats(as.integer(T), ks, tr == "ct", R)
  for (j in seq_along(ks)) {
    cv <- ns$dfgls_critical_values(tr, T, ks[j])
    rej <- colMeans(outer(st[, j], cv, "<"))
    cat(sprintf("%-2s T=%4d k=%d  size at 1/5/10%%: %.4f %.4f %.4f   (MC se at 5%%: %.4f)\n",
                tr, T, ks[j], rej[1], rej[2], rej[3], sqrt(.05*.95/R)))
  }
}

# Large lag orders
set.seed(888); R <- 40000
for (tr in c("c", "ct")) for (T in c(175, 350, 600)) {
  ks <- c(13L, 16L, 20L)
  st <- dfgls_null_stats(as.integer(T), ks, tr == "ct", R)
  for (j in seq_along(ks)) {
    cv <- ns$dfgls_critical_values(tr, T, ks[j]); rej <- colMeans(outer(st[, j], cv, "<"))
    cat(sprintf("%-2s T=%4d k=%2d  size at 1/5/10%%: %.4f %.4f %.4f\n", tr, T, ks[j], rej[1], rej[2], rej[3]))
  }
}
