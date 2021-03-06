# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

asssessCoactivity <- function(A) {
    .Call(`_NetDECODE_asssessCoactivity`, A)
}

constructKstarNN <- function(logPvals, L_C = 5.0, pval_threshold = 0.01) {
    .Call(`_NetDECODE_constructKstarNN`, logPvals, L_C, pval_threshold)
}

symmetrizeNetwork <- function(G) {
    .Call(`_NetDECODE_symmetrizeNetwork`, G)
}

predictActivityScores <- function(A, rows, columns, G) {
    .Call(`_NetDECODE_predictActivityScores`, A, rows, columns, G)
}

