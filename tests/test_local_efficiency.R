#!/usr/bin/env Rscript
# Local efficiency must follow Latora & Marchiori (2001):
#   E_loc(v) = mean over pairs (i,j) in N(v) of 1 / d(i,j)
# where d is a shortest path whose edge lengths are 1/|weight| -- a stronger
# coupling is a SHORTER path, not a longer one.
#
# The old implementation used |weight| directly as distance and then took
# 1/mean(1/d), the harmonic mean of distances, inverting the metric twice.
# Those two forms coincide for a uniform complete neighbour subgraph, which is
# why the bug was not obvious; this fixture uses a 2-hop path so they differ.

suppressMessages({library(here); library(igraph)})
source(here("src", "r", "compat.R"))

# v has neighbours A, B, C.  Among them: A-B (w=4), B-C (w=4), no A-C edge.
edges <- c("v","A", "v","B", "v","C", "A","B", "B","C")
g <- graph_from_edgelist(matrix(edges, ncol = 2, byrow = TRUE), directed = TRUE)
E(g)$weight <- c(1, 1, 1, 4, 4)

local_efficiency_node <- function(g, v) {
  nb <- neighbors(g, v, mode = "all")
  if (length(nb) < 2) return(0)
  sg <- induced_subgraph(g, nb)
  if (ecount(sg) == 0) return(0)
  E(sg)$weight <- 1 / abs(E(sg)$weight)
  d <- distances(sg, weights = E(sg)$weight)
  diag(d) <- NA
  mean(1 / d, na.rm = TRUE)
}

# d(A,B)=1/4, d(B,C)=1/4, d(A,C)=1/2  ->  efficiencies 4, 4, 2  ->  mean = 10/3
expected <- 10 / 3
got <- local_efficiency_node(g, which(V(g)$name == "v"))

# What the old formula produced, for contrast.
sg <- induced_subgraph(g, neighbors(g, which(V(g)$name == "v"), mode = "all"))
E(sg)$weight <- abs(E(sg)$weight)
d_old <- distances(sg, weights = E(sg)$weight)
d_old[d_old == 0] <- NA
old <- 1 / mean(1 / d_old, na.rm = TRUE)

cat(sprintf("  expected (Latora & Marchiori): %.4f\n", expected))
cat(sprintf("  new implementation:            %.4f\n", got))
cat(sprintf("  old implementation:            %.4f\n", old))

ok <- abs(got - expected) < 1e-9
cat(if (ok) "PASS\n" else "FAIL\n")
if (!ok) quit(status = 1)
if (abs(old - expected) < 1e-9) {
  cat("WARNING: fixture does not distinguish old from new\n"); quit(status = 1)
}
cat(sprintf("PASS: fixture separates them by %.4f\n", abs(old - expected)))
