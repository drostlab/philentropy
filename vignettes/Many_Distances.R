## -----------------------------------------------------------------------------
library(philentropy)

## -----------------------------------------------------------------------------
P <- 1:10 / sum(1:10)
Q <- 20:29 / sum(20:29)

## -----------------------------------------------------------------------------
# install.packages("microbenchmark")
microbenchmark::microbenchmark(
  dist(rbind(P, Q), method = "euclidean"),
  distance(rbind(P, Q), method = "euclidean", test.na = FALSE, mute.message = TRUE),
  euclidean(P, Q, FALSE),
  dist_one_one(P, Q, method = "euclidean", testNA = FALSE)
)

## -----------------------------------------------------------------------------
set.seed(2020-08-20)
P <- 1:10 / sum(1:10)
M <- t(replicate(100, sample(1:10, size = 10) / 55))

## -----------------------------------------------------------------------------
# install.packages("microbenchmark")
microbenchmark::microbenchmark(
  as.matrix(dist(rbind(P, M), method = "euclidean"))[1, ][-1],
  distance(rbind(P, M), method = "euclidean", test.na = FALSE, mute.message = TRUE)[1, ][-1],
  dist_one_many(P, M, method = "euclidean", testNA = FALSE)
)

## -----------------------------------------------------------------------------
set.seed(2020-08-20)
M1 <- t(replicate(10, sample(1:10, size = 10) / 55))
M2 <- t(replicate(10, sample(1:10, size = 10) / 55))

## -----------------------------------------------------------------------------
many_dists = function(m1, m2){
  r = matrix(nrow = nrow(m1), ncol = nrow(m2))
  for (i in seq_len(nrow(m1))){
    for (j in seq_len(nrow(m2))){
      x = rbind(m1[i, ], m2[j, ])
      r[i, j] = distance(x, method = "euclidean", mute.message = TRUE)
    }
  }
  r
}

## -----------------------------------------------------------------------------
# install.packages("microbenchmark")
microbenchmark::microbenchmark(
  many_dists(M1, M2),
  dist_many_many(M1, M2, method = "euclidean", testNA = FALSE)
)

