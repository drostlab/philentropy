context("Test implementation of dist_one_one(), dist_one_many(), and dist_many_many() ...")

set.seed(2020-08-20)
P <- 1:10 / sum(1:10)
Q <- 20:29 / sum(20:29)
M1 <- t(replicate(10, sample(1:10, size = 10) / 55))
M2 <- t(replicate(20, sample(1:10, size = 10) / 55))

doo1 <- dist_one_one(P, Q, method = "euclidean", testNA = FALSE)
dom1 <- dist_one_many(P, M1, method = "euclidean", testNA = FALSE)
dmm1 <- dist_many_many(M1, M2, method = "euclidean", testNA = FALSE)

test_that("dist_one_one output structure is correct", {
  expect_type(doo1, "double")
  expect_length(doo1, 1)
})

test_that("dist_one_many output structure is correct", {
  expect_type(dom1, "double")
  expect_length(dom1, nrow(M1))
})

test_that("dist_many_many output structure is correct", {
  expect_type(dmm1, "double")
  expect_equal(dim(dmm1), c(nrow(M1), nrow(M2)))
})

doo2 = euclidean(P, Q, FALSE)

test_that("dist_one_one output is correct", {
  expect_equal(doo1, doo2)
})

dom2 = vector(length = nrow(M1))
for (i in seq_len(nrow(M1))){
  dom2[i] = euclidean(P, M1[i, ], FALSE)
}

test_that("dist_one_many output is correct", {
  expect_equal(dom1, dom2)
})

dmm2 = matrix(nrow = nrow(M1), ncol = nrow(M2))
for (i in seq_len(nrow(M1))){
  for (j in seq_len(nrow(M2))){
    dmm2[i, j] = euclidean(M1[i, ], M2[j, ], FALSE)
  }
}

test_that("dist_many_many output is correct", {
  expect_equal(dmm1, dmm2)
})