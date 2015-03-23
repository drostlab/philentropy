
context("Test implementation of distance measures...")

test_that("distance(method = 'euclidean') computes the correct distance value.", {
        
        expect_identical(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")), sqrt(sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))^2)))
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "euclidean")))
})


test_that("distance(method = 'manhattan') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "manhattan")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "manhattan")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "manhattan")))
})


test_that("distance(method = 'Minkowski') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "minkowski", p = 4)), (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))^4))^0.25)
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "minkowski", p = 4)), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "minkowski", p = 4)))
        
})