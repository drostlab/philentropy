
context("Test implementation of distance measures...")

test_that("distance(method = 'euclidean') computes the correct distance value.", {
        
        expect_identical(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")), sqrt(sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))^2)))
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "euclidean")))
        #expect_error(phylentropy::distance(1:10, 20:29, method = "euclidean"), "Error : Your probability values are not between: [0,1].")
})


test_that("distance(method = 'manhattan') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "manhattan")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "manhattan")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "manhattan")))
})


test_that("distance(method = 'minkowski') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "minkowski", p = 4)), (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))^4))^0.25)
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "minkowski", p = 4)), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "minkowski", p = 4)))
        
})


test_that("distance(method = 'chebyshev') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "chebyshev")), max(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "chebyshev")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "maximum")))
        
})


test_that("distance(method = 'sorensen') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "sorensen")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum((1:10/sum(1:10)) + (20:29/sum(20:29))))
        
})


test_that("distance(method = 'gower') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "gower")), (1/length(1:10)) * sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        
})






