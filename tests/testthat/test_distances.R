
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


test_that("distance(method = 'soergel') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "soergel")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,max)))
        
})


test_that("distance(method = 'kulczynski_d') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "kulczynski_d")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,min)))
        
})


test_that("distance(method = 'canberra') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "canberra")), sum( abs((1:10/sum(1:10)) - (20:29/sum(20:29))) / ((1:10/sum(1:10)) + (20:29/sum(20:29)))))
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "canberra")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "canberra")))
        
})


test_that("distance(method = 'lorentzian') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "lorentzian")), sum( log(1 + abs((1:10/sum(1:10)) - (20:29/sum(20:29))))))
        
})


test_that("distance(method = 'intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "intersection")), sum(apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,min)))
        
})


test_that("distance(method = 'non-intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "non-intersection")), 1 - sum(apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,min)))
        
})


test_that("distance(method = 'wavehedges') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "wavehedges")), sum(abs(1:10/sum(1:10) - 20:29/sum(20:29)) / apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,max)))
        
})


test_that("distance(method = 'czekanowski') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "czekanowski")), sum(abs(1:10/sum(1:10) - 20:29/sum(20:29))) / sum(1:10/sum(1:10) + 20:29/sum(20:29)))
        
})


test_that("distance(method = 'motyka') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "motyka")), sum(apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,max)) / sum(1:10/sum(1:10) + 20:29/sum(20:29)))
        
})


test_that("distance(method = 'kulczynski_s') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "kulczynski_s")), 1 / (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,min))))
        
})


test_that("distance(method = 'tanimoto') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "tanimoto")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,max)))
        
})


test_that("distance(method = 'ruzicka') computes the correct distance value.", {
        
        expect_equal(as.vector(phylentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "ruzicka")), 1 - (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,max))))
        
})



