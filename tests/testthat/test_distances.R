
context("Test implementation of distance measures...")

test_that("distance(method = 'euclidean') computes the correct distance value.", {
        
        expect_identical(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")), sqrt(sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))^2)))
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "euclidean")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "euclidean")))
        #expect_error(philentropy::distance(1:10, 20:29, method = "euclidean"), "Error : Your probability values are not between: [0,1].")
})


test_that("distance(method = 'manhattan') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "manhattan")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "manhattan")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "manhattan")))
})


test_that("distance(method = 'minkowski') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "minkowski", p = 4)), (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))^4))^0.25)
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "minkowski", p = 4)), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "minkowski", p = 4)))
        
})


test_that("distance(method = 'chebyshev') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "chebyshev")), max(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "chebyshev")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "maximum")))
        
})


test_that("distance(method = 'sorensen') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "sorensen")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum((1:10/sum(1:10)) + (20:29/sum(20:29))))
        
})


test_that("distance(method = 'gower') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "gower")), (1/length(1:10)) * sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))))
        
})


test_that("distance(method = 'soergel') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "soergel")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,max)))
        
})


test_that("distance(method = 'kulczynski_d') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "kulczynski_d")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,min)))
        
})


test_that("distance(method = 'canberra') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "canberra")), sum( abs((1:10/sum(1:10)) - (20:29/sum(20:29))) / ((1:10/sum(1:10)) + (20:29/sum(20:29)))))
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "canberra")), as.vector(stats::dist(base::rbind(1:10/sum(1:10),20:29/sum(20:29)), method = "canberra")))
        
})


test_that("distance(method = 'lorentzian') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "lorentzian")), sum( log(1 + abs((1:10/sum(1:10)) - (20:29/sum(20:29))))))
        
})


test_that("distance(method = 'intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "intersection")), sum(apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,min)))
        
})


test_that("distance(method = 'non-intersection') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "non-intersection")), 1 - sum(apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,min)))
        
})


test_that("distance(method = 'wavehedges') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "wavehedges")), sum(abs(1:10/sum(1:10) - 20:29/sum(20:29)) / apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,max)))
        
})


test_that("distance(method = 'czekanowski') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "czekanowski")), sum(abs(1:10/sum(1:10) - 20:29/sum(20:29))) / sum(1:10/sum(1:10) + 20:29/sum(20:29)))
        
})


test_that("distance(method = 'motyka') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "motyka")), sum(apply(base::rbind(1:10/sum(1:10),20:29/sum(20:29)),2,max)) / sum(1:10/sum(1:10) + 20:29/sum(20:29)))
        
})


test_that("distance(method = 'kulczynski_s') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "kulczynski_s")), 1 / (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,min))))
        
})


test_that("distance(method = 'tanimoto') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "tanimoto")), sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,max)))
        
})


test_that("distance(method = 'ruzicka') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "ruzicka")), 1 - (sum(abs((1:10/sum(1:10)) - (20:29/sum(20:29)))) / sum(apply(rbind(1:10/sum(1:10), 20:29/sum(20:29)),2,max))))
        
})


test_that("distance(method = 'inner_product') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "inner_product")), sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) ))
        
})


test_that("distance(method = 'harmonic_mean') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "harmonic_mean")), 2 * sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) / ((1:10/sum(1:10)) + (20:29/sum(20:29))) ))
        
})


test_that("distance(method = 'cosine') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "cosine")), ((sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) )) / (sqrt(sum((1:10/sum(1:10))^2)) * sqrt(sum((20:29/sum(20:29))^2)))))
        
})


test_that("distance(method = 'hassebrook') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "hassebrook")), ((sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) )) / (sum((1:10/sum(1:10))^2) + sum((20:29/sum(20:29))^2) - ((sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) ))))))
        
})



test_that("distance(method = 'jaccard') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "jaccard")), 1 - ((sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) )) / (sum((1:10/sum(1:10))^2) + sum((20:29/sum(20:29))^2) - ((sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) ))))))
        
})


test_that("distance(method = 'dice') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "dice")), 1 - (2 * (sum ( (1:10/sum(1:10)) * (20:29/sum(20:29)) )) / (sum((1:10/sum(1:10))^2) + sum((20:29/sum(20:29))^2) )))
        
})


test_that("distance(method = 'fidelity') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "fidelity")), sum(sqrt(1:10/sum(1:10) * 20:29/sum(20:29))))
        
})


test_that("distance(method = 'bhattacharyya') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "bhattacharyya")), -log(sum(sqrt(1:10/sum(1:10) * 20:29/sum(20:29)))))
        
})


test_that("distance(method = 'hellinger') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "hellinger")), 2 * sqrt(1 - sum(sqrt(1:10/sum(1:10) * 20:29/sum(20:29)))))
        
})



test_that("distance(method = 'matusita') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "matusita")), sqrt(sum((sqrt(1:10/sum(1:10)) - sqrt(20:29/sum(20:29)))^2)))
        
})


test_that("distance(method = 'squared_chord') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "squared_chord")), sum((sqrt(1:10/sum(1:10)) - sqrt(20:29/sum(20:29)))^2))
        
})


test_that("distance(method = 'squared_euclidean') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "squared_euclidean")), sum(((1:10/sum(1:10)) - (20:29/sum(20:29)))^2))
        
})


test_that("distance(method = 'pearson') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "pearson")), sum(((1:10/sum(1:10)) - (20:29/sum(20:29)))^2 / (20:29/sum(20:29))))
        
})


test_that("distance(method = 'neyman') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "neyman")), sum(((1:10/sum(1:10)) - (20:29/sum(20:29)))^2 / (1:10/sum(1:10))))
        
})


test_that("distance(method = 'squared_chi') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "squared_chi")), sum(((1:10/sum(1:10)) - (20:29/sum(20:29)))^2 / ((1:10/sum(1:10)) + (20:29/sum(20:29)))))
        
})


test_that("distance(method = 'prob_symm') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "prob_symm")), 2 * sum(((1:10/sum(1:10)) - (20:29/sum(20:29)))^2 / ((1:10/sum(1:10)) + (20:29/sum(20:29)))))
        
})


test_that("distance(method = 'divergence') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "divergence")), 2 * sum(((1:10/sum(1:10)) - (20:29/sum(20:29)))^2 / ((1:10/sum(1:10)) + (20:29/sum(20:29)))^2))
        
})


test_that("distance(method = 'clark') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "clark")), sqrt(sum((abs((1:10/sum(1:10)) - (20:29/sum(20:29))) / ((1:10/sum(1:10)) + (20:29/sum(20:29))))^2)))
        
})


test_that("distance(method = 'additive_symm') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "additive_symm")), sum((((1:10/sum(1:10)) - (20:29/sum(20:29)))^2 * ((1:10/sum(1:10)) + (20:29/sum(20:29)))) / ((1:10/sum(1:10)) * (20:29/sum(20:29)))))
        
})


test_that("distance(method = 'kullback-leibler') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "kullback-leibler")), sum((1:10/sum(1:10)) * log((1:10/sum(1:10)) / (20:29/sum(20:29)))))
        
})


test_that("distance(method = 'jeffreys') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "jeffreys")), sum(((1:10/sum(1:10)) - (20:29/sum(20:29))) * log((1:10/sum(1:10)) / (20:29/sum(20:29)))))
        
})



test_that("distance(method = 'k_divergence') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "k_divergence")), sum((1:10/sum(1:10)) * log(2 * (1:10/sum(1:10)) / ((1:10/sum(1:10)) + (20:29/sum(20:29))))))
        
})


test_that("distance(method = 'topsoe') computes the correct distance value.", {
        
        expect_equal(as.vector(philentropy::distance(1:10/sum(1:10), 20:29/sum(20:29), method = "topsoe")), sum(((1:10/sum(1:10)) * log(2 * (1:10/sum(1:10)) / ((1:10/sum(1:10)) + (20:29/sum(20:29))))) + ((20:29/sum(20:29)) * log(2 * (20:29/sum(20:29)) / ((1:10/sum(1:10)) + (20:29/sum(20:29)))))))
        
})





