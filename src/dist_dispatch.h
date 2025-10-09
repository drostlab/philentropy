#ifndef philentropy_DISPATCH_H
#define philentropy_DISPATCH_H philentropy_DISPATCH_H

#include <string>
#include <cmath>
#include "distances_internal.h"

// template for dispatching internal distance methods
template <typename InputIt1, typename InputIt2>
double dispatch_dist_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2,
                              const std::string& method,
                              const std::string& unit,
                              double epsilon, double p) {
    if (method == "minkowski") {
        return minkowski_internal(first1, last1, first2, p);
    } else if (method == "euclidean"){
        return euclidean_internal(first1, last1, first2);
    } else if (method == "manhattan") {
        return minkowski_internal(first1, last1, first2, 1.0);
    } else if (method == "chebyshev") {
        return chebyshev_internal(first1, last1, first2);
    } else if (method == "sorensen") {
        return sorensen_internal(first1, last1, first2);
    } else if (method == "gower") {
        return gower_internal(first1, last1, first2);
    } else if (method == "soergel") {
        return soergel_internal(first1, last1, first2);
    } else if (method == "kulczynski_d") {
        return kulczynski_d_internal(first1, last1, first2, epsilon);
    } else if (method == "canberra") {
        return canberra_internal(first1, last1, first2);
    } else if (method == "lorentzian") {
        return lorentzian_internal(first1, last1, first2, unit);
    } else if (method == "intersection") {
        return intersection_dist_internal(first1, last1, first2);
    } else if (method == "non-intersection") {
        return 1.0 - intersection_dist_internal(first1, last1, first2);
    } else if (method == "wavehedges") {
        return wave_hedges_internal(first1, last1, first2);
    } else if (method == "czekanowski") {
        return czekanowski_internal(first1, last1, first2);
    } else if (method == "motyka") {
        return motyka_internal(first1, last1, first2);
    } else if (method == "kulczynski_s") {
        return 1.0 / kulczynski_d_internal(first1, last1, first2, epsilon);
    } else if (method == "tanimoto") {
        return tanimoto_internal(first1, last1, first2);
    } else if (method == "ruzicka") {
        return ruzicka_internal(first1, last1, first2);
    } else if (method == "inner_product") {
        return inner_product_internal(first1, last1, first2);
    } else if (method == "harmonic_mean") {
        return harmonic_mean_internal(first1, last1, first2);
    } else if (method == "cosine") {
        return cosine_internal(first1, last1, first2);
    } else if (method == "hassebrook") {
        return hassebrook_internal(first1, last1, first2);
    } else if (method == "jaccard") {
        return jaccard_internal(first1, last1, first2);
    } else if (method == "dice") {
        return dice_internal(first1, last1, first2);
    } else if (method == "fidelity") {
        return fidelity_internal(first1, last1, first2);
    } else if (method == "bhattacharyya") {
        return bhattacharyya_internal(first1, last1, first2, unit, epsilon);
    } else if (method == "hellinger") {
        return hellinger_internal(first1, last1, first2);
    } else if (method == "matusita") {
        return matusita_internal(first1, last1, first2);
    } else if (method == "squared_chord") {
        return squared_chord_internal(first1, last1, first2);
    } else if (method == "squared_euclidean") {
        return squared_euclidean_internal(first1, last1, first2);
    } else if (method == "pearson") {
        return pearson_internal(first1, last1, first2, epsilon);
    } else if (method == "neyman") {
        return neyman_internal(first1, last1, first2, epsilon);
    } else if (method == "squared_chi") {
        return squared_chi_internal(first1, last1, first2);
    } else if (method == "prob_symm") {
        return prob_symm_internal(first1, last1, first2);
    } else if (method == "divergence") {
        return divergence_internal(first1, last1, first2);
    } else if (method == "clark") {
        return clark_internal(first1, last1, first2);
    } else if (method == "additive_symm") {
        return additive_symm_internal(first1, last1, first2);
    } else if (method == "kullback-leibler") {
        return kullback_leibler_internal(first1, last1, first2, unit, epsilon);
    } else if (method == "jeffreys") {
        return jeffreys_internal(first1, last1, first2, unit, epsilon);
    } else if (method == "k_divergence") {
        return k_divergence_internal(first1, last1, first2, unit);
    } else if (method == "topsoe") {
        return topsoe_internal(first1, last1, first2, unit);
    } else if (method == "jensen-shannon"){
        return jensen_shannon_internal(first1, last1, first2, unit);
    } else if (method == "jensen_difference") {
        return jensen_difference_internal(first1, last1, first2, unit);
    } else if (method == "taneja") {
        return taneja_internal(first1, last1, first2, unit, epsilon);
    } else if (method == "kumar-johnson") {
        return kumar_johnson_internal(first1, last1, first2, epsilon);
    } else if (method == "avg") {
        return avg_internal(first1, last1, first2);
    }
    return NAN;
}

#endif // philentropy_DISPATCH_H