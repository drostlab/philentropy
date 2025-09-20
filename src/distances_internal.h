
#ifndef philentropy_Distances_Internal_H
#define philentropy_Distances_Internal_H philentropy_Distances_Internal_H

#include <Rcpp.h> 
#include <math.h>
#include <iostream>
#include "utils.h"

template <typename InputIt1, typename InputIt2>
double fidelity_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2);

template <typename InputIt1, typename InputIt2>
double squared_chi_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2);

template <typename InputIt1, typename InputIt2>
double k_divergence_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit);

template <typename InputIt1, typename InputIt2>
double topsoe_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit);

template <typename InputIt1, typename InputIt2>
double jensen_shannon_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit);

template <typename InputIt1, typename InputIt2>
double chebyshev_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2);

template <typename InputIt1, typename InputIt2>
double euclidean_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double diff = 0.0;

    while (first1 != last1) {
        diff = *first1 - *first2;
        dist += diff * diff;
        first1++;
        first2++;
    }
    return std::sqrt(dist);
}

template <typename InputIt1, typename InputIt2>
double minkowski_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, double p) {
    double dist = 0.0;
    double diff = 0.0;

    while (first1 != last1) {
        diff = *first1 - *first2;
        dist += std::pow(std::abs(diff), p);
        first1++;
        first2++;
    }
    return std::pow(dist, 1.0 / p);
}


template <typename InputIt1, typename InputIt2>
double lorentzian_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit) {
    double dist = 0.0;
    double diff = 0.0;

    while (first1 != last1) {
        diff = *first1 - *first2;
        dist += std::log(1.0 + std::abs(diff));
        first1++;
        first2++;
    }

    if (unit == "log2")
        return dist / std::log(2.0);
    if (unit == "log10")
        return dist / std::log(10.0);

    return dist;
}


template <typename InputIt1, typename InputIt2>
double sorensen_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist1 = 0.0;
    double dist2 = 0.0;
    double diff = 0.0;
    double sum = 0.0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        sum = *first1 + *first2;
        dist1 += diff;
        dist2 += sum;
        first1++;
        first2++;
    }

    if (dist2 == 0.0) {
        return NAN;
    } else {
        return dist1 / dist2;
    }
}


template <typename InputIt1, typename InputIt2>
double gower_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double diff = 0.0;
    int n = 0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        dist += diff;
        first1++;
        first2++;
        n++;
    }

    return (1.0 / n) * dist;
}


template <typename InputIt1, typename InputIt2>
double soergel_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist1 = 0.0;
    double dist2 = 0.0;
    double diff = 0.0;
    double max_point = 0.0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        if (*first1 >= *first2) {
            max_point = *first1;
        } else {
            max_point = *first2;
        }
        dist1 += diff;
        dist2 += max_point;
        first1++;
        first2++;
    }

    if (dist2 == 0.0) {
        return 0;
    } else {
        return dist1 / dist2;
    }
}


template <typename InputIt1, typename InputIt2>
double kulczynski_d_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, double epsilon) {
    double dist1 = 0.0;
    double dist2 = 0.0;
    double diff = 0.0;
    double min_point = 0.0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        if (*first1 <= *first2) {
            min_point = *first1;
        } else {
            min_point = *first2;
        }
        dist1 += diff;
        if (min_point == 0.0) {
            dist2 += epsilon;
        } else {
            dist2 += min_point;
        }
        first1++;
        first2++;
    }

    if (dist2 == 0.0) {
        return NAN;
    } else {
        return dist1 / dist2;
    }
}


template <typename InputIt1, typename InputIt2>
double tanimoto_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    return soergel_internal(first1, last1, first2);
}


template <typename InputIt1, typename InputIt2>
double ruzicka_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    return 1.0 - soergel_internal(first1, last1, first2);
}


template <typename InputIt1, typename InputIt2>
double inner_product_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;

    while (first1 != last1) {
        dist += *first1 * *first2;
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double harmonic_mean_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double prod = 0.0;
    double sum = 0.0;

    while (first1 != last1) {
        prod = *first1 * *first2;
        sum = *first1 + *first2;
        if (prod != 0.0 && sum != 0.0) {
            dist += prod / sum;
        }
        first1++;
        first2++;
    }

    return 2.0 * dist;
}


template <typename InputIt1, typename InputIt2>
double cosine_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double prod = 0.0;
    double p_square = 0.0;
    double q_square = 0.0;
    double dist = 0.0;

    while (first1 != last1) {
        prod = *first1 * *first2;
        p_square += std::pow(*first1, 2.0);
        q_square += std::pow(*first2, 2.0);
        dist += prod;
        first1++;
        first2++;
    }

    double dist2 = std::sqrt(p_square) * std::sqrt(q_square);

    if (dist2 == 0.0) {
        return NAN;
    } else {
        return dist / dist2;
    }
}


template <typename InputIt1, typename InputIt2>
double hassebrook_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double prod = 0.0;
    double p_square = 0.0;
    double q_square = 0.0;
    double dist = 0.0;

    while (first1 != last1) {
        prod = *first1 * *first2;
        p_square += std::pow(*first1, 2.0);
        q_square += std::pow(*first2, 2.0);
        dist += prod;
        first1++;
        first2++;
    }

    double dist2 = p_square + q_square - dist;

    if (dist2 == 0.0) {
        return 0;
    } else {
        return dist / dist2;
    }
}


template <typename InputIt1, typename InputIt2>
double jaccard_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    return 1.0 - hassebrook_internal(first1, last1, first2);
}


template <typename InputIt1, typename InputIt2>
double dice_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double diff_square = 0.0;
    double p_square = 0.0;
    double q_square = 0.0;
    double dist = 0.0;

    while (first1 != last1) {
        diff_square = std::pow((*first1 - *first2), 2.0);
        p_square += std::pow(*first1, 2.0);
        q_square += std::pow(*first2, 2.0);
        dist += diff_square;
        first1++;
        first2++;
    }

    double dist2 = p_square + q_square;

    if (dist2 == 0.0) {
        return NAN;
    } else {
        return dist / dist2;
    }
}


template <typename InputIt1, typename InputIt2>
double fidelity_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;

    while (first1 != last1) {
        dist += std::sqrt(*first1 * *first2);
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double bhattacharyya_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit, double epsilon) {
    double fid_value = fidelity_internal(first1, last1, first2);
    if (fid_value == 0.0) {
        fid_value += epsilon;
    }

    if (unit == "log2")
        return -std::log(fid_value) / std::log(2.0);
    if (unit == "log10")
        return -std::log(fid_value) / std::log(10.0);

    return -std::log(fid_value);
}


template <typename InputIt1, typename InputIt2>
double hellinger_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    return 2.0 * std::sqrt(1.0 - fidelity_internal(first1, last1, first2));
}


template <typename InputIt1, typename InputIt2>
double matusita_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    return std::sqrt(2.0 - (2.0 * fidelity_internal(first1, last1, first2)));
}


template <typename InputIt1, typename InputIt2>
double squared_chord_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;

    while (first1 != last1) {
        dist += std::pow(std::sqrt(*first1) - std::sqrt(*first2), 2.0);
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double squared_euclidean_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;

    while (first1 != last1) {
        dist += std::pow(*first1 - *first2, 2.0);
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double pearson_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, double epsilon) {
    double dist = 0.0;

    while (first1 != last1) {
        if (*first2 == 0.0) {
            dist += std::pow(*first1 - *first2, 2.0) / epsilon;
        } else {
            dist += std::pow(*first1 - *first2, 2.0) / *first2;
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double neyman_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, double epsilon) {
    double dist = 0.0;

    while (first1 != last1) {
        if (*first1 == 0.0) {
            dist += std::pow(*first1 - *first2, 2.0) / epsilon;
        } else {
            dist += std::pow(*first1 - *first2, 2.0) / *first1;
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double squared_chi_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double PQdiff = 0.0;
    double PQsum = 0.0;

    while (first1 != last1) {
        PQdiff = std::pow(*first1 - *first2, 2.0);
        PQsum = *first1 + *first2;
        if (PQdiff != 0.0 && PQsum != 0.0) {
            dist += PQdiff / PQsum;
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double prob_symm_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    return 2.0 * squared_chi_internal(first1, last1, first2);
}


template <typename InputIt1, typename InputIt2>
double divergence_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double PQdiff = 0.0;
    double PQsum = 0.0;

    while (first1 != last1) {
        PQdiff = std::pow(*first1 - *first2, 2.0);
        PQsum = std::pow(*first1 + *first2, 2.0);
        if (PQdiff != 0.0 && PQsum != 0.0) {
            dist += PQdiff / PQsum;
        }
        first1++;
        first2++;
    }

    return 2.0 * dist;
}


template <typename InputIt1, typename InputIt2>
double clark_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double PQdiff = 0.0;
    double PQsum = 0.0;

    while (first1 != last1) {
        PQdiff = std::abs(*first1 - *first2);
        PQsum = *first1 + *first2;
        if (PQdiff != 0.0 && PQsum != 0.0) {
            dist += std::pow(PQdiff / PQsum, 2.0);
        }
        first1++;
        first2++;
    }

    return std::sqrt(dist);
}


template <typename InputIt1, typename InputIt2>
double additive_symm_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double PQsum = 0.0;
    double PQprod = 0.0;

    while (first1 != last1) {
        PQsum = *first1 + *first2;
        PQprod = *first1 * *first2;
        if (PQsum != 0.0 && PQprod != 0.0) {
            dist += std::pow(*first1 - *first2, 2.0) * (PQsum / PQprod);
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double kullback_leibler_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit, double epsilon) {
    double dist = 0.0;
    double PQratio = 0.0;

    while (first1 != last1) {
        if (*first1 != 0.0 || *first2 != 0.0) {
            if (*first2 == 0.0) {
                PQratio = *first1 / epsilon;
            } else {
                PQratio = *first1 / *first2;
            }

            if (PQratio != 0.0) {
                if (unit == "log2") {
                    dist += *first1 * std::log(PQratio) / std::log(2.0);
                } else if (unit == "log10") {
                    dist += *first1 * std::log(PQratio) / std::log(10.0);
                } else {
                    dist += *first1 * std::log(PQratio);
                }
            }
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double jeffreys_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit, double epsilon) {
    double dist = 0.0;
    double PQrate = 0.0;

    while (first1 != last1) {
        if (*first2 == 0.0) {
            PQrate = *first1 / epsilon;
        } else {
            PQrate = *first1 / *first2;
        }

        if (PQrate != 0.0) {
            if (unit == "log2") {
                dist += (*first1 - *first2) * std::log(PQrate) / std::log(2.0);
            } else if (unit == "log10") {
                dist += (*first1 - *first2) * std::log(PQrate) / std::log(10.0);
            } else {
                dist += (*first1 - *first2) * std::log(PQrate);
            }
        } else {
            if (unit == "log2") {
                dist += (*first1 - *first2) * std::log(epsilon) / std::log(2.0);
            } else if (unit == "log10") {
                dist += (*first1 - *first2) * std::log(epsilon) / std::log(10.0);
            } else {
                dist += (*first1 - *first2) * std::log(epsilon);
            }
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double k_divergence_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit) {
    double dist = 0.0;
    while (first1 != last1) {
        double p_i = *first1;
        double q_i = *first2;

        if (p_i > 0.0) {
            double pq_sum = p_i + q_i;
            double term = 0.0;
            if (pq_sum > 0.0) {
                term = p_i * std::log(2.0 * p_i / pq_sum);
            } else {
                // This handles the 0 * log(0/0) case, which should be 0.
            }

            if (unit == "log2")
                dist += term / std::log(2.0);
            else if (unit == "log10")
                dist += term / std::log(10.0);
            else
                dist += term;
        }
        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double topsoe_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit) {
    double dist = 0.0;

    while (first1 != last1) {
        double p_i = *first1;
        double q_i = *first2;
        double pq_sum = p_i + q_i;

        double term1 = 0.0;
        if (p_i > 0.0 && pq_sum != 0.0) {
            term1 = p_i * std::log(2.0 * p_i / pq_sum);
        }

        double term2 = 0.0;
        if (q_i > 0.0 && pq_sum != 0.0) {
            term2 = q_i * std::log(2.0 * q_i / pq_sum);
        }

        double total_term = term1 + term2;

        if (unit == "log2") {
            dist += total_term / std::log(2.0);
        } else if (unit == "log10") {
            dist += total_term / std::log(10.0);
        } else {
            dist += total_term;
        }

        first1++;
        first2++;
    }
    return dist;
}

template <typename InputIt1, typename InputIt2>
double jensen_shannon_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit) {
    double sum1 = 0.0;
    double sum2 = 0.0;
    double PQsum = 0.0;

    while (first1 != last1) {
        PQsum = *first1 + *first2;
        if (unit == "log2") {
            if (*first1 != 0.0 && PQsum != 0.0) {
                sum1 += *first1 * std::log((2.0 * *first1) / PQsum) / std::log(2.0);
            }
            if (*first2 != 0.0 && PQsum != 0.0) {
                sum2 += *first2 * std::log((2.0 * *first2) / PQsum) / std::log(2.0);
            }
        } else if (unit == "log10") {
            if (*first1 != 0.0 && PQsum != 0.0) {
                sum1 += *first1 * std::log((2.0 * *first1) / PQsum) / std::log(10.0);
            }
            if (*first2 != 0.0 && PQsum != 0.0) {
                sum2 += *first2 * std::log((2.0 * *first2) / PQsum) / std::log(10.0);
            }
        } else {
            if (*first1 != 0.0 && PQsum != 0.0) {
                sum1 += *first1 * std::log((2.0 * *first1) / PQsum);
            }
            if (*first2 != 0.0 && PQsum != 0.0) {
                sum2 += *first2 * std::log((2.0 * *first2) / PQsum);
            }
        }
        first1++;
        first2++;
    }

    return 0.5 * (sum1 + sum2);
}


template <typename InputIt1, typename InputIt2>
double jensen_difference_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit) {
    double dist = 0.0;

    while (first1 != last1) {
        double p_i = *first1;
        double q_i = *first2;
        double pq_sum = p_i + q_i;
        double term = 0.0;

        if (pq_sum == 0.0 && p_i == 0.0 && q_i == 0.0) {
            term = 0.0;
        } else if (p_i == 0.0) { // and q_i > 0
            term = (q_i * std::log(q_i)) / 2.0 - (pq_sum / 2.0) * std::log(pq_sum / 2.0);
        } else if (q_i == 0.0) { // and p_i > 0
            term = (p_i * std::log(p_i)) / 2.0 - (pq_sum / 2.0) * std::log(pq_sum / 2.0);
        } else { // p_i > 0 and q_i > 0
            term = ((p_i * std::log(p_i)) + (q_i * std::log(q_i))) / 2.0 - (pq_sum / 2.0) * std::log(pq_sum / 2.0);
        }

        if (unit == "log2") {
            dist += term / std::log(2.0);
        } else if (unit == "log10") {
            dist += term / std::log(10.0);
        } else { // "log"
            dist += term;
        }

        first1++;
        first2++;
    }

    return dist;
}


template <typename InputIt1, typename InputIt2>
double taneja_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, const std::string& unit, double epsilon) {
    double dist = 0.0;
    double PQsum = 0.0;
    double denominator = 0.0;

    while (first1 != last1) {
        PQsum = *first1 + *first2;
        if (PQsum != 0.0) {
            denominator = (2.0 * std::sqrt(*first1 * *first2));
            if (denominator == 0.0) {
                if (unit == "log2") {
                    dist += (PQsum / 2.0) * std::log(PQsum / epsilon) / std::log(2.0);
                } else if (unit == "log10") {
                    dist += (PQsum / 2.0) * std::log(PQsum / epsilon) / std::log(10.0);
                } else {
                    dist += (PQsum / 2.0) * std::log(PQsum / epsilon);
                }
            } else {
                if (unit == "log2") {
                    dist += (PQsum / 2.0) * std::log(PQsum / denominator) / std::log(2.0);
                } else if (unit == "log10") {
                    dist += (PQsum / 2.0) * std::log(PQsum / denominator) / std::log(10.0);
                } else {
                    dist += (PQsum / 2.0) * std::log(PQsum / denominator);
                }
            }
        }
        first1++;
        first2++;
    }

    return dist;
}

// helper function to compute the canberra distance
template <typename InputIt1, typename InputIt2>
double canberra_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double diff = 0.0;
    double sum = 0.0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        sum = std::abs(*first1) + std::abs(*first2);
        if (sum != 0.0) {
            dist += diff / sum;
        }
        first1++;
        first2++;
    }
    return dist;
}

// helper function to compute the intersection distance
template <typename InputIt1, typename InputIt2>
double intersection_dist_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;

    while (first1 != last1) {
        dist += std::min(*first1, *first2);
        first1++;
        first2++;
    }
    return dist;
}

// helper function to compute the wave hedges distance
template <typename InputIt1, typename InputIt2>
double wave_hedges_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double diff = 0.0;
    double max_val = 0.0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        max_val = std::max(*first1, *first2);
        if (max_val != 0.0) {
            dist += diff / max_val;
        }
        first1++;
        first2++;
    }
    return dist;
}

// helper function to compute the czekanowski distance
template <typename InputIt1, typename InputIt2>
double czekanowski_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double sum_min = 0.0;
    double sum_plus = 0.0;

    while (first1 != last1) {
        sum_min += std::min(*first1, *first2);
        sum_plus += (*first1 + *first2);
        first1++;
        first2++;
    }
    
    if (sum_plus == 0.0) {
        return 0.0;
    }
    
    return 1.0 - ((2.0 * sum_min) / sum_plus);
}


// helper function to compute the motyka distance
template <typename InputIt1, typename InputIt2>
double motyka_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double sum_min = 0.0;
    double sum_plus = 0.0;

    while (first1 != last1) {
        sum_min += std::min(*first1, *first2);
        sum_plus += (*first1 + *first2);
        first1++;
        first2++;
    }
    
    if (sum_plus == 0.0) {
        return 0.0;
    }
    
    return 1 - (sum_min / sum_plus);
}

template <typename InputIt1, typename InputIt2>
double kumar_johnson_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2, double epsilon) {
    double dist = 0.0;
    double divisor = 0.0;

    while (first1 != last1) {
        double p_i = *first1;
        double q_i = *first2;

        divisor = (2.0 * pow(p_i * q_i, 1.5));

        if (divisor == 0.0) {
            dist += pow(pow(p_i, 2.0) - pow(q_i, 2.0), 2.0) / epsilon;
        } else {
            dist += pow(pow(p_i, 2.0) - pow(q_i, 2.0), 2.0) / divisor;
        }

        ++first1;
        ++first2;
    }
    return dist;
}

template <typename InputIt1, typename InputIt2>
double avg_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double pq_diff = 0.0;
    double pq_max = 0.0;

    while (first1 != last1) {
        double p_i = *first1;
        double q_i = *first2;

        pq_diff = fabs(p_i - q_i);

        if (pq_diff > pq_max)
            pq_max = pq_diff;

        dist += pq_diff;

        ++first1;
        ++first2;
    }

    return (dist + pq_max) / 2.0;
}

template <typename InputIt1, typename InputIt2>
double chebyshev_internal(InputIt1 first1, InputIt1 last1, InputIt2 first2) {
    double dist = 0.0;
    double diff = 0.0;

    while (first1 != last1) {
        diff = std::abs(*first1 - *first2);
        if (diff > dist) {
            dist = diff;
        }
        first1++;
        first2++;
    }
    return dist;
}

#endif // philentropy_Distances_Internal_H
