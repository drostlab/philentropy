#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP philentropy_additive_symm_chi_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_as_data_frame(SEXP);
extern SEXP philentropy_as_matrix(SEXP);
extern SEXP philentropy_avg(SEXP, SEXP, SEXP);
extern SEXP philentropy_bhattacharyya(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_canberra(SEXP, SEXP, SEXP);
extern SEXP philentropy_CEcpp(SEXP, SEXP, SEXP);
extern SEXP philentropy_chebyshev(SEXP, SEXP, SEXP);
extern SEXP philentropy_clark_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_cosine_dist(SEXP, SEXP, SEXP);
extern SEXP philentropy_custom_log10(SEXP);
extern SEXP philentropy_custom_log2(SEXP);
extern SEXP philentropy_czekanowski(SEXP, SEXP, SEXP);
extern SEXP philentropy_dice_dist(SEXP, SEXP, SEXP);
extern SEXP philentropy_DistMatrixMinkowskiMAT(SEXP, SEXP, SEXP);
extern SEXP philentropy_DistMatrixWithoutUnitDF(SEXP, SEXP, SEXP);
extern SEXP philentropy_DistMatrixWithoutUnitMAT(SEXP, SEXP, SEXP);
extern SEXP philentropy_DistMatrixWithUnitDF(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_DistMatrixWithUnitMAT(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_divergence_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_Ecpp(SEXP, SEXP);
extern SEXP philentropy_est_prob_empirical(SEXP);
extern SEXP philentropy_euclidean(SEXP, SEXP, SEXP);
extern SEXP philentropy_fidelity(SEXP, SEXP, SEXP);
extern SEXP philentropy_gower(SEXP, SEXP, SEXP);
extern SEXP philentropy_harmonic_mean_dist(SEXP, SEXP, SEXP);
extern SEXP philentropy_hellinger(SEXP, SEXP, SEXP);
extern SEXP philentropy_inner_product(SEXP, SEXP, SEXP);
extern SEXP philentropy_intersection_dist(SEXP, SEXP, SEXP);
extern SEXP philentropy_jaccard(SEXP, SEXP, SEXP);
extern SEXP philentropy_JEcpp(SEXP, SEXP);
extern SEXP philentropy_jeffreys(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_jensen_difference(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_jensen_shannon(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_k_divergence(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_kulczynski_d(SEXP, SEXP, SEXP);
extern SEXP philentropy_kullback_leibler_distance(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_kumar_hassebrook(SEXP, SEXP, SEXP);
extern SEXP philentropy_kumar_johnson(SEXP, SEXP, SEXP);
extern SEXP philentropy_lorentzian(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_manhattan(SEXP, SEXP, SEXP);
extern SEXP philentropy_matusita(SEXP, SEXP, SEXP);
extern SEXP philentropy_MIcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_minkowski(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_motyka(SEXP, SEXP, SEXP);
extern SEXP philentropy_neyman_chi_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_pearson_chi_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_pearson_corr_centred(SEXP, SEXP, SEXP);
extern SEXP philentropy_pearson_corr_uncentred(SEXP, SEXP, SEXP);
extern SEXP philentropy_prob_symm_chi_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_RcppExport_registerCCallable();
extern SEXP philentropy_ruzicka(SEXP, SEXP, SEXP);
extern SEXP philentropy_soergel(SEXP, SEXP, SEXP);
extern SEXP philentropy_sorensen(SEXP, SEXP, SEXP);
extern SEXP philentropy_squared_chi_sq(SEXP, SEXP, SEXP);
extern SEXP philentropy_squared_chord(SEXP, SEXP, SEXP);
extern SEXP philentropy_squared_euclidean(SEXP, SEXP, SEXP);
extern SEXP philentropy_squared_pearson_corr(SEXP, SEXP, SEXP);
extern SEXP philentropy_sum_rcpp(SEXP);
extern SEXP philentropy_taneja(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_tanimoto(SEXP, SEXP, SEXP);
extern SEXP philentropy_topsoe(SEXP, SEXP, SEXP, SEXP);
extern SEXP philentropy_wave_hedges(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"philentropy_additive_symm_chi_sq",         (DL_FUNC) &philentropy_additive_symm_chi_sq,         3},
    {"philentropy_as_data_frame",                (DL_FUNC) &philentropy_as_data_frame,                1},
    {"philentropy_as_matrix",                    (DL_FUNC) &philentropy_as_matrix,                    1},
    {"philentropy_avg",                          (DL_FUNC) &philentropy_avg,                          3},
    {"philentropy_bhattacharyya",                (DL_FUNC) &philentropy_bhattacharyya,                4},
    {"philentropy_canberra",                     (DL_FUNC) &philentropy_canberra,                     3},
    {"philentropy_CEcpp",                        (DL_FUNC) &philentropy_CEcpp,                        3},
    {"philentropy_chebyshev",                    (DL_FUNC) &philentropy_chebyshev,                    3},
    {"philentropy_clark_sq",                     (DL_FUNC) &philentropy_clark_sq,                     3},
    {"philentropy_cosine_dist",                  (DL_FUNC) &philentropy_cosine_dist,                  3},
    {"philentropy_custom_log10",                 (DL_FUNC) &philentropy_custom_log10,                 1},
    {"philentropy_custom_log2",                  (DL_FUNC) &philentropy_custom_log2,                  1},
    {"philentropy_czekanowski",                  (DL_FUNC) &philentropy_czekanowski,                  3},
    {"philentropy_dice_dist",                    (DL_FUNC) &philentropy_dice_dist,                    3},
    {"philentropy_DistMatrixMinkowskiMAT",       (DL_FUNC) &philentropy_DistMatrixMinkowskiMAT,       3},
    {"philentropy_DistMatrixWithoutUnitDF",      (DL_FUNC) &philentropy_DistMatrixWithoutUnitDF,      3},
    {"philentropy_DistMatrixWithoutUnitMAT",     (DL_FUNC) &philentropy_DistMatrixWithoutUnitMAT,     3},
    {"philentropy_DistMatrixWithUnitDF",         (DL_FUNC) &philentropy_DistMatrixWithUnitDF,         4},
    {"philentropy_DistMatrixWithUnitMAT",        (DL_FUNC) &philentropy_DistMatrixWithUnitMAT,        4},
    {"philentropy_divergence_sq",                (DL_FUNC) &philentropy_divergence_sq,                3},
    {"philentropy_Ecpp",                         (DL_FUNC) &philentropy_Ecpp,                         2},
    {"philentropy_est_prob_empirical",           (DL_FUNC) &philentropy_est_prob_empirical,           1},
    {"philentropy_euclidean",                    (DL_FUNC) &philentropy_euclidean,                    3},
    {"philentropy_fidelity",                     (DL_FUNC) &philentropy_fidelity,                     3},
    {"philentropy_gower",                        (DL_FUNC) &philentropy_gower,                        3},
    {"philentropy_harmonic_mean_dist",           (DL_FUNC) &philentropy_harmonic_mean_dist,           3},
    {"philentropy_hellinger",                    (DL_FUNC) &philentropy_hellinger,                    3},
    {"philentropy_inner_product",                (DL_FUNC) &philentropy_inner_product,                3},
    {"philentropy_intersection_dist",            (DL_FUNC) &philentropy_intersection_dist,            3},
    {"philentropy_jaccard",                      (DL_FUNC) &philentropy_jaccard,                      3},
    {"philentropy_JEcpp",                        (DL_FUNC) &philentropy_JEcpp,                        2},
    {"philentropy_jeffreys",                     (DL_FUNC) &philentropy_jeffreys,                     4},
    {"philentropy_jensen_difference",            (DL_FUNC) &philentropy_jensen_difference,            4},
    {"philentropy_jensen_shannon",               (DL_FUNC) &philentropy_jensen_shannon,               4},
    {"philentropy_k_divergence",                 (DL_FUNC) &philentropy_k_divergence,                 4},
    {"philentropy_kulczynski_d",                 (DL_FUNC) &philentropy_kulczynski_d,                 3},
    {"philentropy_kullback_leibler_distance",    (DL_FUNC) &philentropy_kullback_leibler_distance,    4},
    {"philentropy_kumar_hassebrook",             (DL_FUNC) &philentropy_kumar_hassebrook,             3},
    {"philentropy_kumar_johnson",                (DL_FUNC) &philentropy_kumar_johnson,                3},
    {"philentropy_lorentzian",                   (DL_FUNC) &philentropy_lorentzian,                   4},
    {"philentropy_manhattan",                    (DL_FUNC) &philentropy_manhattan,                    3},
    {"philentropy_matusita",                     (DL_FUNC) &philentropy_matusita,                     3},
    {"philentropy_MIcpp",                        (DL_FUNC) &philentropy_MIcpp,                        4},
    {"philentropy_minkowski",                    (DL_FUNC) &philentropy_minkowski,                    4},
    {"philentropy_motyka",                       (DL_FUNC) &philentropy_motyka,                       3},
    {"philentropy_neyman_chi_sq",                (DL_FUNC) &philentropy_neyman_chi_sq,                3},
    {"philentropy_pearson_chi_sq",               (DL_FUNC) &philentropy_pearson_chi_sq,               3},
    {"philentropy_pearson_corr_centred",         (DL_FUNC) &philentropy_pearson_corr_centred,         3},
    {"philentropy_pearson_corr_uncentred",       (DL_FUNC) &philentropy_pearson_corr_uncentred,       3},
    {"philentropy_prob_symm_chi_sq",             (DL_FUNC) &philentropy_prob_symm_chi_sq,             3},
    {"philentropy_RcppExport_registerCCallable", (DL_FUNC) &philentropy_RcppExport_registerCCallable, 0},
    {"philentropy_ruzicka",                      (DL_FUNC) &philentropy_ruzicka,                      3},
    {"philentropy_soergel",                      (DL_FUNC) &philentropy_soergel,                      3},
    {"philentropy_sorensen",                     (DL_FUNC) &philentropy_sorensen,                     3},
    {"philentropy_squared_chi_sq",               (DL_FUNC) &philentropy_squared_chi_sq,               3},
    {"philentropy_squared_chord",                (DL_FUNC) &philentropy_squared_chord,                3},
    {"philentropy_squared_euclidean",            (DL_FUNC) &philentropy_squared_euclidean,            3},
    {"philentropy_squared_pearson_corr",         (DL_FUNC) &philentropy_squared_pearson_corr,         3},
    {"philentropy_sum_rcpp",                     (DL_FUNC) &philentropy_sum_rcpp,                     1},
    {"philentropy_taneja",                       (DL_FUNC) &philentropy_taneja,                       4},
    {"philentropy_tanimoto",                     (DL_FUNC) &philentropy_tanimoto,                     3},
    {"philentropy_topsoe",                       (DL_FUNC) &philentropy_topsoe,                       4},
    {"philentropy_wave_hedges",                  (DL_FUNC) &philentropy_wave_hedges,                  3},
    {NULL, NULL, 0}
};

void R_init_philentropy(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}