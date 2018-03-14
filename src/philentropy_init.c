#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _philentropy_additive_symm_chi_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_as_data_frame(SEXP);
extern SEXP _philentropy_as_matrix(SEXP);
extern SEXP _philentropy_avg(SEXP, SEXP, SEXP);
extern SEXP _philentropy_bhattacharyya(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_canberra(SEXP, SEXP, SEXP);
extern SEXP _philentropy_CEcpp(SEXP, SEXP, SEXP);
extern SEXP _philentropy_chebyshev(SEXP, SEXP, SEXP);
extern SEXP _philentropy_clark_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_cosine_dist(SEXP, SEXP, SEXP);
extern SEXP _philentropy_custom_log10(SEXP);
extern SEXP _philentropy_custom_log2(SEXP);
extern SEXP _philentropy_czekanowski(SEXP, SEXP, SEXP);
extern SEXP _philentropy_dice_dist(SEXP, SEXP, SEXP);
extern SEXP _philentropy_DistMatrixMinkowskiMAT(SEXP, SEXP, SEXP);
extern SEXP _philentropy_DistMatrixWithoutUnitDF(SEXP, SEXP, SEXP);
extern SEXP _philentropy_DistMatrixWithoutUnitMAT(SEXP, SEXP, SEXP);
extern SEXP _philentropy_DistMatrixWithUnitDF(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_DistMatrixWithUnitMAT(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_divergence_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_Ecpp(SEXP, SEXP);
extern SEXP _philentropy_est_prob_empirical(SEXP);
extern SEXP _philentropy_euclidean(SEXP, SEXP, SEXP);
extern SEXP _philentropy_fidelity(SEXP, SEXP, SEXP);
extern SEXP _philentropy_gower(SEXP, SEXP, SEXP);
extern SEXP _philentropy_harmonic_mean_dist(SEXP, SEXP, SEXP);
extern SEXP _philentropy_hellinger(SEXP, SEXP, SEXP);
extern SEXP _philentropy_inner_product(SEXP, SEXP, SEXP);
extern SEXP _philentropy_intersection_dist(SEXP, SEXP, SEXP);
extern SEXP _philentropy_jaccard(SEXP, SEXP, SEXP);
extern SEXP _philentropy_JEcpp(SEXP, SEXP);
extern SEXP _philentropy_jeffreys(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_jensen_difference(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_jensen_shannon(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_k_divergence(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_kulczynski_d(SEXP, SEXP, SEXP);
extern SEXP _philentropy_kullback_leibler_distance(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_kumar_hassebrook(SEXP, SEXP, SEXP);
extern SEXP _philentropy_kumar_johnson(SEXP, SEXP, SEXP);
extern SEXP _philentropy_lorentzian(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_manhattan(SEXP, SEXP, SEXP);
extern SEXP _philentropy_matusita(SEXP, SEXP, SEXP);
extern SEXP _philentropy_MIcpp(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_minkowski(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_motyka(SEXP, SEXP, SEXP);
extern SEXP _philentropy_neyman_chi_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_pearson_chi_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_pearson_corr_centred(SEXP, SEXP, SEXP);
extern SEXP _philentropy_pearson_corr_uncentred(SEXP, SEXP, SEXP);
extern SEXP _philentropy_prob_symm_chi_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_RcppExport_registerCCallable();
extern SEXP _philentropy_ruzicka(SEXP, SEXP, SEXP);
extern SEXP _philentropy_soergel(SEXP, SEXP, SEXP);
extern SEXP _philentropy_sorensen(SEXP, SEXP, SEXP);
extern SEXP _philentropy_squared_chi_sq(SEXP, SEXP, SEXP);
extern SEXP _philentropy_squared_chord(SEXP, SEXP, SEXP);
extern SEXP _philentropy_squared_euclidean(SEXP, SEXP, SEXP);
extern SEXP _philentropy_squared_pearson_corr(SEXP, SEXP, SEXP);
extern SEXP _philentropy_sum_rcpp(SEXP);
extern SEXP _philentropy_taneja(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_tanimoto(SEXP, SEXP, SEXP);
extern SEXP _philentropy_topsoe(SEXP, SEXP, SEXP, SEXP);
extern SEXP _philentropy_wave_hedges(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_philentropy_additive_symm_chi_sq",         (DL_FUNC) &_philentropy_additive_symm_chi_sq,         3},
  {"_philentropy_as_data_frame",                (DL_FUNC) &_philentropy_as_data_frame,                1},
  {"_philentropy_as_matrix",                    (DL_FUNC) &_philentropy_as_matrix,                    1},
  {"_philentropy_avg",                          (DL_FUNC) &_philentropy_avg,                          3},
  {"_philentropy_bhattacharyya",                (DL_FUNC) &_philentropy_bhattacharyya,                4},
  {"_philentropy_canberra",                     (DL_FUNC) &_philentropy_canberra,                     3},
  {"_philentropy_CEcpp",                        (DL_FUNC) &_philentropy_CEcpp,                        3},
  {"_philentropy_chebyshev",                    (DL_FUNC) &_philentropy_chebyshev,                    3},
  {"_philentropy_clark_sq",                     (DL_FUNC) &_philentropy_clark_sq,                     3},
  {"_philentropy_cosine_dist",                  (DL_FUNC) &_philentropy_cosine_dist,                  3},
  {"_philentropy_custom_log10",                 (DL_FUNC) &_philentropy_custom_log10,                 1},
  {"_philentropy_custom_log2",                  (DL_FUNC) &_philentropy_custom_log2,                  1},
  {"_philentropy_czekanowski",                  (DL_FUNC) &_philentropy_czekanowski,                  3},
  {"_philentropy_dice_dist",                    (DL_FUNC) &_philentropy_dice_dist,                    3},
  {"_philentropy_DistMatrixMinkowskiMAT",       (DL_FUNC) &_philentropy_DistMatrixMinkowskiMAT,       3},
  {"_philentropy_DistMatrixWithoutUnitDF",      (DL_FUNC) &_philentropy_DistMatrixWithoutUnitDF,      3},
  {"_philentropy_DistMatrixWithoutUnitMAT",     (DL_FUNC) &_philentropy_DistMatrixWithoutUnitMAT,     3},
  {"_philentropy_DistMatrixWithUnitDF",         (DL_FUNC) &_philentropy_DistMatrixWithUnitDF,         4},
  {"_philentropy_DistMatrixWithUnitMAT",        (DL_FUNC) &_philentropy_DistMatrixWithUnitMAT,        4},
  {"_philentropy_divergence_sq",                (DL_FUNC) &_philentropy_divergence_sq,                3},
  {"_philentropy_Ecpp",                         (DL_FUNC) &_philentropy_Ecpp,                         2},
  {"_philentropy_est_prob_empirical",           (DL_FUNC) &_philentropy_est_prob_empirical,           1},
  {"_philentropy_euclidean",                    (DL_FUNC) &_philentropy_euclidean,                    3},
  {"_philentropy_fidelity",                     (DL_FUNC) &_philentropy_fidelity,                     3},
  {"_philentropy_gower",                        (DL_FUNC) &_philentropy_gower,                        3},
  {"_philentropy_harmonic_mean_dist",           (DL_FUNC) &_philentropy_harmonic_mean_dist,           3},
  {"_philentropy_hellinger",                    (DL_FUNC) &_philentropy_hellinger,                    3},
  {"_philentropy_inner_product",                (DL_FUNC) &_philentropy_inner_product,                3},
  {"_philentropy_intersection_dist",            (DL_FUNC) &_philentropy_intersection_dist,            3},
  {"_philentropy_jaccard",                      (DL_FUNC) &_philentropy_jaccard,                      3},
  {"_philentropy_JEcpp",                        (DL_FUNC) &_philentropy_JEcpp,                        2},
  {"_philentropy_jeffreys",                     (DL_FUNC) &_philentropy_jeffreys,                     4},
  {"_philentropy_jensen_difference",            (DL_FUNC) &_philentropy_jensen_difference,            4},
  {"_philentropy_jensen_shannon",               (DL_FUNC) &_philentropy_jensen_shannon,               4},
  {"_philentropy_k_divergence",                 (DL_FUNC) &_philentropy_k_divergence,                 4},
  {"_philentropy_kulczynski_d",                 (DL_FUNC) &_philentropy_kulczynski_d,                 3},
  {"_philentropy_kullback_leibler_distance",    (DL_FUNC) &_philentropy_kullback_leibler_distance,    4},
  {"_philentropy_kumar_hassebrook",             (DL_FUNC) &_philentropy_kumar_hassebrook,             3},
  {"_philentropy_kumar_johnson",                (DL_FUNC) &_philentropy_kumar_johnson,                3},
  {"_philentropy_lorentzian",                   (DL_FUNC) &_philentropy_lorentzian,                   4},
  {"_philentropy_manhattan",                    (DL_FUNC) &_philentropy_manhattan,                    3},
  {"_philentropy_matusita",                     (DL_FUNC) &_philentropy_matusita,                     3},
  {"_philentropy_MIcpp",                        (DL_FUNC) &_philentropy_MIcpp,                        4},
  {"_philentropy_minkowski",                    (DL_FUNC) &_philentropy_minkowski,                    4},
  {"_philentropy_motyka",                       (DL_FUNC) &_philentropy_motyka,                       3},
  {"_philentropy_neyman_chi_sq",                (DL_FUNC) &_philentropy_neyman_chi_sq,                3},
  {"_philentropy_pearson_chi_sq",               (DL_FUNC) &_philentropy_pearson_chi_sq,               3},
  {"_philentropy_pearson_corr_centred",         (DL_FUNC) &_philentropy_pearson_corr_centred,         3},
  {"_philentropy_pearson_corr_uncentred",       (DL_FUNC) &_philentropy_pearson_corr_uncentred,       3},
  {"_philentropy_prob_symm_chi_sq",             (DL_FUNC) &_philentropy_prob_symm_chi_sq,             3},
  {"_philentropy_RcppExport_registerCCallable", (DL_FUNC) &_philentropy_RcppExport_registerCCallable, 0},
  {"_philentropy_ruzicka",                      (DL_FUNC) &_philentropy_ruzicka,                      3},
  {"_philentropy_soergel",                      (DL_FUNC) &_philentropy_soergel,                      3},
  {"_philentropy_sorensen",                     (DL_FUNC) &_philentropy_sorensen,                     3},
  {"_philentropy_squared_chi_sq",               (DL_FUNC) &_philentropy_squared_chi_sq,               3},
  {"_philentropy_squared_chord",                (DL_FUNC) &_philentropy_squared_chord,                3},
  {"_philentropy_squared_euclidean",            (DL_FUNC) &_philentropy_squared_euclidean,            3},
  {"_philentropy_squared_pearson_corr",         (DL_FUNC) &_philentropy_squared_pearson_corr,         3},
  {"_philentropy_sum_rcpp",                     (DL_FUNC) &_philentropy_sum_rcpp,                     1},
  {"_philentropy_taneja",                       (DL_FUNC) &_philentropy_taneja,                       4},
  {"_philentropy_tanimoto",                     (DL_FUNC) &_philentropy_tanimoto,                     3},
  {"_philentropy_topsoe",                       (DL_FUNC) &_philentropy_topsoe,                       4},
  {"_philentropy_wave_hedges",                  (DL_FUNC) &_philentropy_wave_hedges,                  3},
  {NULL, NULL, 0}
};

void R_init_philentropy(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
