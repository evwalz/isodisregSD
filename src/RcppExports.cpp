// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// new_func
NumericMatrix new_func(NumericMatrix X, double eps);
RcppExport SEXP _isodisregAFSD_new_func(SEXP XSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func(X, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_sd
NumericMatrix new_func_sd(NumericMatrix X);
RcppExport SEXP _isodisregAFSD_new_func_sd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_sd(X));
    return rcpp_result_gen;
END_RCPP
}
// new_func_eps
List new_func_eps(NumericMatrix X);
RcppExport SEXP _isodisregAFSD_new_func_eps(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_eps(X));
    return rcpp_result_gen;
END_RCPP
}
// normal_comp
NumericMatrix normal_comp(NumericMatrix X, double eps);
RcppExport SEXP _isodisregAFSD_normal_comp(SEXP XSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_comp(X, eps));
    return rcpp_result_gen;
END_RCPP
}
// normal_comp_eps
List normal_comp_eps(NumericMatrix X);
RcppExport SEXP _isodisregAFSD_normal_comp_eps(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_comp_eps(X));
    return rcpp_result_gen;
END_RCPP
}
// normal_comp_sd
NumericMatrix normal_comp_sd(NumericMatrix X);
RcppExport SEXP _isodisregAFSD_normal_comp_sd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(normal_comp_sd(X));
    return rcpp_result_gen;
END_RCPP
}
// indx_norm
List indx_norm(NumericMatrix X, NumericMatrix x, double eps);
RcppExport SEXP _isodisregAFSD_indx_norm(SEXP XSEXP, SEXP xSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(indx_norm(X, x, eps));
    return rcpp_result_gen;
END_RCPP
}
// indx_norm_sd
List indx_norm_sd(NumericMatrix X, NumericMatrix x);
RcppExport SEXP _isodisregAFSD_indx_norm_sd(SEXP XSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(indx_norm_sd(X, x));
    return rcpp_result_gen;
END_RCPP
}
// new_func2
List new_func2(NumericMatrix X, NumericMatrix x, double eps);
RcppExport SEXP _isodisregAFSD_new_func2(SEXP XSEXP, SEXP xSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func2(X, x, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func2_sd
List new_func2_sd(NumericMatrix X, NumericMatrix x);
RcppExport SEXP _isodisregAFSD_new_func2_sd(SEXP XSEXP, SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func2_sd(X, x));
    return rcpp_result_gen;
END_RCPP
}
// new_func_single_grid
List new_func_single_grid(NumericMatrix X, NumericMatrix x, NumericVector vec, double eps);
RcppExport SEXP _isodisregAFSD_new_func_single_grid(SEXP XSEXP, SEXP xSEXP, SEXP vecSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_single_grid(X, x, vec, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_single_grid_sd
List new_func_single_grid_sd(NumericMatrix X, NumericMatrix x, NumericVector vec);
RcppExport SEXP _isodisregAFSD_new_func_single_grid_sd(SEXP XSEXP, SEXP xSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_single_grid_sd(X, x, vec));
    return rcpp_result_gen;
END_RCPP
}
// new_func_single_grid_rememeber
List new_func_single_grid_rememeber(NumericMatrix X, NumericMatrix x, NumericVector vec, double eps);
RcppExport SEXP _isodisregAFSD_new_func_single_grid_rememeber(SEXP XSEXP, SEXP xSEXP, SEXP vecSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_single_grid_rememeber(X, x, vec, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_mat
List new_func_mat(NumericMatrix X, NumericMatrix x, NumericVector gridx, NumericVector gridX, double eps);
RcppExport SEXP _isodisregAFSD_new_func_mat(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridX(gridXSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_mat(X, x, gridx, gridX, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_mat_sd
List new_func_mat_sd(NumericMatrix X, NumericMatrix x, NumericVector gridx, NumericVector gridX);
RcppExport SEXP _isodisregAFSD_new_func_mat_sd(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridX(gridXSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_mat_sd(X, x, gridx, gridX));
    return rcpp_result_gen;
END_RCPP
}
// new_func_mat_list
List new_func_mat_list(List X, NumericMatrix x, NumericVector gridx, List gridX, double eps);
RcppExport SEXP _isodisregAFSD_new_func_mat_list(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< List >::type gridX(gridXSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_mat_list(X, x, gridx, gridX, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_mat_list_sd
List new_func_mat_list_sd(List X, NumericMatrix x, NumericVector gridx, List gridX);
RcppExport SEXP _isodisregAFSD_new_func_mat_list_sd(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< List >::type gridX(gridXSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_mat_list_sd(X, x, gridx, gridX));
    return rcpp_result_gen;
END_RCPP
}
// new_func_list
List new_func_list(List X, List x, List gridx, List gridX, double eps);
RcppExport SEXP _isodisregAFSD_new_func_list(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< List >::type gridX(gridXSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_list(X, x, gridx, gridX, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_list_back
List new_func_list_back(List X, List x, List gridx, List gridX, double eps);
RcppExport SEXP _isodisregAFSD_new_func_list_back(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< List >::type gridX(gridXSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_list_back(X, x, gridx, gridX, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_list_sd_back
List new_func_list_sd_back(List X, List x, List gridx, List gridX);
RcppExport SEXP _isodisregAFSD_new_func_list_sd_back(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< List >::type gridX(gridXSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_list_sd_back(X, x, gridx, gridX));
    return rcpp_result_gen;
END_RCPP
}
// new_func_list_sd
List new_func_list_sd(List X, List x, List gridx, List gridX);
RcppExport SEXP _isodisregAFSD_new_func_list_sd(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< List >::type gridX(gridXSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_list_sd(X, x, gridx, gridX));
    return rcpp_result_gen;
END_RCPP
}
// new_func_list_mat
List new_func_list_mat(NumericMatrix X, List x, List gridx, NumericVector gridX, double eps);
RcppExport SEXP _isodisregAFSD_new_func_list_mat(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridX(gridXSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_list_mat(X, x, gridx, gridX, eps));
    return rcpp_result_gen;
END_RCPP
}
// new_func_list_mat_sd
List new_func_list_mat_sd(NumericMatrix X, List x, List gridx, NumericVector gridX);
RcppExport SEXP _isodisregAFSD_new_func_list_mat_sd(SEXP XSEXP, SEXP xSEXP, SEXP gridxSEXP, SEXP gridXSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< List >::type gridx(gridxSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type gridX(gridXSEXP);
    rcpp_result_gen = Rcpp::wrap(new_func_list_mat_sd(X, x, gridx, gridX));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_list_comp_class
List ecdf_list_comp_class(List X, List t, double eps);
RcppExport SEXP _isodisregAFSD_ecdf_list_comp_class(SEXP XSEXP, SEXP tSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_list_comp_class(X, t, eps));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_list_comp_class_sd
List ecdf_list_comp_class_sd(List X, List t);
RcppExport SEXP _isodisregAFSD_ecdf_list_comp_class_sd(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_list_comp_class_sd(X, t));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_list_comp_class_eps
List ecdf_list_comp_class_eps(List X, List t);
RcppExport SEXP _isodisregAFSD_ecdf_list_comp_class_eps(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_list_comp_class_eps(X, t));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp
NumericMatrix ecdf_comp(NumericMatrix X, NumericVector t, double eps);
RcppExport SEXP _isodisregAFSD_ecdf_comp(SEXP XSEXP, SEXP tSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp(X, t, eps));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp_sd
NumericMatrix ecdf_comp_sd(NumericMatrix X, NumericVector t);
RcppExport SEXP _isodisregAFSD_ecdf_comp_sd(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp_sd(X, t));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp_eps
List ecdf_comp_eps(NumericMatrix X, NumericVector t);
RcppExport SEXP _isodisregAFSD_ecdf_comp_eps(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp_eps(X, t));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp_class
List ecdf_comp_class(NumericMatrix X, NumericVector t, double eps);
RcppExport SEXP _isodisregAFSD_ecdf_comp_class(SEXP XSEXP, SEXP tSEXP, SEXP epsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp_class(X, t, eps));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp_class_sd
List ecdf_comp_class_sd(NumericMatrix X, NumericVector t);
RcppExport SEXP _isodisregAFSD_ecdf_comp_class_sd(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp_class_sd(X, t));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_func
NumericMatrix ecdf_func(List X, List t, NumericVector thresholds);
RcppExport SEXP _isodisregAFSD_ecdf_func(SEXP XSEXP, SEXP tSEXP, SEXP thresholdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresholds(thresholdsSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_func(X, t, thresholds));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_func_mat
NumericMatrix ecdf_func_mat(NumericMatrix X, NumericVector t, NumericVector thresholds);
RcppExport SEXP _isodisregAFSD_ecdf_func_mat(SEXP XSEXP, SEXP tSEXP, SEXP thresholdsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresholds(thresholdsSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_func_mat(X, t, thresholds));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp_class_ind
NumericVector ecdf_comp_class_ind(NumericMatrix X);
RcppExport SEXP _isodisregAFSD_ecdf_comp_class_ind(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp_class_ind(X));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_list_comp_class_ind
NumericVector ecdf_list_comp_class_ind(List X, List t);
RcppExport SEXP _isodisregAFSD_ecdf_list_comp_class_ind(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type X(XSEXP);
    Rcpp::traits::input_parameter< List >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_list_comp_class_ind(X, t));
    return rcpp_result_gen;
END_RCPP
}
// ecdf_comp_class_eps
List ecdf_comp_class_eps(NumericMatrix X, NumericVector t);
RcppExport SEXP _isodisregAFSD_ecdf_comp_class_eps(SEXP XSEXP, SEXP tSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type t(tSEXP);
    rcpp_result_gen = Rcpp::wrap(ecdf_comp_class_eps(X, t));
    return rcpp_result_gen;
END_RCPP
}
// pavaCorrect
NumericMatrix pavaCorrect(NumericMatrix y);
RcppExport SEXP _isodisregAFSD_pavaCorrect(SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(pavaCorrect(y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_isodisregAFSD_new_func", (DL_FUNC) &_isodisregAFSD_new_func, 2},
    {"_isodisregAFSD_new_func_sd", (DL_FUNC) &_isodisregAFSD_new_func_sd, 1},
    {"_isodisregAFSD_new_func_eps", (DL_FUNC) &_isodisregAFSD_new_func_eps, 1},
    {"_isodisregAFSD_normal_comp", (DL_FUNC) &_isodisregAFSD_normal_comp, 2},
    {"_isodisregAFSD_normal_comp_eps", (DL_FUNC) &_isodisregAFSD_normal_comp_eps, 1},
    {"_isodisregAFSD_normal_comp_sd", (DL_FUNC) &_isodisregAFSD_normal_comp_sd, 1},
    {"_isodisregAFSD_indx_norm", (DL_FUNC) &_isodisregAFSD_indx_norm, 3},
    {"_isodisregAFSD_indx_norm_sd", (DL_FUNC) &_isodisregAFSD_indx_norm_sd, 2},
    {"_isodisregAFSD_new_func2", (DL_FUNC) &_isodisregAFSD_new_func2, 3},
    {"_isodisregAFSD_new_func2_sd", (DL_FUNC) &_isodisregAFSD_new_func2_sd, 2},
    {"_isodisregAFSD_new_func_single_grid", (DL_FUNC) &_isodisregAFSD_new_func_single_grid, 4},
    {"_isodisregAFSD_new_func_single_grid_sd", (DL_FUNC) &_isodisregAFSD_new_func_single_grid_sd, 3},
    {"_isodisregAFSD_new_func_single_grid_rememeber", (DL_FUNC) &_isodisregAFSD_new_func_single_grid_rememeber, 4},
    {"_isodisregAFSD_new_func_mat", (DL_FUNC) &_isodisregAFSD_new_func_mat, 5},
    {"_isodisregAFSD_new_func_mat_sd", (DL_FUNC) &_isodisregAFSD_new_func_mat_sd, 4},
    {"_isodisregAFSD_new_func_mat_list", (DL_FUNC) &_isodisregAFSD_new_func_mat_list, 5},
    {"_isodisregAFSD_new_func_mat_list_sd", (DL_FUNC) &_isodisregAFSD_new_func_mat_list_sd, 4},
    {"_isodisregAFSD_new_func_list", (DL_FUNC) &_isodisregAFSD_new_func_list, 5},
    {"_isodisregAFSD_new_func_list_back", (DL_FUNC) &_isodisregAFSD_new_func_list_back, 5},
    {"_isodisregAFSD_new_func_list_sd_back", (DL_FUNC) &_isodisregAFSD_new_func_list_sd_back, 4},
    {"_isodisregAFSD_new_func_list_sd", (DL_FUNC) &_isodisregAFSD_new_func_list_sd, 4},
    {"_isodisregAFSD_new_func_list_mat", (DL_FUNC) &_isodisregAFSD_new_func_list_mat, 5},
    {"_isodisregAFSD_new_func_list_mat_sd", (DL_FUNC) &_isodisregAFSD_new_func_list_mat_sd, 4},
    {"_isodisregAFSD_ecdf_list_comp_class", (DL_FUNC) &_isodisregAFSD_ecdf_list_comp_class, 3},
    {"_isodisregAFSD_ecdf_list_comp_class_sd", (DL_FUNC) &_isodisregAFSD_ecdf_list_comp_class_sd, 2},
    {"_isodisregAFSD_ecdf_list_comp_class_eps", (DL_FUNC) &_isodisregAFSD_ecdf_list_comp_class_eps, 2},
    {"_isodisregAFSD_ecdf_comp", (DL_FUNC) &_isodisregAFSD_ecdf_comp, 3},
    {"_isodisregAFSD_ecdf_comp_sd", (DL_FUNC) &_isodisregAFSD_ecdf_comp_sd, 2},
    {"_isodisregAFSD_ecdf_comp_eps", (DL_FUNC) &_isodisregAFSD_ecdf_comp_eps, 2},
    {"_isodisregAFSD_ecdf_comp_class", (DL_FUNC) &_isodisregAFSD_ecdf_comp_class, 3},
    {"_isodisregAFSD_ecdf_comp_class_sd", (DL_FUNC) &_isodisregAFSD_ecdf_comp_class_sd, 2},
    {"_isodisregAFSD_ecdf_func", (DL_FUNC) &_isodisregAFSD_ecdf_func, 3},
    {"_isodisregAFSD_ecdf_func_mat", (DL_FUNC) &_isodisregAFSD_ecdf_func_mat, 3},
    {"_isodisregAFSD_ecdf_comp_class_ind", (DL_FUNC) &_isodisregAFSD_ecdf_comp_class_ind, 1},
    {"_isodisregAFSD_ecdf_list_comp_class_ind", (DL_FUNC) &_isodisregAFSD_ecdf_list_comp_class_ind, 2},
    {"_isodisregAFSD_ecdf_comp_class_eps", (DL_FUNC) &_isodisregAFSD_ecdf_comp_class_eps, 2},
    {"_isodisregAFSD_pavaCorrect", (DL_FUNC) &_isodisregAFSD_pavaCorrect, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_isodisregAFSD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
