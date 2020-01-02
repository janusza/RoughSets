
#include <Rcpp.h>

RcppExport SEXP _RoughSets_compute_indiscernibility(SEXP inputSEXP, SEXP attr_valSEXP, SEXP unique_attr_valSEXP);
RcppExport SEXP _RoughSets_compute_chaos(SEXP inputSEXP, SEXP dec_valSEXP, SEXP unique_dec_valSEXP);

RcppExport void chooseBestCutC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
RcppExport void chooseCutCandidatesC(void *, void *, void *, void *, void *);
RcppExport void computeIndiscernibilityAndChaos(void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"_RoughSets_compute_indiscernibility", (DL_FUNC) &_RoughSets_compute_indiscernibility, 3},
    {"_RoughSets_compute_chaos", (DL_FUNC) &_RoughSets_compute_chaos, 3},
    {NULL, NULL, 0}
};

static const R_CMethodDef CEntries[] = {
    {"chooseBestCutC",                  (DL_FUNC) &chooseBestCutC,                  14},
    {"chooseCutCandidatesC",            (DL_FUNC) &chooseCutCandidatesC,             5},
    {"computeIndiscernibilityAndChaos", (DL_FUNC) &computeIndiscernibilityAndChaos,  8},
    {NULL, NULL, 0}
};

RcppExport void R_init_RoughSets(DllInfo *dll) {
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, TRUE);
}

