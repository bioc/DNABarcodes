#include <R.h>
#include <Rinternals.h>

#include <R_ext/Rdynload.h>

SEXP DNABarcodes_analyse_barcodes(SEXP metricSEXP, SEXP barcodesSEXP, SEXP coresSEXP, SEXP cost_subSEXP, SEXP cost_indelSEXP);
SEXP DNABarcodes_barcode_set_distances(SEXP metricSEXP, SEXP barcodesSEXP, SEXP coresSEXP, SEXP cost_subSEXP, SEXP cost_indelSEXP);
SEXP DNABarcodes_code_falsification(SEXP str_barcodesSEXP, SEXP dSEXP, SEXP metricSEXP, SEXP coresSEXP);
SEXP DNABarcodes_create_dnabarcodes(SEXP nSEXP, SEXP dSEXP, SEXP metricSEXP, SEXP generationSEXP, SEXP filter_tripletsSEXP, SEXP filter_gcSEXP, SEXP filter_self_complementarySEXP, SEXP str_poolSEXP, SEXP iterationsSEXP, SEXP populationSEXP, SEXP coresSEXP, SEXP use_cacheSEXP, SEXP cost_subSEXP, SEXP cost_indelSEXP);
SEXP DNABarcodes_test_distance(SEXP str_seq1SEXP, SEXP str_seq2SEXP);
SEXP DNABarcodes_create_pool(SEXP nSEXP, SEXP filter_tripletsSEXP, SEXP filter_gcSEXP, SEXP filter_self_complementarySEXP, SEXP coresSEXP);
SEXP DNABarcodes_demultiplex(SEXP barcodesSEXP, SEXP readsSEXP, SEXP metricSEXP, SEXP cost_subSEXP, SEXP cost_indelSEXP);
SEXP DNABarcodes_distance(SEXP sequence1SEXP, SEXP sequence2SEXP, SEXP metricSEXP, SEXP cost_subSEXP, SEXP cost_indelSEXP);

R_CallMethodDef callMethods[]  = {
  {"DNABarcodes_analyse_barcodes", (DL_FUNC) &DNABarcodes_analyse_barcodes, 5},
  {"DNABarcodes_barcode_set_distances", (DL_FUNC) &DNABarcodes_barcode_set_distances, 5},
  {"DNABarcodes_code_falsification", (DL_FUNC) &DNABarcodes_code_falsification, 4},
  {"DNABarcodes_create_dnabarcodes", (DL_FUNC) &DNABarcodes_create_dnabarcodes, 14},
  {"DNABarcodes_test_distance", (DL_FUNC) &DNABarcodes_test_distance, 2},
  {"DNABarcodes_create_pool", (DL_FUNC) &DNABarcodes_create_pool, 5},
  {"DNABarcodes_demultiplex", (DL_FUNC) &DNABarcodes_demultiplex, 5},
  {"DNABarcodes_distance", (DL_FUNC) &DNABarcodes_distance, 5},
  {NULL, NULL, 0}
};

void R_init_DNABarcodes(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

