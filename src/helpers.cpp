#include <Rcpp.h>

#include <cmath>

int randWrapper(const int n) { 
  return floor(unif_rand()*n);
}

void check_interrupt_impl(void*) {
  R_CheckUserInterrupt();
}

bool check_interrupt() {
  return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}
