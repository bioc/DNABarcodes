#include <Rcpp.h>

#include <cmath>

void check_interrupt_impl(void*) {
  R_CheckUserInterrupt();
}

bool check_interrupt() {
  return (R_ToplevelExec(check_interrupt_impl, NULL) == FALSE);
}
