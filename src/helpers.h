// helpers.h
#ifndef __HELPERS_H_INCLUDED__ 
#define __HELPERS_H_INCLUDED__

#ifdef _OPENMP
#include <omp.h>
#endif

int randWrapper(const int n);
void check_interrupt_impl(void*);
bool check_interrupt();

#endif 
