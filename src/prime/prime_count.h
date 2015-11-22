#ifndef _HARU_PRIME_H
#define _HARU_PRIME_H

#include <inttypes.h>

int getprimecount (int plan_index, uint64_t start, uint64_t end, uint64_t *count, int nthreads, int inorder);
int getprimecount_cmp_plan (uint64_t start, uint64_t end, uint64_t *count, int ind1, int ind2, int nthreads);

#endif

