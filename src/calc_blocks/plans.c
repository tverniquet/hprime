#include <stdlib.h>
#include <assert.h>

#include "misc.h"
#include "plans.h"

/*
 * Different methods for calculating prime blocks on the 8/30 wheel
 */
#include "plan_register.h"

const struct prime_plan *
get_prime_plan(const int ind) {

   assert((size_t)ind < ARR_SIZEOF(plans));

   return &plans[ind];
}
