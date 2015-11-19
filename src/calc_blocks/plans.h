#ifndef _HARU_PRIME_PLANS_H
#define _HARU_PRIME_PLANS_H

#include <inttypes.h>

struct prime_ctx;
struct prime_thread_ctx;
struct prime_current_block;

/**
 * @FILE
 *
 * The plan essentially lists which functions to call for various sets of
 * sieving primes in order to mark off non-primes in the current block.
 *
 * The idea is to be able to focus on the functions for different sets of
 * seiving primes yet be able to test against a known working solution.
 *
 * So each section must follow the same function definition and will get called
 * by the main function in turn.
 */


/*
 * The function calc_primes will be called with the start/end numbers
 */
struct prime_plan_entry {
   const char *name;
   const uint32_t block_multiplier; /* currently unused */
   const uint32_t start;
   const uint32_t end;
   int (*init)(struct prime_ctx *pctx, uint32_t start_sieve_prime, uint32_t end_sieve_prime, void **ctx);
   int (*free)(void *ctx);
   int (*skip_to)(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx);
   int (*add_sieving_primes)(uint32_t *primes, uint32_t *ind, uint32_t max_ind, void *ctx);
   int (*calc_primes)(struct prime_thread_ctx *ptx, void *ctx);
};


#define MAX_PLAN_ENTRIES 20
struct prime_plan {
   const char *name;
   const uint32_t block_size;
   int num_entries;

   /* This is a bit dodgy but just makes it easier to specify the plan */
   struct prime_plan_entry entries[MAX_PLAN_ENTRIES];
};


/**
 * The prime plans are defined in plans.c. This is supposed to make it
 * easier to choose between plans.
 */
const struct prime_plan * get_prime_plan(const int ind);

#endif
