#ifndef _HARU_CTX_H
#define _HARU_CTX_H

#include <inttypes.h>
#include <pthread.h>
#include <semaphore.h>

#include "misc.h"
#include "wheel.h"

struct prime_plan;

/*
 * The main state for storing information about the run
 */

struct prime_run_info
{
   uint64_t start_num;
   uint64_t adjusted_start_num;
   uint64_t end_num;
   uint64_t adjusted_end_num;
   uint64_t num_blocks;

   uint32_t max_sieve_prime;
   uint32_t added_sieve_primes;
};


struct prime_plan_data
{
   struct timespot timer;
   void           *data;
};


struct prime_plan_info
{
   const struct prime_plan *pp;
   struct prime_plan_data *plan_entry_ctxs;
};


struct prime_results
{
   uint64_t count;
};


/* Threads - initial attempt:
 *
 *  Threads to be striped. This allows each prime range to calculate
 *  the next set of offsets based on a set number.
 *
 *  |__|__|__|__|__|__|__|__|__|__|
 *   t0 t1 t2 t3 t0 t1 t2 t3 t0 t1 ...
 *   ^-----------^------------^
 *
 * XXX First step, have the main function calculate results
 *--
 * Results are calculated sequentially, but can be done by the thread,
 *
 *  main -> set t0 can_process_result (block 0)
 *  main -> wait for t0 can_start_next_result
 *
 *  t0 -> wait for can_process_result()
 *  to -> count block -> add to result
 *  t0 -> send can_start_next_result()
 *
 *  *** or **
 *  t0 -> wait for can_process_result()
 *  t0 -> send can_start_next_result()
 *  to -> count block -> sync_add to result
 *--
 * XXX Second attempt
 *--
 * threads do a set number of blocks each time, they wear the
 * cost of resetting the state between blocks in a run
 * |__|__|__|__|__|__|__|__|__|
 *  t0 t0 t0 t1 t1 t1 t2 t2 t2
 * |--------|
 *  run of 3
 *
 * The next run can be chosen by a sync_add, so active threads
 * get more work
 */

struct prime_ctx;


struct prime_thread_ctx
{
   union {
      char padding[128];
      struct {
         pthread_t hdl;
         struct prime_ctx *main;
         struct prime_current_block current_block;
         sem_t can_start_result;
         int run_num;
         uint32_t thread_index;
      };
   }__attribute__((__packed__));
   union {
      char padding2[64];
      sem_t can_start_result_next;
   }__attribute__((__packed__));
}__attribute__((__packed__));


struct prime_ctx
{
   uint64_t block_num;
   uint64_t process_block_num;
   uint32_t num_threads;
   uint32_t blocks_per_run;
   int thread_i;
   struct prime_run_info      run_info;
   struct prime_plan_info     plan_info;
   struct prime_results       results;
   struct prime_current_block *current_block;
   struct prime_thread_ctx   *threads;
   int run_state;
};



void init_context(struct prime_ctx *pctx, uint64_t start, uint64_t end, int nthreads, const struct prime_plan *pp);
void free_context(struct prime_ctx *pctx);

#endif
