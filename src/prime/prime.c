#define _GNU_SOURCE
#include <sched.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

#include "prime.h"

#include "misc.h"
#include "plans.h"
#include "wheel.h"
#include "ctx.h"
#include "initial.h"


/******************************************************************************
 * STATIC FUNCTIONS
 *****************************************************************************/


/*
 * Applies some special manipulations to the first block
 *
 * Usually a prime is applied after prime^2 which means it won't mark itself
 * (prime*1) as not prime.
 *
 * However the lower_primes methods tended not to do this so the lower primes
 * are explicitely remarked as "prime" here.
 *
 * However, in retrospect this probably should be handled by the calc_block
 * methods since it is sort of their responsibility.
 */
static void
apply_zero_block_mod (struct prime_current_block *pcb)
{
   if (pcb->block_start_num == 0) {
      pcb->block[0] &=~ 0xFE;
      pcb->block[0] |= 1; /* 1 is not prime */
      pcb->block[1] &=~ 0xDF;
      pcb->block[2] &=~ 0xEF;
   }
}


/*
 * Always calculate and work with whole blocks. In the case that only a partial
 * block is needed (or in the special case of the first block), apply a mask to
 * the block to only include the required primes
 */
static void
apply_start_end_sets (struct prime_thread_ctx *ptx)
{
   char *t;
   uint64_t byte;
   uint32_t bit;
   /* If the block is in the middle, then return */
   if (ptx->current_block.block_start_num > ptx->main->run_info.start_num
         &&
       ptx->current_block.block_end_num < ptx->main->run_info.end_num) {
      return;
   }

   if (ptx->current_block.block_start_num < ptx->main->run_info.start_num) {
      byte = num_to_bytes(ptx->main->run_info.start_num);
      bit  = num_to_bit(ptx->main->run_info.start_num);

      t =  ptx->current_block.block + byte - num_to_bytes(ptx->current_block.block_start_num);
      *t |= (1ul << bit) - 1;

      if (t > ptx->current_block.block)
         memset (ptx->current_block.block, 0xff, (t - ptx->current_block.block));
   }

   if (ptx->current_block.block_end_num > ptx->main->run_info.end_num) {
      byte = num_to_bytes(ptx->main->run_info.end_num + 1);
      bit  = num_to_bit(ptx->main->run_info.end_num + 1);

      t =  ptx->current_block.block + byte - num_to_bytes(ptx->current_block.block_start_num);
      *t |= ~0ul << bit;

      if (t - ptx->current_block.block < ptx->current_block.block_size)
         memset (t+1, 0xff, ptx->current_block.block_size - (t - ptx->current_block.block));
   }
}

static void
calc_block_threaded (struct prime_thread_ctx *ptx)
{
   int i;
   struct prime_ctx *pm = ptx->main;
   struct prime_current_block *pcb = &ptx->current_block;

   for (i = 0; i < pm->plan_info.pp->num_entries; i++) {

      if (pcb->sqrt_end_num > pm->plan_info.pp->entries[i].start)
         pm->plan_info.pp->entries[i].calc_primes(ptx, pm->plan_info.plan_entry_ctxs[i].data);

   }
}

/*
 * Run each of the functions associated to the different prime ranges
 * to calculate each block
 */
static void
calc_block (struct prime_ctx *pctx) {
   int i;
   for (i = 0; i < pctx->plan_info.pp->num_entries; i++) {
      mark_time(&pctx->plan_info.plan_entry_ctxs[i].timer);

      if (pctx->current_block->sqrt_end_num > pctx->plan_info.pp->entries[i].start)
         pctx->plan_info.pp->entries[i].calc_primes(&pctx->threads[0], pctx->plan_info.plan_entry_ctxs[i].data);

      add_timediff(&pctx->plan_info.plan_entry_ctxs[i].timer);
   }
}


static void
get_next_block_single (struct prime_current_block *pcb) {
   pcb->block_start_num  += bytes_to_num(pcb->block_size);
   pcb->block_end_num    += bytes_to_num(pcb->block_size);
   pcb->sqrt_end_num      = sqrtl(pcb->block_end_num);
   pcb->block_start_byte += pcb->block_size;
   pcb->block_num++;
   return;
}


static void
get_next_block (struct prime_thread_ctx *ptx, struct prime_current_block *pcb)
{
   struct prime_ctx *pm = ptx->main;

   if (ptx->run_num--) {
      get_next_block_single(pcb);
      return;
   }

   ptx->run_num = pm->blocks_per_run - 1;

   pcb->block_num        = __sync_fetch_and_add(&pm->block_num, pm->blocks_per_run);
   pcb->block_start_byte = pcb->block_num * pcb->block_size;
   pcb->block_start_num  = bytes_to_num(pcb->block_start_byte);
   pcb->block_end_num    = pcb->block_start_num + bytes_to_num(pcb->block_size);
   pcb->sqrt_end_num     = sqrtl(pcb->block_end_num);
}


static void
add_sieve_primes(struct prime_ctx *ctx, uint32_t *primelist, uint32_t size)
{
   int i;
   uint32_t ind = 0;
   for (i = 0; i < ctx->plan_info.pp->num_entries && ind < size; i++)
      if (primelist[ind] >= ctx->plan_info.pp->entries[i].start)
         ctx->plan_info.pp->entries[i].add_sieving_primes(primelist, &ind, size, ctx->plan_info.plan_entry_ctxs[i].data);
}



/*
 * This is probably a case of premature optimisation...
 *
 * Anyway, the plan entries need to store the primes, so instead of getting
 * them for the whole block (~260kB worth) then passing to the plan entries
 * to store again, just get a few at a time.
 *
 * One main benefit is not having to worry how much to malloc to store the
 * primes.
 */
static void
get_some_primes_from_current_block(struct prime_thread_ctx *ptx, uint32_t *primelist, int *size, int max_size, uint32_t *cur_nbytes, uint32_t *cur_bit)
{
   uint32_t prime;
   while ((prime = get_next_prime(&ptx->current_block, cur_nbytes, cur_bit)) != 0) {

      if (prime <= ptx->main->run_info.added_sieve_primes)
         continue;

      if (prime > ptx->main->run_info.max_sieve_prime)
         return;

      primelist[(*size)++] = prime;

      if (*size >= max_size)
         return;
   }
}


static const int primes_at_a_time = 512;
static void
add_sieve_primes_from_current_block(struct prime_thread_ctx *ptx)
{
   uint32_t primelist[512];
   int size = 0;
   uint32_t byte = 0;
   uint32_t bit = 0;

   /* Don't know how much to malloc so just get a few at a time */
   do {
      get_some_primes_from_current_block(ptx, primelist, &size, primes_at_a_time, &byte, &bit);
      add_sieve_primes(ptx->main, primelist, size);
   }
   while ((size -= primes_at_a_time) == 0);

   ptx->main->run_info.added_sieve_primes = ptx->current_block.block_end_num;
}


/* Target num must be aligned to the block.. */
static void
skip_to_block(struct prime_ctx *ctx, uint64_t target_num)
{
   int i;
   for (i = 0; i < ctx->plan_info.pp->num_entries; i++) {
      ctx->plan_info.pp->entries[i].skip_to(&ctx->threads[0], target_num, ctx->plan_info.plan_entry_ctxs[i].data);
   }

   ctx->current_block->block_start_num = target_num;
   ctx->current_block->block_end_num = target_num + bytes_to_num(ctx->current_block->block_size);
   ctx->current_block->sqrt_end_num = sqrtl(ctx->current_block->block_end_num);
   ctx->block_num = target_num / 30 / ctx->current_block->block_size;
   ctx->current_block->block_num = ctx->block_num;
   ctx->threads[0].run_num = 0;
}


/******************************************************************************
 * Thread Main functions
 *****************************************************************************/

static void *
thread_calc_block(void *ctx)
{
   struct prime_thread_ctx *ptx = ctx;
   struct prime_ctx *pm = ptx->main;
   struct prime_current_block *pcb = &ptx->current_block;
   cpu_set_t cpuset;

   /*
    * Set the affinity - this appears to be needed to stop
    * hiccups. Currently no checking is done if nthreads > ncpu
    */
   CPU_ZERO(&cpuset);
   CPU_SET(ptx->thread_index, &cpuset);
   sched_setaffinity(0, sizeof cpuset, &cpuset);

   for (;;) {

      get_next_block(ptx, pcb);

      if (pcb->block_start_num > pm->run_info.end_num)
         break;

      calc_block_threaded(ptx);

      if (pcb->block_start_num == 0)
         apply_zero_block_mod (pcb);

      if (pcb->block_start_num < pm->run_info.start_num || pcb->block_end_num > pm->run_info.end_num)
         apply_start_end_sets(ptx);

      sem_post(&ptx->can_start_result);
      sem_wait(&ptx->can_start_result_next);
   }

   return NULL;
}


struct thread_and_fn {
   struct prime_thread_ctx *ptx;
   int (*fn)(struct prime_thread_ctx *, void *);
   void *th;
};


static void *
thread_calc_block_v2(void *ctx)
{
   struct thread_and_fn *f = ctx;
   struct prime_thread_ctx *ptx = f->ptx;
   struct prime_ctx *pm = ptx->main;
   struct prime_current_block *pcb = &ptx->current_block;
   cpu_set_t cpuset;

   /*
    * Set the affinity - this appears to be needed to stop
    * hiccups. Currently no checking is done if nthreads > ncpu
    */
   CPU_ZERO(&cpuset);
   CPU_SET(ptx->thread_index, &cpuset);
   sched_setaffinity(0, sizeof cpuset, &cpuset);

   for (;;) {

      get_next_block(ptx, pcb);

      if (pcb->block_start_num > pm->run_info.end_num)
         break;

      calc_block_threaded(ptx);

      if (pcb->block_end_num >= pm->run_info.end_num)
         apply_start_end_sets(ptx);

      if (f->fn(f->ptx, f->th))
         break;
   }

   ptx->thread_index = 0;
   return NULL;
}


/******************************************************************************
 * The main "get_next_block" functions
 *****************************************************************************/


/*
 * Calculates the next block in the requested prime range
 *
 * The desired usage is for the calling function to do:
 *
 *    while (get_next_block())
 *      do_something_with_block()
 *
 *
 * There is a question of how and when to calculate the sieveing primes. One
 * mechanism would be to have the calling function do:
 *
 *   calculate_sieving_primes()
 *
 *    while (get_next_block())
 *      do_something_with_block()
 *
 * For whatever reason I decided that if the sieving prime blocks and the
 * output blocks overlap they shouldn't need to be calculated twice.
 *
 * However it ended up being a bit complicated to implement. I handled the
 * tricky-ness with a state machine with the state kept by the calling
 * function. There are essentially 3 phases:
 *  0 - init
 *  1 - calculate sieving primes
 *  2 - calculate requested prime range
 *
 */
static int
calc_next_block_single (struct prime_ctx *ctx)
{
   /* Add initial primes so there is enough to calculate the first block */
   switch (ctx->run_state) {

      /*
       * Case 0 Add the initial primes
       *
       * This will be enough primes to calculate the first block, at which time any
       * further primes will also be added
       */
      case 0 : {
         add_sieve_primes(ctx, initial_primes, ARR_SIZEOF(initial_primes));
         ctx->run_info.added_sieve_primes = initial_primes[ARR_SIZEOF(initial_primes) - 1];

         /* Fallthrough */
      }

      /*
       * Case 1 Calculate sieving primes
       *
       * Process one block at a time until all sieving primes have been added.
       * If the block is in the range of the final primes, then this returns
       * once the block is calculated.
       */
      case 1 : {

         if (ctx->run_state == 1)
            get_next_block_single(&ctx->threads[0].current_block);

         ctx->run_state = 1;

         /* Always start at the 0th block in order to calculate the sieving primes */
         while (ctx->run_info.added_sieve_primes < ctx->run_info.max_sieve_prime) {

            calc_block(ctx);
            apply_zero_block_mod (&ctx->threads[0].current_block);

            add_sieve_primes_from_current_block(&ctx->threads[0]);

            if (ctx->current_block->block_start_num >= ctx->run_info.start_num) {
               apply_zero_block_mod (&ctx->threads[0].current_block);
               apply_start_end_sets(&ctx->threads[0]);
               return 1;
            }

            get_next_block_single(&ctx->threads[0].current_block);
         }

         /*
          * If the next block is to be counted, then continue. Else each plan entry
          * will need to be made aware of the jump
          */
         if (ctx->current_block->block_start_num > ctx->run_info.end_num
                ||
             ctx->current_block->block_end_num < ctx->run_info.start_num) {
            skip_to_block(ctx, ctx->run_info.adjusted_start_num);
         }
         /* Fallthrough */
      }

      /*
       * Case 2 Calculate the next block in the requested range
       *
       * Process one block at a time until all sieving primes have been added.
       * If the block is in the range of the final primes, then this returns
       * once the block is calculated.
       */
      case 2 : {

         if (ctx->run_state == 2)
            get_next_block_single(&ctx->threads[0].current_block);

         ctx->run_state = 2;

         if (ctx->current_block->block_start_num >= ctx->run_info.end_num)
            return 0;

         calc_block(ctx);
         apply_zero_block_mod (&ctx->threads[0].current_block);
         apply_start_end_sets(&ctx->threads[0]);
         return 1;
      }
   }
   /* not reached */
   return -1;
}


/* Thread id that is currently processing the next block */
static int
get_next_thread_id(struct prime_ctx *ctx)
{
   uint32_t i;
   for (;;) {
      for (i = 0; i < ctx->num_threads; i++)
         if (ctx->threads[i].current_block.block_num == ctx->process_block_num)
            return i;
      sched_yield();
   }

}


/*
 * The threaded version of "get_next_block"
 *
 * The general idea is:
 *   - wait for block to be processed
 *   - return (calling process uses the current block)
 *   - when next called, signal to thread that the block is done
 *
 */
static int
calc_next_block_th (struct prime_ctx *ctx)
{
   uint32_t i;
   /*
    * Init first time
    */
   if (ctx->run_state == 0) {

      ctx->blocks_per_run = 1;
      ctx->process_block_num = 0;

      ctx->blocks_per_run = 1;

      /*
       * Calculate all sieveing primes first. (Just makes it easier..)
       */
      ctx->block_num = 0;
      ctx->current_block = &ctx->threads[0].current_block;

      add_sieve_primes(ctx, initial_primes, ARR_SIZEOF(initial_primes));
      ctx->run_info.added_sieve_primes = initial_primes[ARR_SIZEOF(initial_primes) - 1];

      do {
         get_next_block(&ctx->threads[0], &ctx->threads[0].current_block);
         calc_block_threaded(&ctx->threads[0]);
         apply_zero_block_mod (&ctx->threads[0].current_block);
         add_sieve_primes_from_current_block(&ctx->threads[0]);
      }
      while (ctx->threads[0].current_block.block_end_num < ctx->run_info.max_sieve_prime);

      ctx->run_state = 1;

      /*
       * I found I needed to reset all of the entries because of calculating the sieving primes above
       */
      ctx->block_num = ctx->run_info.start_num / 30 / ctx->threads[0].current_block.block_size;
      ctx->process_block_num = ctx->block_num;
      for (i = 0; i < (uint32_t)ctx->plan_info.pp->num_entries; i++)
         ctx->plan_info.pp->entries[i].skip_to(&ctx->threads[0], ctx->block_num * ctx->threads[0].current_block.block_size * 30, ctx->plan_info.plan_entry_ctxs[i].data);

      /*
       * Now start the threads
       */
      for (i = 0; i < ctx->num_threads; i++) {
         ctx->threads[i].run_num = 0;
         ctx->threads[i].current_block.block_num = INT64_MAX;
         pthread_create(&ctx->threads[i].hdl, NULL, thread_calc_block, &ctx->threads[i]);
      }
   }
   else {

      /*
       * Clean up if finished
       */
      if (ctx->current_block->block_end_num >= ctx->run_info.end_num) {
         for (i = 0; i < ctx->num_threads; i++)
            sem_post(&ctx->threads[i].can_start_result_next);

         for (i = 0; i < ctx->num_threads; i++)
            pthread_join(ctx->threads[i].hdl, NULL);
         return 0;
      }
      /*
       * This says we are done processing the block and the thread
       * can move on,
       */
      sem_post(&ctx->threads[ctx->thread_i].can_start_result_next);
      ctx->process_block_num++;
   }


   ctx->thread_i = get_next_thread_id(ctx);
   sem_wait(&ctx->threads[ctx->thread_i].can_start_result);
   ctx->current_block = &ctx->threads[ctx->thread_i].current_block;
   return 1;
}


/******************************************************************************
 *
 * EXTERNAL FUNCTIONS
 *
 *****************************************************************************/

int
calc_next_block (struct prime_ctx *ctx)
{
   return ctx->num_threads == 0 ? calc_next_block_single(ctx)
                                : calc_next_block_th(ctx);
}


int
calc_blocks (struct prime_ctx *ctx, int (*fn)(struct prime_thread_ctx *, void *), void *thunk)
{
   uint32_t i;
   struct thread_and_fn *trfn;

   /*
    * A bit dodgy - just making sure it run a little bit faster for both 10^9 and 10^10 */
   ctx->blocks_per_run = ctx->run_info.end_num > 32*1024*32*1024 ? 64 : 1;

   trfn = calloc(sizeof *trfn, ctx->num_threads);

   ctx->current_block = &ctx->threads[0].current_block;

   add_sieve_primes(ctx, initial_primes, ARR_SIZEOF(initial_primes));
   ctx->run_info.added_sieve_primes = initial_primes[ARR_SIZEOF(initial_primes) - 1];

   ctx->block_num = 0;
   for (i = 0; i < ctx->num_threads; i++)
      ctx->threads[i].run_num = 0;

   for (i = 0; i < ctx->num_threads; i++) {
      trfn[i].fn = fn;
      trfn[i].th = thunk;
      trfn[i].ptx = &ctx->threads[i];
   }
   do {
      get_next_block(&ctx->threads[0], &ctx->threads[0].current_block);
      calc_block_threaded(&ctx->threads[0]);
      apply_zero_block_mod (&ctx->threads[0].current_block);
      add_sieve_primes_from_current_block(&ctx->threads[0]);
   } while (ctx->current_block->block_end_num < ctx->run_info.max_sieve_prime);

   for (i = 1; i < ctx->num_threads; i++)
      pthread_create(&ctx->threads[i].hdl, NULL, thread_calc_block_v2, &trfn[i]);

   fn(&ctx->threads[0], thunk);
   thread_calc_block_v2(&trfn[0]);

   for (i = 1; i < ctx->num_threads; i++) {
      pthread_join(ctx->threads[i].hdl, NULL);
   }

   free(trfn);

   return 0;
}


void
slow_test_primes(struct prime_ctx *pctx)
{
   struct prime_current_block *pcb = pctx->current_block;
   uint64_t start_num = MAX(pctx->run_info.start_num, pcb->block_start_num);
   uint64_t byte_1 = num_to_bytes(start_num);
   uint32_t bit_1 = num_to_bit(start_num);
   uint64_t check_prime;

   int is_prime;

   uint64_t byte_2 = 0;
   uint32_t bit_2 = 0;
   uint64_t check_num;

   for (check_prime = byte_1 * 30 + ind_to_mod[bit_1];
        byte_1 < num_to_bytes(pcb->block_end_num);
        check_prime = get_next_possible_prime(&byte_1, &bit_1)) {

      if (check_prime == 1)
         continue;

      is_prime = (*(pcb->block + byte_1 - num_to_bytes(pcb->block_start_num)) & (1 << bit_1)) == 0;

      byte_2 = 0;
      bit_2 = 0;
      if (is_prime) {
         while ((check_num = get_next_possible_prime(&byte_2, &bit_2)) <= sqrtl(check_prime))
            if ((check_prime % check_num) == 0)
               exit_error("%"PRIu64" [%"PRIu64":%"PRIu64":%"PRIu64"] divides %"PRIu64" [%"PRIu64":%"PRIu64"] by %"PRIu64" [%"PRIu64":%"PRIu64"]\n",
                     check_prime, num_to_bytes(check_prime), num_to_bytes(check_prime) - num_to_bytes(pcb->block_start_num), num_to_bit(check_prime),
                     check_num, num_to_bytes(check_num), num_to_bit(check_num),
                     check_prime / check_num, num_to_bytes(check_prime/check_num), num_to_bit(check_prime/check_num));
      }
      else {
         while ((check_num = get_next_possible_prime(&byte_2, &bit_2)) <= sqrtl(check_prime))
            if ((check_prime % check_num) == 0)
               break;

         if (check_num > sqrtl(check_prime))
            exit_error("%"PRIu64" [%"PRIu64":%"PRIu64":%"PRIu64"] is prime\n", check_prime, num_to_bytes(check_prime), num_to_bytes(check_prime) - num_to_bytes(pcb->block_start_num), num_to_bit(check_prime));
      }
   }
}
