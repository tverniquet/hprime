#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>

#include "ctx.h"
#include "plans.h"
#include "prime.h"


/*
 * Count the number of primes in the block (the number of bits not set)
 */
static uint64_t
count_block (struct prime_current_block *block) {
   uint64_t    *p = (uint64_t *)block->block;
   int          c = block->block_size >> 3;
   uint64_t count;

   count = c * 64;
   while (c--)
      count -= __builtin_popcountl(*p++);

   return count;
}


/*
 * Handle silly "primes between 2 and 5" type situations
 */
static void
adjust_for_early_counts (struct prime_ctx *ctx)
{
   switch (ctx->run_info.start_num) {
      case 0:
      case 1:
      case 2:
         ctx->results.count++;
      case 3:
         ctx->results.count++;
      case 4:
      case 5:
         ctx->results.count++;
   };
   switch (ctx->run_info.end_num) {
      case 0:
      case 1:
         ctx->results.count--;
      case 2:
         ctx->results.count--;
      case 3:
      case 4:
         ctx->results.count--;
   };

}

struct counts {
   union {
      char padding[64];
      struct {
         uint64_t count;
         uint64_t spin;
         sem_t s1;
      };
   };
}__attribute__((__packed__));


static int
count_thr (struct prime_thread_ctx *ptx, void *th)
{
   struct counts *counts = th;
   counts[ptx->thread_index].count += count_block(&ptx->current_block);
   return 0;
}


static void
print_times (struct prime_ctx *ctx)
{
   time_t tot_c_time = 0,
          tot_r_time = 0;
   int i;
   if (ctx->num_threads != 0)
      return;
   /*
    * Printing debug here as I know about the context
    */
   for (i = 0; i < ctx->plan_info.pp->num_entries; i++) {
      tot_c_time += ctx->plan_info.plan_entry_ctxs[i].timer.c_diff;
      tot_r_time += ctx->plan_info.plan_entry_ctxs[i].timer.r_diff;
   }

   fprintf(stderr, "%s\n", ctx->plan_info.pp->name);
   fprintf(stderr, "%15s %10s %3s %10s %3s\n", "   Name       ", "   rtime ", "  %", "   ctime ", "  %");

   for (i = 0; i < ctx->plan_info.pp->num_entries; i++) {
      fprintf(stderr, "%-15s %10lu %2lu%% %10lu %2lu%%\n",
            ctx->plan_info.pp->entries[i].name,
            ctx->plan_info.plan_entry_ctxs[i].timer.r_diff,
            ctx->plan_info.plan_entry_ctxs[i].timer.r_diff * 100 / (tot_r_time?:1),
            ctx->plan_info.plan_entry_ctxs[i].timer.c_diff,
            ctx->plan_info.plan_entry_ctxs[i].timer.c_diff * 100 / (tot_c_time?:1));
   }
}


int
getprimecount (int plan_index, uint64_t start, uint64_t end, uint64_t *count, int nthreads, int inorder) {
   struct prime_ctx ctx;
   int i;

   struct counts *counts;

   init_context(&ctx, start, end, nthreads, get_prime_plan(plan_index));

   if (nthreads == 0 || inorder) {
      while (calc_next_block(&ctx))
         ctx.results.count += count_block(ctx.current_block);

      print_times(&ctx);
   }
   else {
      counts = calloc(sizeof *counts, nthreads);

      calc_blocks(&ctx, count_thr, counts);
      for (i = 0; i < nthreads; i++)
         ctx.results.count += counts[i].count;

      free(counts);
   }

   adjust_for_early_counts(&ctx);
   *count = ctx.results.count;

   free_context(&ctx);
   return 0;
}


/*
 * Debug mode - compare the times and results of the two prime plans
 */
int
getprimecount_cmp_plan (uint64_t start, uint64_t end, uint64_t *count, int ind1, int ind2, int nthreads)
{
   const struct prime_plan *pp_1 = get_prime_plan(ind1);
   const struct prime_plan *pp_2 = get_prime_plan(ind2);
   struct timespot s1, s2;
   int i;
   struct prime_ctx ctx_1, ctx_2;

   init_context(&ctx_1, start, end, nthreads, pp_1);
   init_context(&ctx_2, start, end, nthreads, pp_2);

   /* Can't compare unless block_size is the same */
   assert(ctx_1.run_info.num_blocks == ctx_2.run_info.num_blocks);
   init_time(&s2);

   while (calc_next_block(&ctx_1)) {
      init_time(&s1);
      mark_time(&s1);
      mark_time(&s2);

      fprintf(stderr, "Calculating block %"PRIu64" %3.1f%%: ", ctx_1.current_block->block_start_num, (ctx_1.current_block->block_start_num - ctx_1.run_info.adjusted_start_num) * 100.0 / ((ctx_1.run_info.adjusted_end_num - ctx_1.run_info.adjusted_start_num)+ 1));
      fflush(stderr);

      if (calc_next_block(&ctx_2) == 0)
         exit_error("plan 2 finished earlier");

      ctx_2.results.count += count_block(ctx_2.current_block);

      add_timediff(&s1);
      add_timediff(&s2);
      fprintf(stderr, " "TIME_DIFF_FMT_US" : "TIME_DIFF_FMT_MS"\n", TIME_DIFF_VALUES_US(&s1), TIME_DIFF_VALUES_MS(&s2));
#ifdef SLOW_TEST_PRIMES
      slow_test_primes(&ctx_2);
#endif

      if (memcmp(ctx_1.current_block->block, ctx_2.current_block->block, ctx_1.current_block->block_size) != 0) {
         for (i = 0; i < (int)ctx_1.current_block->block_size; i++) {
            if (ctx_1.current_block->block[i] != ctx_2.current_block->block[i]) {
               fprintf(stderr, "Differs first at byte %ld(%d) %02X vs %02X number from %lu\n", ctx_1.current_block->block_start_num/30 + i, i, ctx_1.current_block->block[i] & 0xFF, ctx_2.current_block->block[i] & 0xFF, ctx_1.current_block->block_start_num + i*30);
               break;
            }
         }

         slow_test_primes(&ctx_2);
         exit_error("Different");
      }
   }

   adjust_for_early_counts(&ctx_2);

   *count = ctx_2.results.count;

   free_context(&ctx_1);
   free_context(&ctx_2);
   return 0;
}
