#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ctx.h"
#include "plans.h"


static void
set_run_info(struct prime_run_info *rinfo, uint64_t start, uint64_t end, uint64_t block_size)
{
   rinfo->start_num = start;
   rinfo->end_num   = end;
   rinfo->adjusted_start_num = FLOOR_TO(start, block_size);
   rinfo->adjusted_end_num   = CEIL_TO(end + 1, block_size);

   rinfo->num_blocks = (rinfo->adjusted_end_num - rinfo->adjusted_start_num) / block_size;

   rinfo->max_sieve_prime = sqrtl(end);
   rinfo->added_sieve_primes = 0;
}


static void
set_initial_block(struct prime_current_block *cb, uint32_t block_size)
{
   cb->block = (char*)aligned_alloc(32*1024, CEIL_TO(block_size + 64*1024 , 1024)) + 32*1024;
   cb->block_size = block_size;
   cb->block_start_num = 0;
   cb->block_end_num = bytes_to_num(block_size);
   cb->sqrt_end_num = sqrtl(cb->block_end_num);
}


static void
free_initial_block(struct prime_current_block *cb)
{
   if (cb->block != NULL)
      free(cb->block - 32*1024);
   bzero(cb, sizeof(*cb));
}


static void
set_plan(struct prime_ctx *pctx, const struct prime_plan *pp)
{
   int i;
   pctx->plan_info.pp = pp;
   pctx->plan_info.plan_entry_ctxs = calloc(ARR_SIZEOF(pp->entries), sizeof (struct prime_plan_data));

   for (i = 0; i < pp->num_entries; i++)
      pp->entries[i].init(pctx, pp->entries[i].start, pp->entries[i].end, &pctx->plan_info.plan_entry_ctxs[i].data);
}


static void
free_plan(struct prime_ctx *pctx, const struct prime_plan *pp)
{
   int i;
   for (i = 0; i < pp->num_entries; i++) {
      pp->entries[i].free(pctx->plan_info.plan_entry_ctxs[i].data);
      pctx->plan_info.plan_entry_ctxs[i].data = NULL;
   }

   FREE(pctx->plan_info.plan_entry_ctxs);
   bzero(&pctx->plan_info, sizeof pctx->plan_info);

}


void
init_context(struct prime_ctx *pctx, uint64_t start, uint64_t end, int nthreads, const struct prime_plan *pp)
{
   int i;
   bzero(pctx, sizeof *pctx);

   set_run_info(&pctx->run_info, start, end, bytes_to_num(pp->block_size));
   if (nthreads == 0) {
      pctx->threads = calloc(sizeof *pctx->threads, 1);
      set_initial_block (&pctx->threads[0].current_block, pp->block_size);
      pctx->current_block = &pctx->threads[0].current_block;
      pctx->threads[0].thread_index = 0;
      pctx->threads[0].main = pctx;
   }
   else {
      pctx->num_threads = nthreads;
      pctx->threads = calloc(sizeof *pctx->threads, nthreads);
      for (i = 0; i < nthreads; i++) {
         set_initial_block (&pctx->threads[i].current_block, pp->block_size);
         pctx->threads[i].thread_index = i;
         pctx->threads[i].main = pctx;
         sem_init(&pctx->threads[i].can_start_result, 0, 0);
         sem_init(&pctx->threads[i].can_start_result_next, 0, 0);
      }
      pctx->current_block = &pctx->threads[0].current_block;
   }
   set_plan(pctx, pp);
}


void
free_context(struct prime_ctx *pctx) {
   uint32_t i;
   for (i = 0; i < pctx->num_threads; i++) {
      free_initial_block(&pctx->threads[i].current_block);
      sem_destroy(&pctx->threads[i].can_start_result_next);
      sem_destroy(&pctx->threads[i].can_start_result);
   }
   if (pctx->num_threads == 0)
      free_initial_block(&pctx->threads[0].current_block);
   free (pctx->threads);
   free_plan(pctx, pctx->plan_info.pp);
   bzero(pctx, sizeof *pctx);
}
