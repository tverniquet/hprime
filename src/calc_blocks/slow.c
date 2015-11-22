#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct slow_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t *primelist;
   uint32_t primelist_count;
};


int
slow_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct slow_ctx *sctx = malloc(sizeof (struct slow_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   sctx->primelist = malloc(sizeof(uint32_t) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;
   return 0;
}


int
slow_free(void *ctx)
{
   struct slow_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(ctx);
   return 0;
}


int
slow_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct slow_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      sctx->primelist[sctx->primelist_count++] = primelist[*ind];
   }
   return 0;
}


int
slow_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;
   (void)ctx;

   return 0;
}


static void
set_bit(struct prime_current_block *pcb, uint64_t num)
{
   pcb->block[(num - pcb->block_start_num)/30] |=  1 << num_to_bit(num);
}


static void
mark_off_prime(struct prime_current_block *pcb, uint32_t sieve_prime)
{
   uint64_t b;
   int bit;

   b = MAX(pcb->block_start_num / sieve_prime, sieve_prime) / 30 * 30;

   bit = 0;
   while (sieve_prime * (b + ind_to_mod[bit]) < pcb->block_start_num)
      if ((bit = (bit + 1) % 8) == 0)
         b += 30;

   for (;;) {
      for (; bit < 8; bit++) {
         if (sieve_prime * (b + ind_to_mod[bit]) >= pcb->block_end_num)
            return;

         set_bit(pcb, sieve_prime * (b + ind_to_mod[bit]));
      }
      b += 30;
      bit = 0;
   }
}


int
slow_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct slow_ctx *sctx = ctx;
   uint32_t i;

   if (sctx->primelist[0] == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   for (i = 0; i < sctx->primelist_count; i++) {
      if (sctx->primelist[i] > ptx->current_block.sqrt_end_num)
         break;

      mark_off_prime(&ptx->current_block, sctx->primelist[i]);
   }
   return 0;
}
