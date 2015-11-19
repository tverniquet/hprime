#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct simple_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t *primelist;
   uint32_t primelist_count;
};


int
simple_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct simple_ctx *sctx = malloc(sizeof (struct simple_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   sctx->primelist = malloc(sizeof(uint32_t) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;

   return 0;
}


int
simple_free(void *ctx)
{
   struct simple_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(ctx);
   return 0;
}


int
simple_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct simple_ctx *sctx = ctx;
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
simple_skip_to(struct prime_thread_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   return 0;
}


static void
mark_off_prime(struct prime_current_block *pcb, uint32_t sieve_prime)
{
   uint32_t sieve_prime_bit;

   int j;
   uint64_t multiplier;
   int64_t  offset;

   sieve_prime_bit   = pp_to_bit(sieve_prime);

   multiplier = FLOOR_TO(MAX(pcb->block_start_num / sieve_prime, sieve_prime), 30);

   for (j = 0; j < 8; j++) {
      offset = (int64_t)((multiplier + ind_to_mod[j]) * sieve_prime)/30 - (int64_t)pcb->block_start_num/30;
      while (offset < 0)
         offset += sieve_prime;
      while (offset < pcb->block_size) {
         *(pcb->block + offset) |= a_x_b_bitmask[sieve_prime_bit][j];
         offset += sieve_prime;
      }
   }
}


int
simple_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct simple_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i;

   if (sctx->primelist_count == 0)
      return 0;

   if (sctx->primelist[0] == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   for (i = 0; i < sctx->primelist_count; i++) {
      if (sctx->primelist[i] > ptx->current_block.sqrt_end_num)
         break;
      mark_off_prime(&ptx->current_block, sctx->primelist[i]);
   }
   return 0;
}
