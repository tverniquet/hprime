#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

/*
 * This should be enough to calculate the next b_bit
 */
struct prime_store
{
   uint32_t res_byte:20;
   uint32_t prime_byte_diff:4;
   uint32_t prime_bit:4;
   uint32_t b_bit:4;
};


struct keep_byte_ctx
{
   uint32_t start_prime_byte;
   uint32_t last_prime_byte;
   uint32_t last_calculated_prime_byte;
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   struct prime_store *primelist;
   uint32_t primelist_count;
   uint32_t calculated_index;
};


int
keep_byte_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct keep_byte_ctx *sctx = malloc(sizeof (struct keep_byte_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   /* Must have only 1 bit per block */
   assert(sctx->start_prime >= pctx->current_block.block_size);
   assert(sctx->end_prime < UINT32_MAX/2);

   sctx->block_size = pctx->current_block.block_size;

   sctx->primelist = malloc(sizeof(struct prime_store) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;
   sctx->calculated_index = 0;
   return 0;
}


int
keep_byte_free(void *ctx)
{
   struct keep_byte_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(ctx);
   return 0;
}


int
keep_byte_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct keep_byte_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      if (sctx->primelist_count == 0) {
         sctx->start_prime_byte = num_to_bytes(primelist[*ind]);
         sctx->last_prime_byte = num_to_bytes(primelist[*ind]);
         sctx->last_calculated_prime_byte = num_to_bytes(primelist[*ind]);

         sctx->primelist[sctx->primelist_count].prime_byte_diff = 0;
      }
      else {
         sctx->primelist[sctx->primelist_count].prime_byte_diff = num_to_bytes(primelist[*ind]) - sctx->last_prime_byte;
         sctx->last_prime_byte = num_to_bytes(primelist[*ind]);
      }
      sctx->primelist[sctx->primelist_count].prime_bit = num_to_bit(primelist[*ind]);
      sctx->primelist_count++;
   }
   return 0;
}


int
keep_byte_skip_to(struct prime_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct keep_byte_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   sctx->calculated_index = 0;
   sctx->last_calculated_prime_byte = sctx->start_prime_byte;
   return 0;
}


static void
check_new_sieve_primes(struct keep_byte_ctx *sctx, struct prime_current_block *pcb)
{
   struct prime_store *ps;
   uint32_t sieve_prime;
   uint64_t multiplier;
   uint64_t res_byte;
   uint32_t b_bit;
   uint32_t prime_byte;

   for ( ; sctx->calculated_index < sctx->primelist_count; sctx->calculated_index++) {

      ps = &sctx->primelist[sctx->calculated_index];

      prime_byte = sctx->last_calculated_prime_byte + ps->prime_byte_diff;
      sieve_prime = bytes_to_num(prime_byte) + ind_to_mod[ps->prime_bit];

      if (sieve_prime > pcb->sqrt_end_num)
         return;

      multiplier = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));
      b_bit = 0;
      res_byte = num_to_bytes((multiplier + 1) * sieve_prime) - a_x_b_bytes[ps->prime_bit][b_bit];
      while (res_byte + a_x_b_bytes[ps->prime_bit][b_bit] < num_to_bytes(pcb->block_start_num))
         get_next_result_set(prime_byte, ps->prime_bit, &res_byte, &b_bit);

      ps->res_byte =  res_byte + 30 - num_to_bytes(pcb->block_start_num);
      ps->b_bit = b_bit;

      sctx->last_calculated_prime_byte += ps->prime_byte_diff;
   }
}


static void
compute_block(struct prime_current_block *pcb, uint32_t prime_byte, struct prime_store *ps)
{
   uint64_t block_byte = num_to_bytes(pcb->block_start_num);
   uint64_t res_byte = block_byte + ps->res_byte - 30;
   uint32_t b_bit = ps->b_bit;

   while ((res_byte + a_x_b_bytes[ps->prime_bit][b_bit]) < block_byte + pcb->block_size) {

      *(pcb->block + (res_byte + a_x_b_bytes[ps->prime_bit][b_bit] - block_byte)) |= a_x_b_bitmask[ps->prime_bit][b_bit];

      get_next_result_set(prime_byte, ps->prime_bit, &res_byte, &b_bit);
   }
   ps->res_byte = ((res_byte + 30) - block_byte) - pcb->block_size;
   ps->b_bit = b_bit;
}


int
keep_byte_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct keep_byte_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i;
   uint32_t prime_byte = sctx->start_prime_byte;

   check_new_sieve_primes(sctx, pcb);

   for (i = 0; i < sctx->calculated_index; i++) {
      prime_byte += sctx->primelist[i].prime_byte_diff;
      compute_block(pcb, prime_byte, &sctx->primelist[i]);
   }

   return 0;
}
