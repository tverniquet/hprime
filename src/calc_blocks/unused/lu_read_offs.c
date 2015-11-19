#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"


struct lu_read_offs_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   uint32_t *primelist;
   uint32_t *offsets;
   int32_t *offset;
   uint32_t primelist_count;
   uint32_t calculated_index;
};


int
lu_read_offs_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct lu_read_offs_ctx *sctx = malloc(sizeof (struct lu_read_offs_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   /* Must have only 1 bit per block */
   assert(sctx->start_prime >= pctx->current_block.block_size);
   assert(sctx->end_prime < UINT32_MAX/2);

   sctx->block_size = pctx->current_block.block_size;

   sctx->primelist = malloc(sizeof(uint32_t) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->offset  = malloc(sizeof(uint32_t) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->offsets = malloc(sizeof(uint32_t ) * 8 * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;
   sctx->calculated_index = 0;
   return 0;
}


int
lu_read_offs_free(void *ctx)
{
   struct lu_read_offs_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(sctx->offsets);
   FREE(ctx);
   return 0;
}


int
lu_read_offs_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct lu_read_offs_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      sctx->primelist[sctx->primelist_count] = primelist[*ind];
      sctx->primelist_count++;
   }
   return 0;
}


int
lu_read_offs_skip_to(struct prime_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct lu_read_offs_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   sctx->calculated_index = 0;
   return 0;
}


static void
check_new_sieve_primes(struct lu_read_offs_ctx *sctx, struct prime_current_block *pcb)
{
   int bit;
   uint32_t sieve_prime;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;

   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);
   /*const unsigned char *rbit;*/

   for ( ; sctx->calculated_index < sctx->primelist_count; sctx->calculated_index++) {

      sieve_prime = sctx->primelist[sctx->calculated_index];
      sctx->offset[sctx->calculated_index] = 0;
      /*rbit = a_x_b_bits[pp_to_bit(sieve_prime)];*/

      if (sieve_prime > pcb->sqrt_end_num)
         return;

      multiplier_base = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));

      for (bit = 0; bit < 8; bit++) {
         multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * sieve_prime);

         while (multiplier_byte < block_start_byte)
            multiplier_byte += sieve_prime;

         while (multiplier_byte - block_start_byte >= sieve_prime)
            multiplier_byte -= sieve_prime;

         /*sctx->offsets[8*sctx->calculated_index + rbit[bit]] = multiplier_byte - block_start_byte;*/
         sctx->offsets[8*sctx->calculated_index + bit] = multiplier_byte - block_start_byte;
      }
   }
}


static void
compute_block(struct prime_current_block *pcb, uint32_t sieve_prime, uint32_t *offsets, int32_t *offset)
{
   int i;
   /*const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};*/
   const unsigned char *bits = a_x_b_bitmask[pp_to_bit(sieve_prime)];
   int32_t off;

   for (i = 0; i < 8; i++) {
      off = (offsets[i] - *offset) +  ((((int32_t)offsets[i] - *offset)>>31)&sieve_prime);
      if (off < (int32_t)pcb->block_size)
         *(pcb->block + off) |= bits[i];
   }

   *offset += pcb->block_size;
   *offset -= (((int32_t)sieve_prime - (int32_t)*offset) >> 31) & sieve_prime;
}


int
lu_read_offs_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct lu_read_offs_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i;

   check_new_sieve_primes(sctx, pcb);

   for (i = 0; i < sctx->calculated_index; i++)
      compute_block(pcb, sctx->primelist[i], &sctx->offsets[i<<3], &sctx->offset[i]);

   return 0;
}
