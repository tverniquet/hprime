#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct lu_skip_a_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   uint16_t *primeskip[8];
   uint32_t *offsets[8];
   uint32_t primelist_first[8];
   uint32_t primelist_latest[8];
   uint32_t primelist_calculated[8];
   uint32_t primelist_count[8];
   uint32_t calculated_index[8];
};


int
lu_skip_a_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct lu_skip_a_ctx *sctx = malloc(sizeof (struct lu_skip_a_ctx));
   int i;
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime < UINT32_MAX);

   sctx->block_size = pctx->current_block.block_size;

   for (i = 0; i < 8; i++) {
      sctx->primeskip[i] = malloc(sizeof(uint16_t) * (sctx->end_prime - start_prime) / 4 / 8 + 1000);
      sctx->offsets[i] = malloc(sizeof(uint16_t ) * 8 * (sctx->end_prime - start_prime) / 4 / 8 + 1000);
      sctx->primelist_count[i] = 0;
      sctx->calculated_index[i] = 0;
   }
   return 0;
}


int
lu_skip_a_free(void *ctx)
{
   struct lu_skip_a_ctx *sctx = ctx;
   int i;
   for (i = 0; i < 8; i++) {
      FREE(sctx->primeskip[i]);
      FREE(sctx->offsets[i]);
   }
   FREE(ctx);
   return 0;
}


int
lu_skip_a_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct lu_skip_a_ctx *sctx = ctx;
   int bit;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      bit = pp_to_bit(primelist[*ind]);

      if (sctx->primelist_count[bit] == 0) {
         sctx->primelist_first[bit] = primelist[*ind];
         sctx->primeskip[bit][sctx->primelist_count[bit]] = 0;
      }
      else {
         sctx->primeskip[bit][sctx->primelist_count[bit]] = num_to_bytes(primelist[*ind] - sctx->primelist_latest[bit]);
      }

      sctx->primelist_latest[bit] = primelist[*ind];
      sctx->primelist_count[bit]++;
   }
   return 0;
}


int
lu_skip_a_skip_to(struct prime_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct lu_skip_a_ctx *sctx = ctx;
   int i;
   /* Recalculate all of the offsets depending on the new block */
   for (i = 0; i < 8; i++)
      sctx->calculated_index[i] = 0;
   return 0;
}


/*
 * The first time we don't know how many times to mark off the prime in the
 * block (ie don't really know 'n'). Anyway, keeping this separate to the
 * calc_block function allows the "calc_block" function to make more
 * assumptions and optimise better.
 */
static void
compute_block_first_time(struct prime_current_block *pcb, uint32_t sieve_prime, uint32_t *offsets)
{
   char *bms[8];
   int i;
   const unsigned char *bits = a_x_b_bitmask[pp_to_bit(sieve_prime)];
   int done = 0;

   for (i = 0; i < 8; i++)
      bms[i] = pcb->block + offsets[i];

   while ( ! done) {
      for (i = 0; i < 8; i++) {
         if (bms[i]  >= pcb->block + pcb->block_size)
            done = 1;
         else {
            *bms[i] |= bits[i];
            bms[i]  += sieve_prime;
         }
      }
   }
   for (i = 0; i < 8; i++)
      offsets[i] = bms[i] - (pcb->block + pcb->block_size);
}


static void
check_new_sieve_primes(struct lu_skip_a_ctx *sctx, struct prime_current_block *pcb)
{
   int a_bit;
   int bit;
   uint32_t sieve_prime;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;

   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);

   for (a_bit = 0; a_bit < 8; a_bit++) {
      for ( ; sctx->calculated_index[a_bit] < sctx->primelist_count[a_bit]; sctx->calculated_index[a_bit]++) {

         if (sctx->calculated_index[a_bit] == 0) {
            sieve_prime = sctx->primelist_calculated[a_bit] = sctx->primelist_first[a_bit];
         }
         else {
            sieve_prime = sctx->primelist_calculated[a_bit] + bytes_to_num(sctx->primeskip[a_bit][sctx->calculated_index[a_bit]]);
         }

         if (sieve_prime > pcb->sqrt_end_num)
            break;

         multiplier_base = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));

         for (bit = 0; bit < 8; bit++) {
            multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * sieve_prime);
            while (multiplier_byte < block_start_byte)
               multiplier_byte += sieve_prime;

            while (multiplier_byte - block_start_byte >= UINT32_MAX)
               multiplier_byte -= sieve_prime;

            sctx->offsets[a_bit][8*sctx->calculated_index[a_bit] + bit] = multiplier_byte - block_start_byte;
         }
         compute_block_first_time(pcb, sieve_prime, &sctx->offsets[a_bit][sctx->calculated_index[a_bit]*8]);
         sctx->primelist_calculated[a_bit] = sieve_prime;
      }
   }
}


static void
compute_block(struct prime_current_block *pcb, uint32_t sieve_prime, uint32_t *offsets, uint32_t n, int a_bit)
{
   char *bms[8];
   int i;
   const unsigned char *bits = a_x_b_bitmask[a_bit];

   for (i = 0; i < 8; i++)
      bms[i] = (char *)pcb->block + offsets[i];

   while (n--) {
      for (i = 0; i < 8; i++)
         *bms[i]   |= bits[i];
      for (i = 0; i < 8; i++)
         bms[i]   += sieve_prime;
   }

   for (i = 0; i < 8; i++) {
      if (bms[i]  < (char *)pcb->block + pcb->block_size) {
         *bms[i] |= bits[i];
         bms[i]  += sieve_prime;
      }
   }

   for (i = 0; i < 8; i++)
      offsets[i] = bms[i] - ((char *)pcb->block + pcb->block_size);
}


int
lu_skip_a_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct lu_skip_a_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i, a_bit;
   uint32_t sieve_prime;

   for (a_bit = 0; a_bit < 8; a_bit++) {
      if (sctx->primelist_count[a_bit] == 0)
         continue;

      sieve_prime = sctx->primelist_first[a_bit];

      for (i = 0; i < sctx->calculated_index[a_bit]; i++) {
         sieve_prime += bytes_to_num(sctx->primeskip[a_bit][i]);
         compute_block(pcb, sieve_prime, &sctx->offsets[a_bit][i*8], 0, a_bit);
      }
   }
   check_new_sieve_primes(sctx, pcb);

   return 0;
}
