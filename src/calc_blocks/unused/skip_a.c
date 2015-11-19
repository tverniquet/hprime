#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct prime_and_n
{
   uint16_t prime;
   uint16_t n;
}__attribute__((__packed__));

struct skip_a_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   struct prime_and_n *primelist[8];
   uint16_t *offsets[8];
   uint32_t primelist_count[8];
   uint32_t calculated_index[8];
};


int
skip_a_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct skip_a_ctx *sctx = malloc(sizeof (struct skip_a_ctx));
   int i;
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime < UINT16_MAX);

   sctx->block_size = pctx->current_block.block_size;

   for (i = 0; i < 8; i++) {
      sctx->primelist[i] = malloc(sizeof(struct prime_and_n) * (sctx->end_prime - start_prime) / 4 / 8 + 1000);
      sctx->offsets[i] = malloc(sizeof(uint16_t ) * 8 * (sctx->end_prime - start_prime) / 4 / 8 + 1000);
      sctx->primelist_count[i] = 0;
      sctx->calculated_index[i] = 0;
   }
   return 0;
}


int
skip_a_free(void *ctx)
{
   struct skip_a_ctx *sctx = ctx;
   int i;
   for (i = 0; i < 8; i++) {
      FREE(sctx->primelist[i]);
      FREE(sctx->offsets[i]);
   }
   FREE(ctx);
   return 0;
}


int
skip_a_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct skip_a_ctx *sctx = ctx;
   int bit;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      bit = pp_to_bit(primelist[*ind]);
      sctx->primelist[bit][sctx->primelist_count[bit]].prime = primelist[*ind];
      sctx->primelist[bit][sctx->primelist_count[bit]].n = sctx->block_size / primelist[*ind];
      sctx->primelist_count[bit]++;
   }
   return 0;
}


int
skip_a_skip_to(struct prime_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct skip_a_ctx *sctx = ctx;
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
compute_block_first_time(struct prime_current_block *pcb, uint32_t sieve_prime, uint16_t *offsets)
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
check_new_sieve_primes(struct skip_a_ctx *sctx, struct prime_current_block *pcb)
{
   int a_bit;
   int bit;
   uint32_t sieve_prime;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;

   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);

   for (a_bit = 0; a_bit < 8; a_bit++) {
      for ( ; sctx->calculated_index[a_bit] < sctx->primelist_count[a_bit]; sctx->calculated_index[a_bit]++) {

         sieve_prime = sctx->primelist[a_bit][sctx->calculated_index[a_bit]].prime;

         if (sieve_prime > pcb->sqrt_end_num)
            break;

         multiplier_base = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));

         for (bit = 0; bit < 8; bit++) {
            multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * sieve_prime);
            while (multiplier_byte < block_start_byte)
               multiplier_byte += sieve_prime;

            while (multiplier_byte - block_start_byte >= UINT16_MAX)
               multiplier_byte -= sieve_prime;

            sctx->offsets[a_bit][8*sctx->calculated_index[a_bit] + bit] = multiplier_byte - block_start_byte;
         }
         compute_block_first_time(pcb, sctx->primelist[a_bit][sctx->calculated_index[a_bit]].prime, &sctx->offsets[a_bit][sctx->calculated_index[a_bit]*8]);
      }
   }
}


static void
compute_block(struct prime_current_block *pcb, uint32_t sieve_prime, uint16_t *offsets, uint32_t n, int a_bit)
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
skip_a_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct skip_a_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i, a_bit;

   if (sctx->primelist_count[1] > 0 && sctx->primelist[1][0].prime == 7)
      memset(pcb->block, 0, pcb->block_size);


   for (a_bit = 0; a_bit < 8; a_bit++) {
      if (sctx->primelist_count[a_bit] == 0)
         continue;

      for (i = 0; i < sctx->calculated_index[a_bit]; i++)
         compute_block(pcb, sctx->primelist[a_bit][i].prime, &sctx->offsets[a_bit][i*8], sctx->primelist[a_bit][i].n, a_bit);
   }
   check_new_sieve_primes(sctx, pcb);

   return 0;
}
