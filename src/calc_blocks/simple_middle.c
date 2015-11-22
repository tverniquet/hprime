#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> /* XXX Debug */

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct prime_and_n
{
   uint16_t prime;
   uint16_t n;
}__attribute__((__packed__));

struct simple_middle_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   struct prime_and_n *primelist;
   uint16_t *offsets;
   uint32_t primelist_count;
   uint32_t calculated_index;
};


int
simple_middle_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct simple_middle_ctx *sctx = malloc(sizeof (struct simple_middle_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime < UINT16_MAX);

   sctx->block_size = pctx->threads[0].current_block.block_size;

   sctx->primelist = malloc(sizeof(struct prime_and_n) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->offsets = malloc(sizeof(uint16_t ) * 8 * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;
   sctx->calculated_index = 0;
   return 0;
}


int
simple_middle_free(void *ctx)
{
   struct simple_middle_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(sctx->offsets);
   FREE(ctx);
   return 0;
}


int
simple_middle_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct simple_middle_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      sctx->primelist[sctx->primelist_count].prime = primelist[*ind];
      sctx->primelist[sctx->primelist_count].n = sctx->block_size / primelist[*ind];
      sctx->primelist_count++;
   }
   return 0;
}


int
simple_middle_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;

   struct simple_middle_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   sctx->calculated_index = 0;
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
   const unsigned char *rbits = a_x_b_bits[pp_to_bit(sieve_prime)];
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
      offsets[rbits[i]] = bms[i] - (pcb->block + pcb->block_size);
}


static void
check_new_sieve_primes(struct simple_middle_ctx *sctx, struct prime_current_block *pcb)
{
   int bit;
   uint32_t sieve_prime;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;

   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);

   for ( ; sctx->calculated_index < sctx->primelist_count; sctx->calculated_index++) {

      sieve_prime = sctx->primelist[sctx->calculated_index].prime;

      if (sieve_prime > pcb->sqrt_end_num)
         return;

      multiplier_base = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));

      for (bit = 0; bit < 8; bit++) {
         multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * sieve_prime);
         while (multiplier_byte < block_start_byte)
            multiplier_byte += sieve_prime;

         while (multiplier_byte - block_start_byte >= UINT16_MAX)
            multiplier_byte -= sieve_prime;

         sctx->offsets[8*sctx->calculated_index + bit] = multiplier_byte - block_start_byte;
      }
      compute_block_first_time(pcb, sctx->primelist[sctx->calculated_index].prime, &sctx->offsets[sctx->calculated_index*8]);
   }
}


static void
compute_block(struct prime_current_block *pcb, uint32_t sieve_prime, uint16_t *offsets, uint32_t n)
{
   char *bms[8];
   int i;
   /*const unsigned char *bits = a_x_b_bitmask[pp_to_bit(sieve_prime)];*/
   const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

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
simple_middle_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct simple_middle_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i;

   if (sctx->primelist_count == 0)
      return 0;

   if (sctx->primelist[0].prime == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   for (i = 0; i < sctx->calculated_index; i++)
      compute_block(&ptx->current_block, sctx->primelist[i].prime, &sctx->offsets[i*8], sctx->primelist[i].n);

   check_new_sieve_primes(sctx, &ptx->current_block);

   return 0;
}
