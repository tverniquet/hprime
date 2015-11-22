#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> /* XXX Debug */

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct prime_and_offset
{
   uint16_t prime;
   uint16_t offset;
}__attribute__((__packed__));

struct read_offs_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   struct prime_and_offset *primelist;
   uint16_t *offsets;
   uint32_t primelist_count;
   uint32_t calculated_index;
   uint32_t top10_index[10];
   uint32_t top10_count[10];
};


int
read_offs_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct read_offs_ctx *sctx = calloc(1, sizeof (struct read_offs_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime <= pctx->current_block->block_size);
   assert(sctx->end_prime <= UINT16_MAX);

   sctx->primelist = malloc(sizeof(struct prime_and_offset) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->offsets = malloc(sizeof(uint16_t ) * 8 * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;
   sctx->calculated_index = 0;
   return 0;
}


int
read_offs_free(void *ctx)
{
   struct read_offs_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(sctx->offsets);
   FREE(ctx);
   return 0;
}


int
read_offs_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct read_offs_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      sctx->primelist[sctx->primelist_count].prime = primelist[*ind];
      sctx->primelist[sctx->primelist_count].offset = primelist[*ind];
      sctx->primelist_count++;
   }
   return 0;
}


int
read_offs_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;

   struct read_offs_ctx *sctx = ctx;

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
   /* For simplicity, write the offsets the first time, with an 'a' offset of 0,
    * then apply the offset for all the other times.
    */
   for (i = 0; i < 8; i++)
      offsets[rbits[i]] = bms[i] - (pcb->block + pcb->block_size);
}


static void
check_new_sieve_primes(struct read_offs_ctx *sctx, struct prime_current_block *pcb)
{
   int bit;
   uint32_t sieve_prime;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;

   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);
   int i;

   for ( ; sctx->calculated_index < sctx->primelist_count; sctx->calculated_index++) {

      sieve_prime = sctx->primelist[sctx->calculated_index].prime;
      sctx->primelist[sctx->primelist_count].offset = sieve_prime;

      if (sieve_prime > pcb->sqrt_end_num)
         return;

      for (i = 10; i > 0; i--)
         if (sieve_prime <= pcb->block_size / i)
            sctx->top10_index[i-1] = sctx->calculated_index + 1;

      sctx->top10_count[9] = sctx->top10_index[9];
      for (i = 8; i >= 0; i--)
         if (sctx->top10_index[i] > 0)
            sctx->top10_count[i] = sctx->top10_index[i] - sctx->top10_index[i+1];

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


/* the idea is that the inline remove the switch */
static inline void __attribute__((always_inline))
compute_block_n(struct prime_current_block *pcb, struct prime_and_offset *po, uint16_t *offsets, int n)
{
   char *bms[8];
   int i;
   const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   uint32_t sieve_prime = po->prime;
   int32_t offset = po->offset - sieve_prime;
#if 0
   char *end = pcb->block + pcb->block_size;
#endif

   for (i = 0; i < 8; i++) {
      bms[i] = (char *)pcb->block + (offsets[i] + offset) + ((int32_t)(offsets[i] + offset)>>31 & sieve_prime);
      *bms[i] |= bits[i];
   }

   switch (n) {
      case 10:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 9:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 8:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 7:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 6:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 5:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 4:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
      case 3:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];

         /* This will overrun the buffer and cause more cache misses (when viewed with cachegrind)
          * HOWEVER, it seems faster than any check that I can come up with to conditionally set
          * the last bit.
          */
#if 1
      case 2:
         for (i = 0; i < 8; i++)
            *(bms[i] += sieve_prime) |= bits[i];
#endif
   }

#if 0
   for (i = 0; i < 8; i++)
      *(bms[i] + (((bms[i] - end)>>31) & sieve_prime)) |= bits[i];
#endif
   offset += n*sieve_prime - pcb->block_size;
   offset += (offset>>31) & sieve_prime;

   po->offset = offset;
}


static void
compute_block_low(struct prime_current_block *pcb, struct prime_and_offset *po, uint16_t *offsets)
{
   char *bms[8];
   int i;
   const unsigned char bits[] = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   int32_t offset;
   int32_t sieve_prime = po->prime;
   int n;

   offset = po->offset - sieve_prime;

   /*
    * po->offset is relative to the start of the block. It is positive but
    * smaller than sieve_prime, meaning offset = po->offset - sieve_prime is
    * not positive.
    *
    * If (offset + offsets[i]) is still negative then sieve_prime is added
    * which will make it positive. The idea with the masking was to avoid
    * branching (and the >>31 produces the mask for -ve numbers).
    */
   for (i = 0; i < 8; i++) {
      bms[i] = (char *)pcb->block + (offsets[i] + offset) + ((int32_t)(offsets[i] + offset)>>31 & sieve_prime);
      *bms[i] |= bits[i];
   }

   n = pcb->block_size / sieve_prime;
   offset += sieve_prime - (pcb->block_size % sieve_prime);
   offset += (offset>>31) & sieve_prime;

   while (n--)
      for (i = 0; i < 8; i++)
         *(bms[i] += sieve_prime) |= bits[i];

   po->offset = offset;
}



int
read_offs_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct read_offs_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   uint32_t i;
   uint16_t *offs;
   struct prime_and_offset *po;

   if (sctx->primelist_count == 0)
      return 0;

   if (sctx->primelist[0].prime == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   i = 0;

   offs = &sctx->offsets[0] - 8;
   po = sctx->primelist;

   /*
    * These didn't like being putting in a loop, presumably because of all of
    * the inlining, so I unrolled them manually
    */

   i = sctx->top10_count[9];
   while (i--)
      compute_block_low(&ptx->current_block, po++, (offs += 8));

   i = sctx->top10_count[8];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 10);

   i = sctx->top10_count[7];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 9);

   i = sctx->top10_count[6];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 8);

   i = sctx->top10_count[5];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 7);

   i = sctx->top10_count[4];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 6);

   i = sctx->top10_count[3];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 5);

   i = sctx->top10_count[2];
   while (i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 4);

   i = sctx->top10_count[1];
   while(i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 3);

   i = sctx->top10_count[0];
   while(i--)
      compute_block_n(&ptx->current_block, po++, (offs += 8), 2);

   check_new_sieve_primes(sctx, &ptx->current_block);

   return 0;
}
