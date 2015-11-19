#if 0
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

uint32_t prime;

struct upper_prime
{
   uint32_t a_byte:28;  /* a */
   uint32_t a_bit:4;
   uint32_t res_byte:28;
   uint32_t b_bit:4;
}__attribute__((__packed__));


struct upper_prime_bin
{
   struct upper_prime *up;
   int count;
   int size;
};

struct upper_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   uint32_t *primelist;
   uint32_t primelist_count;
   uint32_t calculated_index;

   struct upper_prime_bin *up_bins;
   int      nbins;

   uint64_t start_num;
   uint64_t end_num;
};


int
upper_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct upper_ctx *sctx = malloc(sizeof (struct upper_ctx));
   int nblocks;
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   /* TODO: Assert */
   sctx->block_size = pctx->current_block.block_size;

   sctx->primelist = malloc(sizeof(uint32_t) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->up_bins = malloc(sizeof(struct upper_prime_bin ) * 8 * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->nbins = 1;
   nblocks = (pctx->run_info.adjusted_end_num - pctx->run_info.adjusted_start_num) / pctx->current_block.block_size + 1;
   while (sctx->nbins <= (int)((pctx->run_info.max_sieve_prime / pctx->current_block.block_size) + 1) && sctx->nbins < nblocks)
      sctx->nbins <<= 1;

   sctx->up_bins = calloc (sctx->nbins, sizeof(struct upper_prime_bin));

   sctx->primelist_count = 0;
   sctx->calculated_index = 0;

   sctx->start_num = pctx->run_info.adjusted_start_num;
   sctx->end_num = pctx->run_info.adjusted_end_num;
   return 0;
}


int
upper_free(void *ctx)
{
   struct upper_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(ctx);
   return 0;
}


int
upper_skip_to(struct prime_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct upper_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   sctx->calculated_index = 0;
   return 0;
}


static void
add_upper_prime_to_bin(struct upper_prime_bin *upb, struct upper_prime *up)
{
   if (upb->count >= upb->size) {
      upb->size += 1024; /* TODO */
      upb->up = realloc(upb->up, upb->size * sizeof(struct upper_prime));
   }
   memcpy(&upb->up[upb->count], up, sizeof(struct upper_prime));
   upb->count++;
}


static void
add_to_upper_primes(struct upper_ctx *sctx, uint32_t a)
{
   int i;
   uint64_t off;
   uint64_t inA = num_to_bytes(sctx->start_num);
   uint64_t res_byte;
   uint64_t base = bytes_to_num(num_to_bytes(MAX(sctx->start_num/a, a)));
   int a_ind = pp_to_bit(a);
   const unsigned char *bits = a_x_b_bitmask[a_ind];
   struct upper_prime up;


   for (i = 0; i < 8; i++) {
      if ((off = num_to_bytes((base + ind_to_mod[i])*a)) < inA)
         continue;
      if (bytes_to_num(off) >= sctx->end_num)
         return;
      break;
   }

   off -= inA;

   /* Work out which block it is in.. */
   while (off > sctx->block_size) {
      off -= sctx->block_size;
      inA += sctx->block_size;
   }

   up.a_byte   = num_to_bytes(a);
   up.a_bit    = a_ind;
   up.b_bit    = i;
   res_byte = a_x_b_bytes[a_ind][i];
   if (res_byte > off)
      off += sctx->block_size;
   up.res_byte = (off - res_byte) % sctx->block_size;

   add_upper_prime_to_bin (&sctx->up_bins[ (inA / sctx->block_size) & (sctx->nbins - 1) ], &up);
}


int
upper_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct upper_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue
      if (primelist[*ind] > sctx->end_prime)
         return 0;
      add_to_upper_primes(sctx, primelist[*ind]);
   }
   return 0;
}


int
upper_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct upper_ctx *sctx = ctx;

   if (sctx->primelist_count == 0 || sctx->primelist[0] > pcb->sqrt_end_num)
      return 0;

   i = (pr->cb.A / BLOCK_NUMBERS) & (pr->sp.nbins - 1);
   upb = &pr->sp.up_bins[i];

   for (up = upb->up; up < upb->up + upb->count; up++) {

      *((char *)pr->cb.bitmap + up->offset) |= up->bitmask;

      offset = (int64_t)up->offset + up->a;
      j = i + offset / BLOCK_SIZE;
      up->offset = offset % BLOCK_SIZE;

      if (pr->cb.A + offset * 30 > pr->max)
         continue;

      add_upper_prime_to_bin(pr, &pr->sp.up_bins[j & (pr->sp.nbins - 1)], up);
   }
   upb->count = 0;
   /* Mark off multiples of 'a' in the block */
   uint32_t i;

   check_new_sieve_primes(sctx, pcb);

   for (i = 0; i < sctx->calculated_index; i++)
      compute_block(pcb, sctx->primelist[i], &sctx->offsets[i*8]);

   return 0;
}
#endif
