#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h> /* XXX Debug */
#include <immintrin.h>
#include <avx2intrin.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"


/*
 * Was trying to make this 4 bytes. It can fit but it was taking too long to extract
 */
struct prime_and_offset
{
   uint8_t ind;
   uint8_t a_byte_diff;
   uint32_t offset;
}__attribute__((__packed__));


struct prime_list
{
   struct prime_and_offset *primes;
   uint32_t start_a_byte;
   uint32_t cur_a_byte;
   uint32_t ind_a_byte;
   uint32_t count;
   uint32_t index;
   int      last_blockno;
};


struct lu_calc_offs_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   int nthreads;
   struct prime_list **plist;
};


int
lu_calc_offs_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct lu_calc_offs_ctx *sctx = calloc(sizeof (struct lu_calc_offs_ctx), 1);
   int i;
   int k;
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime < 900000);

   sctx->block_size = pctx->current_block->block_size;
   sctx->nthreads = pctx->num_threads?:1;

   sctx->plist = calloc((sizeof *sctx->plist), sctx->nthreads);
   for (k = 0; k < sctx->nthreads; k++) {
      sctx->plist[k] = calloc((sizeof *sctx->plist[k]), 8);
      for (i = 0; i < 8; i++) {
         sctx->plist[k][i].primes = malloc(sizeof(struct prime_and_offset) * (sctx->end_prime - start_prime) / 4 / 10 + 1000);
         sctx->plist[k][i].start_a_byte = sctx->start_prime/30;
         sctx->plist[k][i].cur_a_byte = sctx->start_prime/30;
         sctx->plist[k][i].ind_a_byte = sctx->start_prime/30;
      }
   }
   return 0;
}


int
lu_calc_offs_free(void *ctx)
{
   struct lu_calc_offs_ctx *sctx = ctx;
   int i;
   int k;
   for (k = 0; k < sctx->nthreads; k++)
      for (i = 0; i < 8; i++)
         FREE(sctx->plist[k][i].primes);
   /* TODO: more freeing */
   FREE(ctx);
   return 0;
}


int
lu_calc_offs_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct lu_calc_offs_ctx *sctx = ctx;
   struct prime_list *pl;

   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      int k;

      for (k = 0; k < sctx->nthreads; k++) {
         pl = &sctx->plist[k][pp_to_bit(primelist[*ind])];
         pl->primes[pl->count].a_byte_diff = primelist[*ind]/30 - pl->cur_a_byte;
         pl->primes[pl->count].offset = 0;
         pl->cur_a_byte = primelist[*ind]/30;
         pl->count++;
      }
   }
   return 0;
}


int
lu_calc_offs_skip_to(struct prime_thread_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct lu_calc_offs_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   int i;
   for (i = 0; i < 8; i++) {
      sctx->plist[pctx->thread_index][i].index = 0;
      sctx->plist[pctx->thread_index][i].ind_a_byte = sctx->start_prime / 30;
      sctx->plist[pctx->thread_index][i].last_blockno = INT32_MAX;
   }
   return 0;
}


static void
check_new_sieve_primes_v2(struct lu_calc_offs_ctx *sctx, struct prime_current_block *pcb, int thread_id)
{
   uint32_t sieve_prime;

   uint64_t block_start_byte = pcb->block_start_byte;
   int i;
   int32_t offset;
   struct prime_list *pl;
   struct prime_and_offset *po;
   const unsigned char *bytes;
   uint32_t a_byte;
   for (i = 0; i < 8; i++) {
      pl = &sctx->plist[thread_id][i];
      a_byte = pl->ind_a_byte;
      for ( ; pl->index < pl->count; pl->index++) {
         po = &pl->primes[pl->index];
         a_byte += po->a_byte_diff;
         sieve_prime = (int32_t)a_byte * 30 + ind_to_mod[i];

         if (sieve_prime > pcb->sqrt_end_num)
            break;

         bytes = a_x_b_byte_diffs[i];

         offset = 0 - (block_start_byte % sieve_prime);

         if ((offset += a_byte) >= 0)
            po->ind = 0;
         else if ((offset += 6*a_byte + bytes[1]) >= 0)
            po->ind = 1;
         else if ((offset += 4*a_byte + bytes[2]) >= 0)
            po->ind = 2;
         else if ((offset += 2*a_byte + bytes[3]) >= 0)
            po->ind = 3;
         else if ((offset += 4*a_byte + bytes[4]) >= 0)
            po->ind = 4;
         else if ((offset += 2*a_byte + bytes[5]) >= 0)
            po->ind = 5;
         else if ((offset += 4*a_byte + bytes[6]) >= 0)
            po->ind = 6;
         else if ((offset += 6*a_byte + bytes[7]) >= 0)
            po->ind = 7;
         else if ((offset += 2*a_byte + 1) > 0)
            po->ind = 0;

         po->offset = offset;
         pl->ind_a_byte = a_byte;
      }
   }
}


static inline int __attribute__((always_inline))
check_set (struct prime_current_block *pcb, struct prime_and_offset *po, int32_t *off, int32_t a_byte_x_2, const unsigned char *bits, const unsigned char *bytes, const int ind, const int a_x)
{
      if (*off >= (int32_t)pcb->block_size) {
         po->ind = ind;
         po->offset = *off - pcb->block_size;
         return 1;
      }

      *(pcb->block + *off) |= bits[ind];
      /* I was trying to avoid mutliplication and using more registers...?? */
      switch (a_x) {
         case 6:
            *off += a_byte_x_2;
         case 4 :
            *off += a_byte_x_2;
         case 2 :
            *off += a_byte_x_2;
      }
      *off += ind == 7 ? 1 : bytes[ind+1];
      return 0;
}


/*
 * The main reason for the inline is to turn the bits and bytes into immediates. They differ
 * for the different a_bits
 *
 * Unlike primes < 32*1024, these can't be done out of order to keep the same bit pattern
 */
static inline void __attribute__((always_inline))
compute_block(struct prime_current_block *pcb, struct prime_and_offset *po, int a_byte, int a_bit)
{
   const unsigned char *bits = a_x_b_bitmask[a_bit];
   const unsigned char *bytes = a_x_b_byte_diffs[a_bit];

   int32_t a_byte_x_2 = a_byte << 1;
   int32_t off;

   off = po->offset;

   switch (po->ind) {
      case 0 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 0, 6) ) return;
      case 1 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 1, 4) ) return;
      case 2 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 2, 2) ) return;
      case 3 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 3, 4) ) return;
      case 4 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 4, 2) ) return;
      case 5 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 5, 4) ) return;
      case 6 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 6, 6) ) return;
      case 7 : if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 7, 2) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 0, 6) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 1, 4) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 2, 2) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 3, 4) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 4, 2) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 5, 4) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 6, 6) ) return;
               if ( check_set(pcb, po, &off, a_byte_x_2, bits, bytes, 7, 2) ) return;
   }
}


static void
do_for_bit(struct lu_calc_offs_ctx *sctx, struct prime_current_block *pcb, int bit, int thread_id)
{
   uint32_t i;
   int32_t a_byte = sctx->plist[thread_id][bit].start_a_byte;
   i = sctx->plist[thread_id][bit].index;
   struct prime_and_offset *po = &sctx->plist[thread_id][bit].primes[0];

   while (i--) {
      a_byte += po->a_byte_diff;
      compute_block(pcb, po, a_byte, bit);
      po++;
   }
}


int
lu_calc_offs_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct lu_calc_offs_ctx *sctx = ctx;
   int skip = ptx->current_block.block_num - sctx->plist[ptx->thread_index][0].last_blockno;
   int i;

   sctx->plist[ptx->thread_index][0].last_blockno = ptx->current_block.block_num;

   /* TODO: currently cannat handle skip rate > 1 */
   if (skip < 0 || skip > 1) {
      for (i = 0; i < 8; i++) {
         sctx->plist[ptx->thread_index][i].index = 0;
         sctx->plist[ptx->thread_index][i].ind_a_byte = sctx->start_prime / 30;
      }
   }

   /* Mark off multiples of 'a' in the block */
   check_new_sieve_primes_v2(sctx, &ptx->current_block, ptx->thread_index);

   do_for_bit(sctx, &ptx->current_block, 0, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 1, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 2, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 3, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 4, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 5, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 6, ptx->thread_index);
   do_for_bit(sctx, &ptx->current_block, 7, ptx->thread_index);

   return 0;
}
