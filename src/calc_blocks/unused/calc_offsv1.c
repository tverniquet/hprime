#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>
#include <avx2intrin.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

#if __AVX2__
#define BPV   32
#define WPV   16

#define FV v32_16ui

#define DECLV(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P) {A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P}
#else
#define BPV   16
#define WPV    8

#define FV v16_8ui

#define DECLV(A,B,C,D,E,F,G,H,I,J,K,L,M,N,O,P) {A,B,C,D,E,F,G,H}
#endif

/*
 * "Stuff" contains items grouped into 16 (2-byte) ints so they can be read
 * directly into 32 byte vectors.
 *
 * |abyte1|abyte2|...|abyte15|  n1  |  n2  |...|  n15  |prime1|prime2|...|prime15|
 * |-------------------------|-------------------------|--------------------------|
 *         v0                        v1                        v2
 *
 * NOTE: Actually only abyte or prime is needed since these can be converted to
 * the others. This would drop the amount of storage from:
 *  3 * 2 * ~3500 primes = 21K
 *  1 * 2 * ~3500 primes =  7K
 *
 *  However, currently it is faster to store these directly.
 */
struct prime_list
{
   unsigned char *stuff;
   int      stuff_c;
};

struct prime_offs {
   unsigned char *offs[8];
   int            offs_i[8];
};


struct calc_offsv1_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   struct prime_list primes[8];
   struct prime_offs *thread_offs;
   int nthreads;
};


int
calc_offsv1_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct calc_offsv1_ctx *sctx = calloc(sizeof (struct calc_offsv1_ctx), 1);
   int i;
   int j;
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime <= pctx->current_block->block_size && sctx->end_prime <= (1<<15));

   sctx->nthreads = pctx->num_threads ?:1;

   sctx->thread_offs = calloc(sizeof *sctx->thread_offs, sctx->nthreads);

   for (i = 0; i < 8; i++) {
      sctx->primes[i].stuff = aligned_alloc(32, (3*sizeof(uint16_t)) * 3500);
   }
   for (j = 0; j < sctx->nthreads; j++) {
      for (i = 0; i < 8; i++) {
         sctx->thread_offs[j].offs[i] = aligned_alloc(32, (sizeof(uint16_t) * 3500));
      }
   }

   return 0;
}


int
calc_offsv1_free(void *ctx)
{
   struct calc_offsv1_ctx *sctx = ctx;
   int i;
   int j;
   for (i = 0; i < 8; i++) {
      free(sctx->primes[i].stuff);
   }
   for (j = 0; j < sctx->nthreads; j++) {
      for (i = 0; i < 8; i++) {
         free(sctx->thread_offs[j].offs[i]);
      }
   }
   free(sctx->thread_offs);
   FREE(ctx);
   return 0;
}


int
calc_offsv1_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct calc_offsv1_ctx *sctx = ctx;
   struct prime_list *pl;

   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      pl = &sctx->primes[pp_to_bit(primelist[*ind])];

      *((uint16_t *)pl->stuff + (pl->stuff_c % WPV) + (pl->stuff_c/WPV*3*WPV)) = num_to_bytes(primelist[*ind]);
      *((uint16_t *)pl->stuff + (pl->stuff_c % WPV) + (pl->stuff_c/WPV*3*WPV)+WPV) = 32*1024/primelist[*ind];
      *((uint16_t *)pl->stuff + (pl->stuff_c % WPV) + (pl->stuff_c/WPV*3*WPV)+2*WPV) = primelist[*ind];

      pl->stuff_c++;
   }
   return 0;
}


int
calc_offsv1_skip_to(struct prime_thread_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct calc_offsv1_ctx *sctx = ctx;
   /* Recalculate all of the offsets depending on the new block */
   int i;
   for (i = 0; i < 8; i++)
      sctx->thread_offs[pctx->thread_index].offs_i[i] = 0;
   return 0;
}


/*
 * The first time we don't know how many times to mark off the prime in the
 * block (ie don't really know 'n'). Anyway, keeping this separate to the
 * calc_block function allows the "calc_block" function to make more
 * assumptions and optimise better.
 */
static void
compute_block_first_time(struct prime_current_block *pcb, uint32_t sieve_prime)
{
   char *bms[8];
   int i;
   const unsigned char *bits = a_x_b_bitmask[pp_to_bit(sieve_prime)];
   int bit;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;
   uint16_t offsets[8];
   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);

   int done = 0;

   multiplier_base = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));

   for (bit = 0; bit < 8; bit++) {
      multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * sieve_prime);
      while (multiplier_byte < block_start_byte)
         multiplier_byte += sieve_prime;

      while (multiplier_byte - block_start_byte >= UINT16_MAX)
         multiplier_byte -= sieve_prime;

      offsets[bit] = multiplier_byte - block_start_byte;
   }

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
}


static void
check_new_sieve_primes_a(struct calc_offsv1_ctx *sctx, struct prime_current_block *pcb, struct prime_offs *offs)
{
   uint32_t sieve_prime;

   int i;
   struct prime_list *pl;
   for (i = 0; i < 8; i++) {
      pl = &sctx->primes[i];
      for ( ; offs->offs_i[i] < pl->stuff_c; offs->offs_i[i]++) {
         sieve_prime = *((uint16_t *)pl->stuff + 2*WPV + (offs->offs_i[i]%WPV) + (offs->offs_i[i]/WPV * 3*WPV));

         if (sieve_prime > pcb->sqrt_end_num)
            break;

         compute_block_first_time(pcb, sieve_prime);
         *((int16_t *)offs->offs[i] + offs->offs_i[i]) = sieve_prime - ((num_to_bytes(pcb->block_start_num) + sctx->nthreads * pcb->block_size) % sieve_prime);
      }
   }
}


static inline void __attribute__((always_inline))
compute_block_offs_low(struct prime_current_block *pcb __attribute__((unused)), uint32_t sieve_prime, char *bms[8], int n)
{
   const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   int i;

   for (i = 0; i < 8; i++)
      *bms[i] |= bits[i];

   while (n--)
      for (i = 0; i < 8; i++)
         *(bms[i] += sieve_prime) |= bits[i];
}


static inline void __attribute__((always_inline))
get_offs_a(struct calc_offsv1_ctx *ctx __attribute__((unused)), FV *ap, FV *pp, FV *np, FV *op, FV *curoffs, const int bit, int nthreads)
{
   FV a_bytes = *ap;
   FV ns = *np;
   FV primes = *pp;
   FV v_offsets_0;
   FV v_offsets_1;
   FV v_offsets_2;
   FV v_offsets_3;
   FV v_offsets_4;
   FV v_offsets_5;
   FV v_offsets_6;
   FV v_offsets_7;
   FV ones = DECLV(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1);

   /* XXX hardcoded for now */
   FV bs = DECLV(32*1024,32*1024,32*1024,32*1024,32*1024,32*1024,32*1024,32*1024,
                  32*1024,32*1024,32*1024,32*1024,32*1024,32*1024,32*1024,32*1024);

   switch(bit) {
      case 0: {
         v_offsets_0 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_1 = v_offsets_0 + a_bytes + a_bytes + a_bytes;
         v_offsets_2 = v_offsets_1 + a_bytes + a_bytes;
         v_offsets_3 = v_offsets_2 + a_bytes;
         v_offsets_4 = v_offsets_3 + a_bytes + a_bytes;
         v_offsets_5 = v_offsets_4 + a_bytes;
         v_offsets_6 = v_offsets_5 + a_bytes + a_bytes;
         v_offsets_7 = v_offsets_6 + a_bytes + a_bytes + a_bytes;
      }
      break;
      case 1: {
         v_offsets_1 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_5 = v_offsets_1 + a_bytes + a_bytes + a_bytes + ones;
         v_offsets_4 = v_offsets_5 + a_bytes + a_bytes + ones;
         v_offsets_0 = v_offsets_4 + a_bytes + ones;
         v_offsets_7 = v_offsets_0 + a_bytes + a_bytes;
         v_offsets_3 = v_offsets_7 + a_bytes + ones;
         v_offsets_2 = v_offsets_3 + a_bytes + a_bytes + ones;
         v_offsets_6 = v_offsets_2 + a_bytes + a_bytes + a_bytes + ones;
      }
      break;
      case 2: {
         FV twos = DECLV(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);
         v_offsets_2 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_4 = v_offsets_2 + a_bytes + a_bytes + a_bytes + twos;
         v_offsets_0 = v_offsets_4 + a_bytes + a_bytes + twos;
         v_offsets_6 = v_offsets_0 + a_bytes;
         v_offsets_1 = v_offsets_6 + a_bytes + a_bytes + twos;
         v_offsets_7 = v_offsets_1 + a_bytes;
         v_offsets_3 = v_offsets_7 + a_bytes + a_bytes + twos;
         v_offsets_5 = v_offsets_3 + a_bytes + a_bytes + a_bytes + twos;
      }
      break;
      case 3: {
         FV threes = DECLV(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3);
         v_offsets_3 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_0 = v_offsets_3 + a_bytes + a_bytes + a_bytes + threes;
         v_offsets_6 = v_offsets_0 + a_bytes + a_bytes + ones;
         v_offsets_5 = v_offsets_6 + a_bytes + ones;
         v_offsets_2 = v_offsets_5 + a_bytes + a_bytes + ones + ones;
         v_offsets_1 = v_offsets_2 + a_bytes + ones;
         v_offsets_7 = v_offsets_1 + a_bytes + a_bytes + ones;
         v_offsets_4 = v_offsets_7 + a_bytes + a_bytes + a_bytes + threes;
      }
      break;
      case 4: {
         FV threes = DECLV(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3);
         v_offsets_4 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_7 = v_offsets_4 + a_bytes + a_bytes + a_bytes + threes;
         v_offsets_1 = v_offsets_7 + a_bytes + a_bytes + threes;
         v_offsets_2 = v_offsets_1 + a_bytes + ones;
         v_offsets_5 = v_offsets_2 + a_bytes + a_bytes + ones + ones;
         v_offsets_6 = v_offsets_5 + a_bytes + ones;
         v_offsets_0 = v_offsets_6 + a_bytes + a_bytes + threes;
         v_offsets_3 = v_offsets_0 + a_bytes + a_bytes + a_bytes + threes;
      }
      break;
      case 5: {
         FV twos = DECLV(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);
         v_offsets_5 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_3 = v_offsets_5 + a_bytes + a_bytes + a_bytes + twos + twos;
         v_offsets_7 = v_offsets_3 + a_bytes + a_bytes + twos;
         v_offsets_1 = v_offsets_7 + a_bytes + twos;
         v_offsets_6 = v_offsets_1 + a_bytes + a_bytes + twos;
         v_offsets_0 = v_offsets_6 + a_bytes + twos;
         v_offsets_4 = v_offsets_0 + a_bytes + a_bytes + twos;
         v_offsets_2 = v_offsets_4 + a_bytes + a_bytes + a_bytes + twos + twos;
      }
      break;
      case 6: {
         FV threes = DECLV(3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3);
         FV twos = DECLV(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);
         v_offsets_6 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_2 = v_offsets_6 + a_bytes + a_bytes + a_bytes + threes + twos;
         v_offsets_3 = v_offsets_2 + a_bytes + a_bytes + threes;
         v_offsets_7 = v_offsets_3 + a_bytes + ones;
         v_offsets_0 = v_offsets_7 + a_bytes + a_bytes + twos + twos;
         v_offsets_4 = v_offsets_0 + a_bytes + ones;
         v_offsets_5 = v_offsets_4 + a_bytes + a_bytes + threes;
         v_offsets_1 = v_offsets_5 + a_bytes + a_bytes + a_bytes + threes + twos;
      }
      break;
      case 7: {
         FV twos = DECLV(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2);
         v_offsets_7 = *op + a_bytes;
         a_bytes <<= ones;
         v_offsets_6 = v_offsets_7 + a_bytes + a_bytes + a_bytes + twos + twos + twos;
         v_offsets_5 = v_offsets_6 + a_bytes + a_bytes + twos + twos;
         v_offsets_4 = v_offsets_5 + a_bytes + twos;
         v_offsets_3 = v_offsets_4 + a_bytes + a_bytes + twos + twos;
         v_offsets_2 = v_offsets_3 + a_bytes + twos;
         v_offsets_1 = v_offsets_2 + a_bytes + a_bytes + twos + twos;
         v_offsets_0 = v_offsets_1 + a_bytes + a_bytes + a_bytes + twos + twos + twos;
      }
      break;
      default:
      /* never reach but stop compiler warning */
         return;
   }

   /* TODO */
   switch(nthreads) {
      case 8:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 7:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 6:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 5:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 4:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 3:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 2:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
      case 1:
      case 0:
         *op += ns * primes - bs;
         *op += (*op > bs) & primes;
   }

   v_offsets_0 -= (v_offsets_0 >= primes) & primes;
   v_offsets_1 -= (v_offsets_1 >= primes) & primes;
   v_offsets_2 -= (v_offsets_2 >= primes) & primes;
   v_offsets_3 -= (v_offsets_3 >= primes) & primes;
   v_offsets_4 -= (v_offsets_4 >= primes) & primes;
   v_offsets_5 -= (v_offsets_5 >= primes) & primes;
   v_offsets_6 -= (v_offsets_6 >= primes) & primes;
   v_offsets_7 -= (v_offsets_7 >= primes) & primes;

   curoffs[0] = v_offsets_0;
   curoffs[1] = v_offsets_1;
   curoffs[2] = v_offsets_2;
   curoffs[3] = v_offsets_3;
   curoffs[4] = v_offsets_4;
   curoffs[5] = v_offsets_5;
   curoffs[6] = v_offsets_6;
   curoffs[7] = v_offsets_7;
}


static void
do_bit_a(struct calc_offsv1_ctx *ctx, struct prime_current_block *pcb, int nthreads, struct prime_offs *po, const int bit)
{
   unsigned char buf[BPV*8] __attribute__((aligned(BPV)));
   FV *ap = (FV *)ctx->primes[bit].stuff;
   FV *op = (FV *)po->offs[bit];
   uint16_t *offs;

   char *bms[8];
   uint16_t *primep = (uint16_t *)ctx->primes[bit].stuff + WPV*2;
   uint16_t *np = (uint16_t *)ctx->primes[bit].stuff + WPV;

   int k = po->offs_i[bit]/WPV + 1;
   int i = po->offs_i[bit] % WPV + 1;
   int i0;
   int i1;

   if (po->offs_i[bit] == 0)
      return;


   while (k--) {

      if (k == 0 && i-- == 0)
         return;

      get_offs_a(ctx, ap, ap+2, ap+1, op, (FV *)buf, bit, nthreads);
      offs = (uint16_t *)buf;

      /*
       * For primes > 16k, overrunning the buffer is more likely to
       * result in a page fault so the idea of treating this separately
       * was to reduce the cache misses, which it does, however it also
       * increases the number of instructions.. Overall it doesn't
       * have much effect.
       */
      if (np[0] == 1) {
         const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
         uint32_t sieve_prime;
         int d;
         for (i0 = 0; i0 < WPV; i0++, offs++) {
            if (k == 0 && i-- == 0)
               return;
            sieve_prime = *primep++;
            d = pcb->block_size - sieve_prime;
            for (i1 = 0; i1 < 8; i1++) {
               bms[i1] = pcb->block + *(offs + i1*WPV);
               *bms[i1] |= bits[i1];
               bms[i1] += *(offs + i1*WPV) > d ? 0 : sieve_prime;
               *bms[i1] |= bits[i1];
            }
         }
         np += WPV;
      }
      else {
         for (i0 = 0; i0 < WPV; i0++, offs++) {
            if (k == 0 && i-- == 0)
               return;
            for (i1 = 0; i1 < 8; i1++)
               bms[i1] = pcb->block + *(offs + i1*WPV);
            compute_block_offs_low(pcb, (uint32_t)*primep++, bms, *np++);
         }
      }


      op++;
      ap += 3;
      primep += 2*WPV;
      np += 2*WPV;
   }
}


int
calc_offsv1_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct calc_offsv1_ctx *sctx = ctx;

   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 0);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 1);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 2);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 3);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 4);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 5);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 6);
   do_bit_a(sctx, &ptx->current_block, sctx->nthreads, &sctx->thread_offs[ptx->thread_index], 7);

   check_new_sieve_primes_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index]);

   return 0;
}
