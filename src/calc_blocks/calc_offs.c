#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
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
   int            last_blockno;
};


struct calc_offs_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   struct prime_list primes[8];
   struct prime_offs *thread_offs;
   int nthreads;
   int blocks_per_run;
};


int
calc_offs_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct calc_offs_ctx *sctx = calloc(sizeof (struct calc_offs_ctx), 1);
   int i;
   int j;
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->end_prime <= pctx->current_block->block_size && sctx->end_prime <= (1<<15));

   sctx->nthreads = pctx->num_threads ?:1;
   sctx->blocks_per_run = pctx->blocks_per_run ?:1;

   sctx->thread_offs = calloc(sizeof *sctx->thread_offs, sctx->nthreads);

   for (i = 0; i < 8; i++) {
      sctx->primes[i].stuff = aligned_alloc(32, (3*sizeof(uint16_t)) * 3500);
   }
   for (j = 0; j < sctx->nthreads; j++) {
      for (i = 0; i < 8; i++) {
         sctx->thread_offs[j].offs[i] = aligned_alloc(32, (sizeof(uint16_t) * 3500));
      }
      sctx->thread_offs[j].last_blockno = INT32_MAX;
   }

   return 0;
}


int
calc_offs_free(void *ctx)
{
   struct calc_offs_ctx *sctx = ctx;
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
calc_offs_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct calc_offs_ctx *sctx = ctx;
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
calc_offs_skip_to(struct prime_thread_ctx *pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   struct calc_offs_ctx *sctx = ctx;
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
   char *bmp;
   int   i;
   int offsets[8];
   const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

   for (i = 0; i < 8; i++)
      offsets[sieve_prime * ind_to_mod[i] % 30 * 8 / 30] = (sieve_prime * ind_to_mod[i]) / 30;

   bmp = pcb->block + (uint64_t)sieve_prime * (sieve_prime / 30) - pcb->block_start_byte;

   for (; bmp < pcb->block + pcb->block_size; bmp += sieve_prime)
      for (i = 0; i < 8; i++)
         *(bmp + offsets[i]) |= bits[i];
}


static void
check_new_sieve_primes_a(struct calc_offs_ctx *sctx, struct prime_current_block *pcb, struct prime_offs *offs, int skip, int mode)
{
   uint32_t sieve_prime;

   int i;
   struct prime_list *pl;
   for (i = 0; i < 8; i++) {
      pl = &sctx->primes[i];
      for ( ; offs->offs_i[i] < pl->stuff_c; offs->offs_i[i]++) {
         sieve_prime = *((uint16_t *)pl->stuff + 2*WPV + (offs->offs_i[i]%WPV) + (offs->offs_i[i]/WPV * 3*WPV));

         if (mode == 0) {
            if (sieve_prime * sieve_prime >= pcb->block_start_num)
               break;

            *((int16_t *)offs->offs[i] + offs->offs_i[i]) = sieve_prime - (pcb->block_start_byte - skip * pcb->block_size) % sieve_prime;
         }
         else {
            if (sieve_prime > pcb->sqrt_end_num)
               break;
            compute_block_first_time(pcb, sieve_prime);
            *((int16_t *)offs->offs[i] + offs->offs_i[i]) = sieve_prime - pcb->block_start_byte % sieve_prime;
         }
      }
   }
}


static inline void __attribute__((always_inline))
adjust_offsets(FV *curoffs, FV *ns, FV *primes, FV *bs, FV *o, int ind)
{
   *o -= (*o >= *primes) & *primes;
   curoffs[ind] = *o;
   *o += *ns * *primes;
   *o -= (*o >= *bs) & *primes;
   curoffs[8 + ind] = *o;
}


static inline void __attribute__((always_inline))
get_offs_a(struct calc_offs_ctx *ctx __attribute__((unused)), FV *ap, FV *pp, FV *np, FV *op, FV *curoffs, int skip, const int bit)
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

   /*
    * Calculate (and store) the new offsets based on how many blocks have been
    * skipped
    */
   while (skip--) {
      *op -= bs - ns * primes;
      *op += (*op > bs) & primes;  /* unsigned so this detects the wrap */
   }


   /*
    * The idea is that each section can be inlined in a new version of this
    * function.
    */
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


   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_0, 0);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_1, 1);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_2, 2);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_3, 3);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_4, 4);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_5, 5);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_6, 6);
   adjust_offsets(curoffs, &ns, &primes, &bs, &v_offsets_7, 7);
}


static inline void __attribute__((always_inline))
mark_16_a (struct prime_current_block *pcb, uint16_t *offs_orig, uint16_t *primes, uint16_t *ns, int nconsume)
{
   int i;
   int j;
   int n;
   int a;
   char *bmp;
   uint16_t *offs = offs_orig;
   const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   for (i = 0; i < nconsume; i++, offs++) {
      n = ns[i] + 1;
      a = primes[i];
      for (bmp = pcb->block; --n; bmp += a)
         for (j = 0; j < 8; j++)
            *(bmp + *(offs + j*WPV)) |= bits[j];
   }
   offs = offs_orig + WPV*8;
   bmp = pcb->block;
   for (i = 0; i < nconsume; i++, offs++)
      for (j = 0; j < 8; j++)
         *(bmp + *(offs + j*WPV)) |= bits[j];
}


/*
 * using ((always_inline)) OR ((no_inline)) seems to causes it to slow down
 *
 * I think that the switch inside of get_offs_a is promoted outside of the --k loop.
 *
 */
static void
do_bit_a(struct calc_offs_ctx *ctx, struct prime_current_block *pcb, struct prime_offs *po, int skip, const int bit)
{
   unsigned char buf[BPV*16] __attribute__((aligned(BPV)));
   FV *ap = (FV *)ctx->primes[bit].stuff;
   FV *op = (FV *)po->offs[bit];
   uint16_t *offs;

   uint16_t *primep = (uint16_t *)ctx->primes[bit].stuff + WPV*2;
   uint16_t *np = (uint16_t *)ctx->primes[bit].stuff + WPV;

   int k = po->offs_i[bit]/WPV + 1;
   int i = po->offs_i[bit] % WPV;

   if (po->offs_i[bit] == 0)
      return;


   while (--k) {

      get_offs_a(ctx, ap, ap+2, ap+1, op, (FV *)buf, skip, bit);
      offs = (uint16_t *)buf;

      /* This switch seems to help the branch prediction */
      switch(np[0]) {
         case 1:  mark_16_a (pcb, offs, primep, np, WPV); break;
         case 2:  mark_16_a (pcb, offs, primep, np, WPV); break;
         case 3:  mark_16_a (pcb, offs, primep, np, WPV); break;
         case 4:  mark_16_a (pcb, offs, primep, np, WPV); break;
         default: mark_16_a (pcb, offs, primep, np, WPV); break;
      }
      primep += WPV;
      np += WPV;

      op++;
      ap += 3;
      primep += 2*WPV;
      np += 2*WPV;
   }

   if (i) {
      get_offs_a(ctx, ap, ap+2, ap+1, op, (FV *)buf, skip, bit);
      offs = (uint16_t *)buf;
      mark_16_a (pcb, offs, primep, np, i);
   }
}


int
calc_offs_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct calc_offs_ctx *sctx = ctx;
   struct prime_offs *po = &sctx->thread_offs[ptx->thread_index];
   int skip = ptx->current_block.block_num - po->last_blockno;

   po->last_blockno = ptx->current_block.block_num;

   if (skip < 0 || skip > 8) {
      memset(po->offs_i, 0, sizeof po->offs_i);
      skip = 0;
   }

   check_new_sieve_primes_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 0);

   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 0);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 1);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 2);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 3);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 4);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 5);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 6);
   do_bit_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 7);

   check_new_sieve_primes_a(sctx, &ptx->current_block, &sctx->thread_offs[ptx->thread_index], skip, 1);

   return 0;
}
