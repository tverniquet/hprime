#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <immintrin.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct load_unaligned_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;

   unsigned char pbuf1_7[8*32]   __attribute__((aligned(32)));
   unsigned char pbuf1_11[12*32] __attribute__((aligned(32)));
   unsigned char pbuf1_13[14*32] __attribute__((aligned(32)));

   unsigned char pbuf2[5][2*32] __attribute__((aligned(32)));
   unsigned char pbuf3[7][4*32] __attribute__((aligned(32)));
   unsigned char pbuf4[6][6*32] __attribute__((aligned(32)));
};


#define PRIMES_1_to_16  7,11,13
#define PRIMES_16_to_32 17,19,23,29,31
#define PRIMES_32_to_64 37,41,43,47,53,59,61
#define PRIMES_64_to_96 67,71,73,79,83,89

static inline void __attribute__((always_inline))
init_early_prime_a (uint64_t A, const uint32_t a, uint64_t *buf, const uint32_t num)
{
   int i = 0;
   uint64_t base = A/a/30*30;
   uint64_t cur;
   memset(buf, 0, num * sizeof(uint64_t));

   for (;;) {
      for (i = 0; i < 8; i++) {
         cur = (base + ind_to_mod[i]) * a;

         if (cur > A + num*8*30)
            return;

         if (cur < A)
            continue;

         buf[(cur-A)/30/8] |= 1ul << (((cur - A)/30%8)*8 + ((cur-A) % 30) * 8 / 30);
      }
      base += 30;
   }
}


int
load_unaligned_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct load_unaligned_ctx *sctx = malloc(sizeof (struct load_unaligned_ctx));
   int i;
   *ctx = sctx;

   assert(end_prime <= 96);

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);


   i = 0;
   init_early_prime_a(0ul, 7, (uint64_t *)sctx->pbuf1_7, 8*4);
   init_early_prime_a(0ul, 11, (uint64_t *)sctx->pbuf1_11, 12*4);
   init_early_prime_a(0ul, 13, (uint64_t *)sctx->pbuf1_13, 14*4);

   i = 0;
#define INIT_BIN_2(PRIME) \
      init_early_prime_a(0ul, PRIME, (uint64_t *)&sctx->pbuf2[i++][0], 2*4);

   DO_FOR( INIT_BIN_2, PRIMES_16_to_32)

   i = 0;
#define INIT_BIN_3(PRIME) \
      init_early_prime_a(0ul, PRIME, (uint64_t *)&sctx->pbuf3[i++][0], 4*4);

   DO_FOR( INIT_BIN_3, PRIMES_32_to_64)

   i = 0;
#define INIT_BIN_4(PRIME) \
      init_early_prime_a(0ul, PRIME, (uint64_t *)&sctx->pbuf4[i++][0], 6*4);

   DO_FOR( INIT_BIN_4, PRIMES_64_to_96)

   return 0;
}


int
load_unaligned_free(void *ctx)
{
   struct load_unaligned_ctx *sctx = ctx;
   bzero(sctx, sizeof *sctx);
   FREE(ctx);
   return 0;
}


int
load_unaligned_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct load_unaligned_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;
   }
   return 0;
}


int
load_unaligned_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;
   (void)ctx;

   return 0;
}


/*
 * DECL_LOAD_UNALIGNED(VAL, PRIME)
 *    VAL    is used to differentiate different variable names, and is also used as the offset to perform the load (VAL * 32)
 *    PRIME  is used to refer to the table to load from, and also to calculate the mod which gives the initial offset
 *
 * NOTE: I could not get the gcc vectors to perform an unaligned load, so reverted to _mm256 intrinsics
 */
#define DECL_LOAD_UNALIGNED(VAL, PRIME) v32ui l##VAL = (v32ui)_mm256_loadu_si256((__m256i *)(ctx->pbuf1_##PRIME + VAL * 32 + (pcb->block_start_num / 30) % PRIME));
#define DECL_7(VAL)  DECL_LOAD_UNALIGNED(VAL, 7)
#define DECL_11(VAL) DECL_LOAD_UNALIGNED(VAL, 11)
#define DECL_13(VAL) DECL_LOAD_UNALIGNED(VAL, 13)

#define APPLY_SET(VAL) *out++ = l##VAL;
#define APPLY_OR(VAL) *out++ |= l##VAL;

static void
do_early_prime_group_1_7 (struct load_unaligned_ctx *ctx, struct prime_current_block *pcb)
{
   int c = CEIL_DIV(pcb->block_size, 32*7);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   DO_FOR(DECL_7,0,1,2,3,4,5,6)

   while (c--) {
      DO_FOR(APPLY_SET,0,1,2,3,4,5,6)
   }
}


static void
do_early_prime_group_1_11 (struct load_unaligned_ctx *ctx, struct prime_current_block *pcb)
{
   int c = CEIL_DIV(pcb->block_size, 32*11);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   DO_FOR(DECL_11,0,1,2,3,4,5,6,7,8,9,10)

   while (c--) {
      DO_FOR(APPLY_OR,0,1,2,3,4,5,6,7,8,9,10)
   }
}


static void
do_early_prime_group_1_13 (struct load_unaligned_ctx *ctx, struct prime_current_block *pcb)
{
   int c = CEIL_DIV(pcb->block_size, 32*13);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   DO_FOR(DECL_13,0,1,2,3,4,5,6,7,8,9,10,11,12)

   while (c--) {
      DO_FOR(APPLY_OR,0,1,2,3,4,5,6,7,8,9,10,11,12)
   }
}


static void
do_early_prime_group_2(struct load_unaligned_ctx *sctx, struct prime_current_block *pcb)
{
   int n = CEIL_DIV(pcb->block_size, 32);
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 512);
   int i = 0;

#define DEC_X_V32_B(PRIME) \
   int c##PRIME = (pcb->block_start_num / 30) % PRIME; \
   unsigned char lp_##PRIME [2*32] __attribute__((aligned(32))); \
   memcpy(lp_##PRIME, &sctx->pbuf2[i++][0], sizeof lp_##PRIME);

   DO_FOR(DEC_X_V32_B, PRIMES_16_to_32)

#define OR_PRIME(PRIME) \
      *out |= (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_##PRIME+c##PRIME)); \
      c##PRIME += (32 - PRIME); \
      c##PRIME -= c##PRIME >= PRIME ? PRIME : 0; \

   while (n--) {
      DO_FOR(OR_PRIME, PRIMES_16_to_32)
      out++;
   }
}


static void
do_early_prime_group_3(struct load_unaligned_ctx *sctx, struct prime_current_block *pcb)
{
   int n = CEIL_DIV(pcb->block_size, 64);
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 512);
   int i = 0;

#define DEC_X_V32_C(PRIME) \
   int c##PRIME = (pcb->block_start_num / 30) % PRIME; \
   unsigned char lp_##PRIME [4*32] __attribute__((aligned(32))); \
   memcpy(lp_##PRIME, &sctx->pbuf3[i++][0], sizeof lp_##PRIME);

   DO_FOR(DEC_X_V32_C, PRIMES_32_to_64)

#define OR_PRIME_C(PRIME) \
      *out        |= (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_##PRIME + c##PRIME)); \
      *(out + 1)  |= (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_##PRIME + c##PRIME + 32)); \
      c##PRIME += (64 - PRIME); \
      c##PRIME -= c##PRIME >= PRIME ? PRIME : 0; \

   while (n--) {
      DO_FOR(OR_PRIME_C, PRIMES_32_to_64)
      out+=2;
   }
}


static void
do_early_prime_group_4(struct load_unaligned_ctx *sctx, struct prime_current_block *pcb)
{
   int n = CEIL_DIV(pcb->block_size,96);
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 512);
   int i = 0;

#define DEC_X_V32_D(PRIME) \
   int c##PRIME = (pcb->block_start_num / 30) % PRIME; \
   unsigned char lp_##PRIME [6*32] __attribute__((aligned(32))); \
   memcpy(lp_##PRIME, &sctx->pbuf4[i++][0], sizeof lp_##PRIME);

   DO_FOR(DEC_X_V32_D, PRIMES_64_to_96)

#define OR_PRIME_D(PRIME) \
      *out        |= (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_##PRIME + c##PRIME)); \
      *(out + 1)  |= (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_##PRIME + c##PRIME + 32)); \
      *(out + 2)  |= (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_##PRIME + c##PRIME + 64)); \
      c##PRIME += (96 - PRIME); \
      c##PRIME -= c##PRIME >= PRIME ? PRIME : 0; \

   while (n--) {
      DO_FOR(OR_PRIME_D, PRIMES_64_to_96)
      out+=3;
   }
}



int
load_unaligned_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct load_unaligned_ctx *sctx = ctx;

   /* Note the conditional checks means that another method be used to, say, calculate primes 32-64 */

   if (sctx->start_prime <=7 && sctx->end_prime >= 7)
      do_early_prime_group_1_7 (sctx, &ptx->current_block);

   if (sctx->start_prime <= 11 && sctx->end_prime >= 11)
      do_early_prime_group_1_11 (sctx, &ptx->current_block);

   if (sctx->start_prime <= 13 && sctx->end_prime >= 13)
      do_early_prime_group_1_13 (sctx, &ptx->current_block);

   if (sctx->start_prime < 32 && sctx->end_prime >= 17)
      do_early_prime_group_2(sctx, &ptx->current_block);

   if (sctx->start_prime < 64 && sctx->end_prime > 32)
      do_early_prime_group_3(sctx, &ptx->current_block);

   if (sctx->start_prime < 96 && sctx->end_prime > 64)
      do_early_prime_group_4(sctx, &ptx->current_block);

   return 0;
}
