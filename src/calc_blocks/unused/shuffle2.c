#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <immintrin.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct shuffle2_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint64_t *low_prime_buf;
   uint64_t lp_17[8];
   uint64_t lp_19[8];
   uint64_t lp_23[8];
   uint64_t lp_29[8];
   uint64_t lp_31[8];
   uint64_t lp_37[16];
   uint64_t lp_41[16];
   uint64_t lp_43[16];
   uint64_t lp_47[16];
};


static inline void __attribute__((always_inline)) init_early_prime_a (uint64_t A, const uint32_t a, uint64_t *buf, const uint32_t num);


int
shuffle2_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct shuffle2_ctx *sctx = malloc(sizeof (struct shuffle2_ctx));
   *ctx = sctx;

   assert(end_prime < 64);
   assert(start_prime < end_prime);

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);
   sctx->low_prime_buf = aligned_alloc(32, 100*4*sizeof(uint64_t));
   init_early_prime_a(0ul, 17, sctx->lp_17, 8);
   init_early_prime_a(0ul, 19, sctx->lp_19, 8);
   init_early_prime_a(0ul, 23, sctx->lp_23, 8);
   init_early_prime_a(0ul, 29, sctx->lp_29, 8);
   init_early_prime_a(0ul, 31, sctx->lp_31, 8);

   init_early_prime_a(0ul, 37, sctx->lp_37, 16);
   init_early_prime_a(0ul, 41, sctx->lp_41, 16);
   init_early_prime_a(0ul, 43, sctx->lp_43, 16);
   init_early_prime_a(0ul, 47, sctx->lp_47, 16);
   return 0;
}


int
shuffle2_free(void *ctx)
{
   struct shuffle2_ctx *sctx = ctx;
   free(sctx->low_prime_buf);
   bzero(sctx, sizeof *sctx);
   FREE(ctx);
   return 0;
}


int
shuffle2_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct shuffle2_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;
   }
   return 0;
}


int
shuffle2_skip_to(struct prime_thread_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   return 0;
}


#define CEIL_i64_A(PRIME) ( (PRIME + 7)/8 <= 2 ? 2 \
                           :(PRIME + 7)/8 <= 4 ? 4 \
                           : 8)


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


static inline void __attribute__((always_inline))
shift_a (const uint32_t a, uint64_t *arr)
{
   const int numa = ((a + 7) / 8);     /* The ceil of the number of i64 per x 30 cycle */
   const int num  = ((a + 7) / 8);     /*CEIL_i64_A(a); */    /* Further adjusted to be a power of 2 */
   const int i64shift = num - numa;    /* Number of i64 of padding */
   const int byteshift = numa * 8 - a; /* After each x 30 cycle, the number of bytes left over */
   int i;

   for (i = 0; i < num; i++)
      arr[i] = (arr[(i + i64shift) % num] >> (byteshift*8)) | (arr[(i + i64shift + 1) % num] << (64 - (byteshift*8)));
}




#define DOALL_1(MACRO, ...) MACRO(7, 1, __VA_ARGS__)
#define DOALL_2(MACRO, ...) MACRO(11, 2, __VA_ARGS__) MACRO(13, 2, __VA_ARGS__)
#define DOALL_4(MACRO, ...) MACRO(17, 4, __VA_ARGS__) MACRO(19, 4, __VA_ARGS__) MACRO(23, 4, __VA_ARGS__) MACRO(29, 4, __VA_ARGS__) MACRO(31, 4, __VA_ARGS__)
#define DOALL_8(MACRO, ...) MACRO(37, 8, __VA_ARGS__) MACRO(41, 8, __VA_ARGS__) MACRO(43, 8, __VA_ARGS__) MACRO(47, 8, __VA_ARGS__) MACRO(53, 8, __VA_ARGS__) MACRO(59, 8, __VA_ARGS__) MACRO(61, 8, __VA_ARGS__)

#define DO_ALL(MACRO, ...) DOALL_2(MACRO, __VA_ARGS__) DOALL_4(MACRO, __VA_ARGS__) DOALL_8(MACRO, __VA_ARGS__)


#define ALL_1 7
#define ALL_2 11,13
#define ALL_3 17,19,23
#define ALL_4 29,31
#define ALL_5 37
#define ALL_6 41,43,47
#define ALL_7 53
#define ALL_8 59,61
#define ALL   ALL_1 , ALL_2 , ALL_3 , ALL_4 , ALL_5 , ALL_6 , ALL_7 , ALL_8


#define DO_FOR_GROUP(X, MACRO) DO_FOR(MACRO, ALL_##X)
#define DO_FOR_ALL(MACRO) DO_FOR(MACRO, ALL)

#define DEC_X_V16(PRIME) unsigned char lp_##PRIME [8*CEIL_i64_A(PRIME)] __attribute__((aligned(16))); v16ui *lp##PRIME = (v16ui *)lp_##PRIME;
#define DEC_X_V32(PRIME) unsigned char lp_##PRIME [8*2*CEIL_i64_A(PRIME)] __attribute__((aligned(32))); v32ui *lp##PRIME = (v32ui *)lp_##PRIME;

#define DEC_X_V32_B(PRIME) unsigned char lp_##PRIME [8*8] __attribute__((aligned(32))); v32ui l##PRIME; int c##PRIME = 0;

#define DEC_X_V32_C(PRIME) \
   unsigned char lp_##PRIME [16*8] __attribute__((aligned(32))); \
   v32ui    l##PRIME##_1 __attribute__((aligned(32)));\
   v32ui    l##PRIME##_2 __attribute__((aligned(32)));\
   int      c##PRIME = 0;

#define DEC_X_V32_F(PRIME) unsigned char lp_##PRIME [8*8*2] __attribute__((aligned(32))); v32ui *lp##PRIME = (v32ui *)lp_##PRIME;
#define INIT_X_32_F(PRIME) init_early_prime_a(pcb->block_start_num, PRIME, (uint64_t *)lp_##PRIME, 2*8);
#define DEC_X(PRIME) uint64_t lp_##PRIME [CEIL_i64_A(PRIME)];

#define INIT_X(PRIME) init_early_prime_a(pcb->block_start_num, PRIME, (uint64_t *)lp_##PRIME, CEIL_i64_A(PRIME));

#define INIT_X_32(PRIME) init_early_prime_a(pcb->block_start_num, PRIME, (uint64_t *)lp_##PRIME, 2*8);

#define OR_X(PRIME)     |(lp_##PRIME[i % (PRIME + 7 / 8)])
#define SHIFT_X(PRIME) shift_a(PRIME, lp_##PRIME);



static inline void __attribute__((always_inline))
do_early_prime_group_1_7 (struct shuffle2_ctx *ctx, struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8 + (4*7 - 1))/(4*7);
   v32ui *lp = (v32ui *)__builtin_assume_aligned(ctx->low_prime_buf, 32);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   register v32ui l0;
   register v32ui l1;
   register v32ui l2;
   register v32ui l3;
   register v32ui l4;
   register v32ui l5;
   register v32ui l6;

   init_early_prime_a(pcb->block_start_num, 7, ctx->low_prime_buf, 4*7);

   l0 = lp[0];
   l1 = lp[1];
   l2 = lp[2];
   l3 = lp[3];
   l4 = lp[4];
   l5 = lp[5];
   l6 = lp[6];

   while (c--) {
      *out++ = l0;
      *out++ = l1;
      *out++ = l2;
      *out++ = l3;
      *out++ = l4;
      *out++ = l5;
      *out++ = l6;
   }
}


static inline void __attribute__((always_inline))
do_early_prime_group_1_11 (struct shuffle2_ctx *ctx, struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8 + (4*11 - 1))/(4*11);
   v32ui *lp = (v32ui *)__builtin_assume_aligned(ctx->low_prime_buf, 32);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   register v32ui l0;
   register v32ui l1;
   register v32ui l2;
   register v32ui l3;
   register v32ui l4;
   register v32ui l5;
   register v32ui l6;
   register v32ui l7;
   register v32ui l8;
   register v32ui l9;
   register v32ui l10;

   init_early_prime_a(pcb->block_start_num, 11, ctx->low_prime_buf, 4*11);

   l0 = lp[0];
   l1 = lp[1];
   l2 = lp[2];
   l3 = lp[3];
   l4 = lp[4];
   l5 = lp[5];
   l6 = lp[6];
   l7 = lp[7];
   l8 = lp[8];
   l9 = lp[9];
   l10 = lp[10];

   while (c--) {
      *out++ |= l0;
      *out++ |= l1;
      *out++ |= l2;
      *out++ |= l3;
      *out++ |= l4;
      *out++ |= l5;
      *out++ |= l6;
      *out++ |= l7;
      *out++ |= l8;
      *out++ |= l9;
      *out++ |= l10;
   }
}


static inline void __attribute__((always_inline))
do_early_prime_group_1_13 (struct shuffle2_ctx *ctx, struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8 + (4*13 - 1))/(4*13);
   v32ui *lp = (v32ui *)__builtin_assume_aligned(ctx->low_prime_buf, 32);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   register v32ui l0;
   register v32ui l1;
   register v32ui l2;
   register v32ui l3;
   register v32ui l4;
   register v32ui l5;
   register v32ui l6;
   register v32ui l7;
   register v32ui l8;
   register v32ui l9;
   register v32ui l10;
   register v32ui l11;
   register v32ui l12;

   init_early_prime_a(pcb->block_start_num, 13, ctx->low_prime_buf, 4*13);

   l0 = lp[0];
   l1 = lp[1];
   l2 = lp[2];
   l3 = lp[3];
   l4 = lp[4];
   l5 = lp[5];
   l6 = lp[6];
   l7 = lp[7];
   l8 = lp[8];
   l9 = lp[9];
   l10 = lp[10];
   l11 = lp[11];
   l12 = lp[12];

   while (c--) {
      *out++ |= l0;
      *out++ |= l1;
      *out++ |= l2;
      *out++ |= l3;
      *out++ |= l4;
      *out++ |= l5;
      *out++ |= l6;
      *out++ |= l7;
      *out++ |= l8;
      *out++ |= l9;
      *out++ |= l10;
      *out++ |= l11;
      *out++ |= l12;
   }
}


static inline void __attribute__((always_inline))
do_early_prime_group_1_type (struct shuffle2_ctx *ctx, struct prime_current_block *pcb, const int a, const int mode)
{
   int c = (pcb->block_size/8 + (4*a - 1))/(4*a);
   int i;
   v32ui *lp = (v32ui *)__builtin_assume_aligned(ctx->low_prime_buf, 32);
   v32ui *out =  (v32ui *)__builtin_assume_aligned(pcb->block,  512);

   init_early_prime_a(pcb->block_start_num, a, ctx->low_prime_buf, 4*a);
   while (c--) {
      for (i = 0; i < a; i++) {
         if (mode)
            *out++ = lp[i];
         else
            *out++ |= lp[i];
      }
   }
}


/* seem to go better inside a function? */
#if 1
static void do_early_prime_7(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_7(ctx, pcb); }
static void do_early_prime_11(struct shuffle2_ctx * ctx, struct prime_current_block *pcb) { do_early_prime_group_1_11(ctx, pcb); }
static void do_early_prime_13(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_13(ctx, pcb); }

/* Only quicker doing individually for the others all the ymms fit in registers */
/*static void do_early_prime_17(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_type(ctx, pcb, 17, 0); }*/
/*static void do_early_prime_19(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_type(ctx, pcb, 19, 0); }*/
/*static void do_early_prime_23(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_type(ctx, pcb, 23, 0); }*/
/*static void do_early_prime_29(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_type(ctx, pcb, 19, 0); }*/
/*static void do_early_prime_31(struct shuffle2_ctx *ctx, struct prime_current_block *pcb) { do_early_prime_group_1_type(ctx, pcb, 31, 0); }*/
#else
static void
do_early_prime_group_1(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8)/2;
   int i;
   uint64_t *bm = (uint64_t *)pcb->block;
   v16ui *out = (v16ui *)pcb->block;

   DO_FOR(DEC_X, ALL_1, ALL_2)
   DO_FOR(INIT_X, ALL_1, ALL_2)

   v16ui l7_shuffle = { 2 , 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1 };
   v16ui l11_shuffle = { 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4 };
   v16ui l13_shuffle = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2 };

   /* 5 shuffle xxm's + 5 * 2 other xxms == 15 xxms so this is maxes out */

   while (c--) {
      *out++ = lp7[0] | lp11[0] | lp13[0];
      lp7[0]  = __builtin_shuffle(lp7[0],  l7_shuffle);
      lp11[0] = __builtin_shuffle(lp11[0], l11_shuffle);
      lp13[0] = __builtin_shuffle(lp13[0], l13_shuffle);
   }
}
#endif


/*
 * NOTE:
 * Shuffle of YMM registers only performs a shuffle within the 2 128 bit parts
 * and not across "lanes".
 *
 * THIS IS NOT WHAT I WANTED  T-T
 *
 */

#define SHUFFLE_IND(PRIME) \
   (0+64%PRIME)%64, \
   (1+64%PRIME)%64, \
   (2+64%PRIME)%64, \
   (3+64%PRIME)%64, \
   (4+64%PRIME)%64, \
   (5+64%PRIME)%64, \
   (6+64%PRIME)%64, \
   (7+64%PRIME)%64, \
   (8+64%PRIME)%64, \
   (9+64%PRIME)%64, \
   (10+64%PRIME)%64, \
   (11+64%PRIME)%64, \
   (12+64%PRIME)%64, \
   (13+64%PRIME)%64, \
   (14+64%PRIME)%64, \
   (15+64%PRIME)%64, \
   (16+64%PRIME)%64, \
   (17+64%PRIME)%64, \
   (18+64%PRIME)%64, \
   (19+64%PRIME)%64, \
   (20+64%PRIME)%64, \
   (21+64%PRIME)%64, \
   (22+64%PRIME)%64, \
   (23+64%PRIME)%64, \
   (24+64%PRIME)%64, \
   (25+64%PRIME)%64, \
   (26+64%PRIME)%64, \
   (27+64%PRIME)%64, \
   (28+64%PRIME)%64, \
   (29+64%PRIME)%64, \
   (30+64%PRIME)%64, \
   (31+64%PRIME)%64, \


#if 0
static void
do_early_prime_group_2_all_again(struct shuffle2_ctx *sctx, struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8)/4;
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 512);

   DO_FOR(DEC_X_V32_B, ALL_3, ALL_4)

   memcpy(lp_17, sctx->lp_17, sizeof lp_17);
   memcpy(lp_19, sctx->lp_19, sizeof lp_19);
   memcpy(lp_23, sctx->lp_23, sizeof lp_23);
   memcpy(lp_29, sctx->lp_29, sizeof lp_29);
   memcpy(lp_31, sctx->lp_31, sizeof lp_31);

   c17 = (pcb->block_start_num / 30) % 17;
   c19 = (pcb->block_start_num / 30) % 19;
   c23 = (pcb->block_start_num / 30) % 23;
   c29 = (pcb->block_start_num / 30) % 29;
   c31 = (pcb->block_start_num / 30) % 31;

   while (c--) {
      l17 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_17+c17));
      l19 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_19+c19));
      l23 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_23+c23));
      l29 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_29+c29));
      l31 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_31+c31));
      *out++ |= l17 | l19 | l23 | l29 | l31;
      c17 += (32 - 17);
      c17 -= c17>=17?17:0;
      c19 += (32 - 19);
      c19 -= c19>=19?19:0;
      c23 += (32 - 23);
      c23 -= c23>=23?23:0;
      c29 += (32 - 29);
      c29 -= c29>=29?29:0;
      c31 += (32 - 31);
      c31 -= c31>=31?31:0;
   }
}
#endif


static void
do_early_prime_group_2_all(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8)/8;
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 512);

   DO_FOR(DEC_X_V32_F, ALL_3, ALL_4)
   DO_FOR(INIT_X_32_F, ALL_3, ALL_4)

   v32ui l17_shuffle = { SHUFFLE_IND(17) };
   v32ui l19_shuffle = { SHUFFLE_IND(19) };
   v32ui l23_shuffle = { SHUFFLE_IND(23) };
   v32ui l29_shuffle = { SHUFFLE_IND(29) };
   v32ui l31_shuffle = { SHUFFLE_IND(31) };

   while (c--) {
      *out++ |= lp17[0] | lp19[0] | lp23[0] | lp29[0] | lp31[0];
      *out++ |= lp17[1] | lp19[1] | lp23[1] | lp29[1] | lp31[1];

      lp17[0] = __builtin_shuffle(lp17[0], lp17[1], l17_shuffle);
      lp17[1] = __builtin_shuffle(lp17[1], lp17[0], l17_shuffle);
      lp19[0] = __builtin_shuffle(lp19[0], lp19[1], l19_shuffle);
      lp19[1] = __builtin_shuffle(lp19[1], lp19[0], l19_shuffle);
      lp23[0] = __builtin_shuffle(lp23[0], lp23[1], l23_shuffle);
      lp23[1] = __builtin_shuffle(lp23[1], lp23[0], l23_shuffle);
      lp29[0] = __builtin_shuffle(lp29[0], lp29[1], l29_shuffle);
      lp29[1] = __builtin_shuffle(lp29[1], lp29[0], l29_shuffle);
      lp31[0] = __builtin_shuffle(lp31[0], lp31[1], l31_shuffle);
      lp31[1] = __builtin_shuffle(lp31[1], lp31[0], l31_shuffle);
   }
}


static void
do_early_prime_group_6(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/ 8 + 7)/8;
   int i;
   uint64_t *bm = (uint64_t *)pcb->block;

   DO_FOR(DEC_X, ALL_8)
   DO_FOR(INIT_X, ALL_8)

   while (c--) {
      for (i = 0; i < 8; i++)
         *bm++ |= 0ul DO_FOR(OR_X, ALL_8);
      DO_FOR(SHIFT_X, ALL_8)
   }
}


static void
do_early_prime_group_2_c(struct shuffle2_ctx *sctx, struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8)/8;
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 512);

   DO_FOR(DEC_X_V32_C, ALL_5, ALL_6)

   memcpy(lp_37, sctx->lp_37, sizeof lp_37);
   memcpy(lp_41, sctx->lp_41, sizeof lp_41);
   memcpy(lp_43, sctx->lp_43, sizeof lp_43);
   memcpy(lp_47, sctx->lp_47, sizeof lp_47);

   c37 = (pcb->block_start_num / 30) % 37;
   c41 = (pcb->block_start_num / 30) % 41;
   c43 = (pcb->block_start_num / 30) % 43;
   c47 = (pcb->block_start_num / 30) % 47;

   while (c--) {
      l37_1 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_37+c37));
      l37_2 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_37+c37+32));
      l41_1 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_41+c41));
      l41_2 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_41+c41+32));
      l43_1 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_43+c43));
      l43_2 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_43+c43+32));
      l47_1 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_47+c47));
      l47_2 = (v32ui)_mm256_loadu_si256((__m256i *)((char *)lp_47+c47+32));

      *out++ |= l37_1 | l41_1 | l43_1 | l47_1;
      *out++ |= l37_2 | l41_2 | l43_2 | l47_2;

      c37 += (64 - 37);
      c37 -= c37>=37?37:0;
      c41 += (64 - 41);
      c41 -= c41>=41?41:0;
      c43 += (64 - 43);
      c43 -= c43>=43?43:0;
      c47 += (64 - 47);
      c47 -= c47>=47?47:0;
   }
}


static void
do_early_prime_group_7(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8 + 6)/7;
   int i;
   uint64_t *bm = (uint64_t *)pcb->block;

   DO_FOR(DEC_X, ALL_7)
   DO_FOR(INIT_X, ALL_7)

   while (c--) {
      for (i = 0; i < 7; i++)
         *bm++ |= 0ul DO_FOR(OR_X, ALL_7);
      DO_FOR(SHIFT_X, ALL_7)
   }
}


int
shuffle2_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct shuffle2_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   if (sctx->start_prime <= 7 && sctx->end_prime >= 7)
      do_early_prime_7(sctx, pcb);
   if (sctx->start_prime <= 11 && sctx->end_prime >= 11)
      do_early_prime_11(sctx, pcb);
   if (sctx->start_prime <= 13 && sctx->end_prime >= 13)
      do_early_prime_13(sctx, pcb);
   if (sctx->start_prime <= 17 && sctx->end_prime >= 31)
      do_early_prime_group_2_all(pcb);
   if (sctx->start_prime <= 37 && sctx->end_prime >= 47)
      do_early_prime_group_2_c(sctx, pcb);
   if (sctx->start_prime <= 53 && sctx->end_prime >= 53)
      do_early_prime_group_7(pcb);
   if (sctx->start_prime <= 59 && sctx->end_prime >= 61)
      do_early_prime_group_6(pcb);

   return 0;
}
