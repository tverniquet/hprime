#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <assert.h>
#include <nmmintrin.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct shuffle_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
};


int
shuffle_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct shuffle_ctx *sctx = malloc(sizeof (struct shuffle_ctx));
   *ctx = sctx;

   assert(end_prime < 64);
   assert(start_prime < end_prime);

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   return 0;
}


int
shuffle_free(void *ctx)
{
   struct shuffle_ctx *sctx = ctx;
   bzero(sctx, sizeof *sctx);
   FREE(ctx);
   return 0;
}


int
shuffle_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct shuffle_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;
   }
   return 0;
}


int
shuffle_skip_to(struct prime_thread_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
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
         if (cur > A + num*8*30) {
            /*for (i = 0; i < (int)num; i++) {*/
               /*fprintf(stdout, "%d init %d is %016lX\n", a, i, buf[i]);*/
            /*}*/
            return;
         }
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

#define DEC_X_V16(PRIME) uint64_t lp_##PRIME [CEIL_i64_A(PRIME)] __attribute__((aligned(16))); v16ui *lp##PRIME = (v16ui *)lp_##PRIME;
#define DEC_X(PRIME) uint64_t lp_##PRIME [CEIL_i64_A(PRIME)];

#define INIT_X(PRIME) init_early_prime_a(pcb->block_start_num, PRIME, lp_##PRIME, CEIL_i64_A(PRIME));
#define OR_X(PRIME)     |(lp_##PRIME[i % (PRIME + 7 / 8)])
#define SHIFT_X(PRIME) shift_a(PRIME, lp_##PRIME);


typedef unsigned char v16ui __attribute__((vector_size(16)));
typedef uint16_t       v8ui __attribute__((vector_size(16)));
typedef uint64_t       v2ui __attribute__((vector_size(16)));
typedef uint64_t       v4ui __attribute__((vector_size(32)));



#if 0
static void print_v(v16ui *a)
{
   unsigned char t2[16] __attribute__((aligned(16)));
   *(v16ui *)t2 = *a;
   fprintf(stderr, "%02X %02X %02X %02X  %02X %02X %02X %02X  %02X %02X %02X %02X  %02X %02X %02X %02X  \n", t2[0], t2[1], t2[2], t2[3], t2[4], t2[5], t2[6], t2[7], t2[8], t2[9], t2[10], t2[11], t2[12], t2[13], t2[14], t2[15]);
}
#endif

/* I use clang as the syntax checker but compile with gcc. But clang doesn't seem to know about shuffle */
#ifdef __clang__
#define DO_NONE(A) (void)A,
#define __builtin_shuffle(A, ...) (DO_FOR(DO_NONE, __VA_ARGS__) A)
#define __builtin_assume_aligned(A,B) A
#endif


/* TODO: It is silly 'init'ing the primes every time when they just
 * continue on..?
 *
 * The idea here is to avoid the shifting work by using more registers to
 */
static inline void __attribute__((always_inline))
do_early_prime_group_1_type (struct prime_current_block *pcb, const int a, const int mode)
{
   uint64_t l_p[2*a];
   int c = (pcb->block_size/8 + (2*a - 1))/(2*a);
   int i;
   v16ui *lp = (v16ui *)__builtin_assume_aligned(l_p, 16);
   v16ui *out =  (v16ui *)__builtin_assume_aligned(pcb->block,  512);

   init_early_prime_a(pcb->block_start_num, a, l_p, 2*a);
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
static void do_early_prime_7(struct prime_current_block *pcb) { do_early_prime_group_1_type(pcb, 7, 1); }
static void do_early_prime_11(struct prime_current_block *pcb) { do_early_prime_group_1_type(pcb, 11, 0); }
static void do_early_prime_13(struct prime_current_block *pcb) { do_early_prime_group_1_type(pcb, 13, 0); }
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

static void
do_early_prime_group_2(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8)/4;
   v16ui *out = (v16ui *)pcb->block;

   DO_FOR(DEC_X_V16, ALL_3, ALL_4)
   DO_FOR(INIT_X, ALL_3, ALL_4)

   v16ui l17_shuffle = { 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30 };
   v16ui l19_shuffle = { 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28 };
   v16ui l23_shuffle = {  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24 };
   v16ui l29_shuffle = {  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18 };
   v16ui l31_shuffle = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16 };

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

/* This one is a bit odd - maybe make 5,7 also '8' */
static void
do_early_prime_group_4(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8 + 4)/5;
   int i;
   uint64_t *bm = (uint64_t *)pcb->block;
   DO_FOR(DEC_X, ALL_5)
   DO_FOR(INIT_X, ALL_5)

   while (c--) {
      for (i = 0; i < 5; i++)
         *bm++ |= 0ul DO_FOR(OR_X, ALL_5);
      DO_FOR(SHIFT_X, ALL_5)
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

static void
do_early_prime_group_5(struct prime_current_block *pcb)
{
   int c = (pcb->block_size/8 + 5)/6;
   int i;
   uint64_t *bm = (uint64_t *)pcb->block;

   DO_FOR(DEC_X, ALL_6)
   DO_FOR(INIT_X, ALL_6)

   while (c--) {
      for (i = 0; i < 6; i++)
         *bm++ |= 0ul DO_FOR(OR_X, ALL_6);
      DO_FOR(SHIFT_X, ALL_6)
   }
}


int
shuffle_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct shuffle_ctx *sctx = ctx;
   /* Mark off multiples of 'a' in the block */
   if (sctx->start_prime <= 7 && sctx->end_prime >= 7)
      do_early_prime_7(pcb);
   if (sctx->start_prime <= 11 && sctx->end_prime >= 11)
      do_early_prime_11(pcb);
   if (sctx->start_prime <= 13 && sctx->end_prime >= 13)
      do_early_prime_13(pcb);
   if (sctx->start_prime <= 17 && sctx->end_prime >= 31)
      do_early_prime_group_2(pcb);
   if (sctx->start_prime <= 37 && sctx->end_prime >= 37)
      do_early_prime_group_4(pcb);
   if (sctx->start_prime <= 41 && sctx->end_prime >= 47)
      do_early_prime_group_5(pcb);
   if (sctx->start_prime <= 59 && sctx->end_prime >= 61)
      do_early_prime_group_6(pcb);
   if (sctx->start_prime <= 53 && sctx->end_prime >= 53)
      do_early_prime_group_7(pcb);

   return 0;
}
