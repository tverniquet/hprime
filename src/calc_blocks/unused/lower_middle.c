#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <immintrin.h>

#include <stdio.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"


#define PRIMES_64_to_128 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127
#define USED_PRIMES 67, 71, 73, 79, 83, 89, 97


/*
 *  Target 64 - 128
 *
 * Using the abit trick:
 *
 *  These primes set one bit ~ every 2 cache lines
 *
 *  can fit more than 128 bytes into ymm registers...
 *
 *  So load 4 x output into ymm registers
 *
 *  Iterate through the 64-128 primes setting the bit for the next two blocks
 *
 *  (will need to check if they need to set a second bit, eg if 67 sets a bit
 *  at the start of one 64 block it will need to set another one within 128 bytes)
 *
 *  So the trick will be making the iterating through the primes fast.
 *
 *  Have 13 primes between 64-128
 *
 *  67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127
 *
 * For each of 8 bits need 7 bits to store offset of bit - so 64bit int can store
 * one block.
 *
 * So each prime can sit in a register..
 *
 *
 * XXXXXXX
 *
 *
 * Instead:
 *
 * for each of the 8 output bits:
 *
 *   each bit is now a byte
 *   so for eg, set every 67th bit
 *   so use the unaligned load trick?? 67+32 bytes...
 *   Put each prime map in a register
 *
 *   while have more output
 *      for each prime
 *         or with the output mask
 *         shift
 *      for each 4x8 parts of the output mask
 *         explode into a ymm according to the output bit
 *         output the result
 */


struct lower_middle_ctx
{
   uint32_t start_prime;
   uint32_t end_prime;
   uint32_t block_size;
   unsigned char *primebuf[13];

   uint64_t offsets[13*8*8] __attribute__((aligned(32)));
   uint64_t offsets_v2[13*4*4] __attribute__((aligned(32)));
   uint32_t offsets_v3[8];
};


/* Will need prime * 8 bits (so prime bytes) to capture full pattern... */


/*
 * hmm,
 *
 * given 128 bits, setting 2 x 67th bits:
 *
 * x...............................x......
 *
 * There is at most 2 bits set in a 128bit sequence
 * So just shift the starting register by the set amount and it should be fine..
 * if offset += (128-67)
 * offset -= offset >= 67 ? 67 : 0;
 *
 * next = orig >> offset;
 *
 */

#if 0
static void
set_every_nth_bit (unsigned char *buf, const uint32_t a, int nbytes)
{
   int i;
   for (i = 0; i < nbytes*8; i+=a) {
      buf[i/8] |= 1 << (i%8);
   }
}
#endif


static void
set_starting_byte (const uint32_t A, const uint32_t a, uint64_t *buf)
{
   int nums[] = {1,7,11,13,17,19,23,29};
   uint64_t s;
   int starting_bytes[8];
   int i;
   int j;
   uint64_t b_num = A/a/30*30;

   for (i = 0; i < 8; i++) {
      s = (b_num + nums[i]) * a;
      if (s < A)
         s += (a*30);
      starting_bytes[((s%30)*8/30)] = (s / 30) - (A / 30);
   }

   for (i = 0; i < 8; i++)
      for (j = 0; j < 8; j++)
         buf[i*8 + j] = ((j*a)-(j*64) + starting_bytes[i])%a;
}


static void
set_starting_v2 (const uint32_t A, const uint32_t a, uint64_t *buf)
{
   int nums[] = {1,7,11,13,17,19,23,29};
   uint64_t s;
   int starting_bytes[8];
   int i;
   int j;
   uint64_t b_num = A/a/30*30;

   for (i = 0; i < 8; i++) {
      s = (b_num + nums[i]) * a;
      if (s < A)
         s += (a*30);
      starting_bytes[((s%30)*8/30)] = (s / 30) - (A / 30);
   }

   for (i = 0; i < 4; i++) {
      for (j = 0; j < 2; j++)
         buf[i*4 + j] = ((j*a)-(j*64) + starting_bytes[2*i])%a;
      for (j = 0; j < 2; j++)
         buf[i*4 + 2 + j] = ((j*a)-(j*64) + starting_bytes[2*i+1])%a;
   }
}


static void
init_offsets(uint64_t A, uint32_t a, uint32_t *offsets)
{
   int bit;
   uint64_t multiplier_base = bytes_to_num(num_to_bytes(MAX(A / a, a)));
   uint64_t multiplier_byte;

   for (bit = 0; bit < 8; bit++) {
      multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * a);

      while (multiplier_byte < num_to_bytes(A))
         multiplier_byte += a;

      while (multiplier_byte - num_to_bytes(A) >= UINT16_MAX)
         multiplier_byte -= a;

      offsets[bit] = multiplier_byte - num_to_bytes(A);
   }
}




int
lower_middle_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct lower_middle_ctx *sctx = malloc(sizeof (struct lower_middle_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   assert(sctx->start_prime >= 64);
   assert(sctx->end_prime <= 128);

   sctx->block_size = pctx->current_block.block_size;

   sctx->primebuf[0] = aligned_alloc(32, 67+32);
   sctx->primebuf[1] = aligned_alloc(32, 67+32);

   int k = 0;

   init_offsets(0, 67, sctx->offsets_v3);

#define LM_INIT_X(PRIME) \
   set_starting_v2(0, PRIME, &sctx->offsets_v2[k++*32]);

   DO_FOR(LM_INIT_X, USED_PRIMES)

   set_starting_byte(0, 67, &sctx->offsets[k++*64]);
   set_starting_byte(0, 71, &sctx->offsets[k++*64]);
   set_starting_byte(0, 73, &sctx->offsets[k++*64]);
   set_starting_byte(0, 79, &sctx->offsets[k++*64]);
   set_starting_byte(0, 83, &sctx->offsets[k++*64]);
   set_starting_byte(0, 89, &sctx->offsets[k++*64]);

   return 0;
}


int
lower_middle_free(void *ctx)
{
   FREE(ctx);
   return 0;
}


int
lower_middle_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct lower_middle_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;
   }
   return 0;
}


int
lower_middle_skip_to(struct prime_ctx *__attribute__((unused))pctx, uint64_t __attribute__((unused))target_num, void *__attribute__((unused))ctx)
{
   return 0;
}



/*
 *
 *  0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32   33
 *|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|--|  |--|
 * |                       |                       |                          |                          |
 * 0                       3                       6                          1                          4
 *
 *|----------------------|------------------------|-----------------------|-----------------------|
 *
 * - 0                      -64                        -128                     -192
 * +67                       67                        67                       67
 *  67                      [3]                        -61                      -125
 * +67                       67                        67                       67
 *  134                      70                        [6]                      -58
 * +67                       67                        67                       67
 *  201                     137                        73                       [9]
 * +67                       67                        67                       67
 * % 256
 * [12]                     204                       140                       76
 * 
 * 
 *
 *
 *
 * xxx
 *
 * if A == 32*1024*30 == 983040
 *
 * offset[8];
 *
 * b_base = A / 67 / 30 * 30 == 14670
 * for each b_bit = 0 .. 7
 *    res = (b_base + ind_to_mod[b_bit]) * 67;
 *    diff = 983040 - res;
 *    if (diff < 0)
 *    diff += 67;
 *    offset[res_bit] = diff;
 *
 * 
 * for (i = 0; i < 8; i++) {
 *   bit_to_offset[67*i%8] = 67*i/8;
 * }
 *
 * NOW each bit represents a byte. The first bit is:
 *   offset[res_bit] % 8;
 *
 * the location of the byte that has the bit set:
 *   bit_to_offset[bit];
 *
 * The number of bytes beforehand is:
 *   diff / 8;
 *
 * so the load_offset is
 *
 */
#if 0
static void
compute_67_71_73_79(struct lower_middle_ctx *sctx, struct prime_current_block *pcb)
{
   v32ui *out;
   v32_4ui x4 = {0xFF,0xFF,0xFF,0xFF};
   v32_4ui x7 = {1,1,1,1};
   v32_4ui x8;
   v32ui b0 = (v32ui) { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3 };
   v32ui b1 = (v32ui) { 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 };

   v32ui b;
   v32ui c;

   v32ui d =  {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   int i;

   /* 67 */
   v32_4ui x67_3 = {7*67,7*67,7*67,7*67};
   v32_4ui x67_5 = {67,67,67,67};
   v32_4ui x67_6 = {43,43,43,43};

   v32_4ui x67_1;
   v32_4ui x67_2;

   /* 71 */
   v32_4ui x71_3 = {7*71,7*71,7*71,7*71};
   v32_4ui x71_5 = {71,71,71,71};
   v32_4ui x71_6 = {15,15,15,15};

   v32_4ui x71_1;
   v32_4ui x71_2;

   /* 73 */
   v32_4ui x73_3 = {7*73,7*73,7*73,7*73};
   v32_4ui x73_5 = {73,73,73,73};
   v32_4ui x73_6 = {1,1,1,1};

   v32_4ui x73_1;
   v32_4ui x73_2;

   /* 79 */
   v32_4ui x79_3 = {6*79,6*79,6*79,6*79};
   v32_4ui x79_5 = {79,79,79,79};
   v32_4ui x79_6 = {38,38,38,38};

   v32_4ui x79_1;
   v32_4ui x79_2;

   for (i = 0; i < 8; i++) {
      x67_1 = *(v32_4ui *)&sctx->offsets[i*8];
      x67_2 = *(v32_4ui *)&sctx->offsets[i*8+4];
      x71_1 = *(v32_4ui *)&sctx->offsets[64 + i*8];
      x71_2 = *(v32_4ui *)&sctx->offsets[64 + i*8+4];
      x73_1 = *(v32_4ui *)&sctx->offsets[128 + i*8];
      x73_2 = *(v32_4ui *)&sctx->offsets[128 + i*8+4];
      x79_1 = *(v32_4ui *)&sctx->offsets[192 + i*8];
      x79_2 = *(v32_4ui *)&sctx->offsets[192 + i*8+4];

      v32ui g = {1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i};


      out = (v32ui *)__builtin_assume_aligned(pcb->block, 32);

      /* There are 2x32=64=512 bytes per run */
      int j = 64;
      while (j--) {
         x8 = (x7 << x67_1) | (x7 << x71_1) | (x7 << x73_1) | (x7 << x79_1);

         b = b0;
         int k;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }
         x8 = (x7 << x67_2) | (x7 << x71_2) | (x7 << x73_2) | (x7 << x79_2);
         b = b0;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }

         x67_1 = (x67_1 + x67_3 + (x67_5 & (x67_1 < x67_6))) & x4;
         x67_2 = (x67_2 + x67_3 + (x67_5 & (x67_2 < x67_6))) & x4;
         x71_1 = (x71_1 + x71_3 + (x71_5 & (x71_1 < x71_6))) & x4;
         x71_2 = (x71_2 + x71_3 + (x71_5 & (x71_2 < x71_6))) & x4;
         x73_1 = (x73_1 + x73_3 + (x73_5 & (x73_1 < x73_6))) & x4;
         x73_2 = (x73_2 + x73_3 + (x73_5 & (x73_2 < x73_6))) & x4;
         x79_1 = (x79_1 + x79_3 + (x79_5 & (x79_1 < x79_6))) & x4;
         x79_2 = (x79_2 + x79_3 + (x79_5 & (x79_2 < x79_6))) & x4;
      }

      *(v32_4ui *)&sctx->offsets[i*8] = x67_1;
      *(v32_4ui *)&sctx->offsets[i*8+4] = x67_2;
      *(v32_4ui *)&sctx->offsets[64 + i*8] = x71_1;
      *(v32_4ui *)&sctx->offsets[64 + i*8+4] = x71_2;
      *(v32_4ui *)&sctx->offsets[128 + i*8] = x73_1;
      *(v32_4ui *)&sctx->offsets[128 + i*8+4] = x73_2;
      *(v32_4ui *)&sctx->offsets[192 + i*8] = x79_1;
      *(v32_4ui *)&sctx->offsets[192 + i*8+4] = x79_2;

   }
}
static void
compute_67_71_73(struct lower_middle_ctx *sctx, struct prime_current_block *pcb)
{
   v32ui *out;
   v32_4ui x4 = {0xFF,0xFF,0xFF,0xFF};
   v32_4ui x7 = {1,1,1,1};
   v32_4ui x8;
   v32ui b0 = (v32ui) { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3 };
   v32ui b1 = (v32ui) { 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 };

   v32ui b;
   v32ui c;

   v32ui d =  {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   int i;

   /* 67 */
   v32_4ui x67_3 = {7*67,7*67,7*67,7*67};
   v32_4ui x67_5 = {67,67,67,67};
   v32_4ui x67_6 = {43,43,43,43};

   v32_4ui x67_1;
   v32_4ui x67_2;

   /* 71 */
   v32_4ui x71_3 = {7*71,7*71,7*71,7*71};
   v32_4ui x71_5 = {71,71,71,71};
   v32_4ui x71_6 = {15,15,15,15};

   v32_4ui x71_1;
   v32_4ui x71_2;

   /* 73 */
   v32_4ui x73_3 = {7*73,7*73,7*73,7*73};
   v32_4ui x73_5 = {73,73,73,73};
   v32_4ui x73_6 = {1,1,1,1};

   v32_4ui x73_1;
   v32_4ui x73_2;

   for (i = 0; i < 8; i++) {
      x67_1 = *(v32_4ui *)&sctx->offsets[i*8];
      x67_2 = *(v32_4ui *)&sctx->offsets[i*8+4];
      x71_1 = *(v32_4ui *)&sctx->offsets[64 + i*8];
      x71_2 = *(v32_4ui *)&sctx->offsets[64 + i*8+4];
      x73_1 = *(v32_4ui *)&sctx->offsets[128 + i*8];
      x73_2 = *(v32_4ui *)&sctx->offsets[128 + i*8+4];

      v32ui g = {1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i};


      out = (v32ui *)__builtin_assume_aligned(pcb->block, 32);

      /* There are 2x32=64=512 bytes per run */
      int j = 64;
      while (j--) {
         x8 = (x7 << x67_1) | (x7 << x71_1) | (x7 << x73_1);

         b = b0;
         int k;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }
         x8 = (x7 << x67_2) | (x7 << x71_2) | (x7 << x73_2);
         b = b0;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }

         x67_1 = (x67_1 + x67_3 + (x67_5 & (x67_1 < x67_6))) & x4;
         x67_2 = (x67_2 + x67_3 + (x67_5 & (x67_2 < x67_6))) & x4;
         x71_1 = (x71_1 + x71_3 + (x71_5 & (x71_1 < x71_6))) & x4;
         x71_2 = (x71_2 + x71_3 + (x71_5 & (x71_2 < x71_6))) & x4;
         x73_1 = (x73_1 + x73_3 + (x73_5 & (x73_1 < x73_6))) & x4;
         x73_2 = (x73_2 + x73_3 + (x73_5 & (x73_2 < x73_6))) & x4;
      }

      *(v32_4ui *)&sctx->offsets[i*8] = x67_1;
      *(v32_4ui *)&sctx->offsets[i*8+4] = x67_2;
      *(v32_4ui *)&sctx->offsets[64 + i*8] = x71_1;
      *(v32_4ui *)&sctx->offsets[64 + i*8+4] = x71_2;
      *(v32_4ui *)&sctx->offsets[128 + i*8] = x73_1;
      *(v32_4ui *)&sctx->offsets[128 + i*8+4] = x73_2;

   }
}

static void
compute_67_71(struct lower_middle_ctx *sctx, struct prime_current_block *pcb)
{
   v32ui *out;
   v32_4ui x4 = {0xFF,0xFF,0xFF,0xFF};
   v32_4ui x7 = {1,1,1,1};
   v32_4ui x8;
   v32ui b0 = (v32ui) { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3 };
   v32ui b1 = (v32ui) { 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 };

   v32ui b;
   v32ui c;

   v32ui d =  {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   int i;

   /* 67 */
   v32_4ui x67_3 = {7*67,7*67,7*67,7*67};
   v32_4ui x67_5 = {67,67,67,67};
   v32_4ui x67_6 = {43,43,43,43};

   v32_4ui x67_1;
   v32_4ui x67_2;

   /* 71 */
   v32_4ui x71_3 = {7*71,7*71,7*71,7*71};
   v32_4ui x71_5 = {71,71,71,71};
   v32_4ui x71_6 = {15,15,15,15};

   v32_4ui x71_1;
   v32_4ui x71_2;


   for (i = 0; i < 8; i++) {
      x67_1 = *(v32_4ui *)&sctx->offsets[i*8];
      x67_2 = *(v32_4ui *)&sctx->offsets[i*8+4];
      x71_1 = *(v32_4ui *)&sctx->offsets[64 + i*8];
      x71_2 = *(v32_4ui *)&sctx->offsets[64 + i*8+4];

      v32ui g = {1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i};


      out = (v32ui *)__builtin_assume_aligned(pcb->block, 32);

      /* There are 2x32=64=512 bytes per run */
      int j = 64;
      while (j--) {
         x8 = (x7 << x67_1) | (x7 << x71_1);

         b = b0;
         int k;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }
         x8 = (x7 << x67_2) | (x7 << x71_2);
         b = b0;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }

         x67_1 = (x67_1 + x67_3 + (x67_5 & (x67_1 < x67_6))) & x4;
         x67_2 = (x67_2 + x67_3 + (x67_5 & (x67_2 < x67_6))) & x4;
         x71_1 = (x71_1 + x71_3 + (x71_5 & (x71_1 < x71_6))) & x4;
         x71_2 = (x71_2 + x71_3 + (x71_5 & (x71_2 < x71_6))) & x4;
      }

      *(v32_4ui *)&sctx->offsets[i*8] = x67_1;
      *(v32_4ui *)&sctx->offsets[i*8+4] = x67_2;
      *(v32_4ui *)&sctx->offsets[64 + i*8] = x71_1;
      *(v32_4ui *)&sctx->offsets[64 + i*8+4] = x71_2;

   }
}

static void
compute_all(struct lower_middle_ctx *sctx, struct prime_current_block *pcb)
{
   v32ui *out;
   v32_4ui x1;
   v32_4ui x2;
   v32_4ui x3 = {7*67,7*67,7*67,7*67};
   v32_4ui x4 = {0xFF,0xFF,0xFF,0xFF};
   v32_4ui x5 = {67,67,67,67};
   v32_4ui x6 = {43,43,43,43};
   v32_4ui x7 = {1,1,1,1};
   v32_4ui x8;
   v32ui b0 = (v32ui) { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3 };
   v32ui b1 = (v32ui) { 4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4 };

   v32ui b;
   v32ui c;

   v32ui d =  {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   int i;

   for (i = 0; i < 8; i++) {
      x1 = *(v32_4ui *)&sctx->offsets[i*8];
      x2 = *(v32_4ui *)&sctx->offsets[i*8+4];

      v32ui g = {1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,
                 1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i,1<<i};


      out = (v32ui *)__builtin_assume_aligned(pcb->block, 32);

      /* There are 2x32=64=512 bytes per run */
      int j = 64;
      while (j--) {
         x8 = x7 << x1;

         b = b0;
         int k;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }
         x8 = x7 << x2;
         b = b0;
         for (k = 0; k < 8; k++) {
            /* note 2 operations to do the shuffle */
            c = __builtin_shuffle ((v32ui)x8, b);
            *out++ |= (((c&d) == d) & g);
            b += b1;
         }

         x1 = (x1 + x3 + (x5 & (x1 < x6))) & x4;
         x2 = (x2 + x3 + (x5 & (x2 < x6))) & x4;
      }

      *(v32_4ui *)&sctx->offsets[i*8] = x1;
      *(v32_4ui *)&sctx->offsets[i*8+4] = x2;
   }
}


static void
compute_all_v2(struct lower_middle_ctx *sctx, struct prime_current_block *pcb)
{
   v32ui *out;
   int i;
   int m;
   v32ui b0 = (v32ui){ 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1 };
   v32ui b1 = (v32ui) { 2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2 };
   v32ui d =  {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};
   v32_4ui x4 = {0x7F,0x7F,0x7F,0x7F};
   v32_4ui x7 = {1,1,1,1};

#define CAv2_DECX(PRIME) \
   v32_4ui x##PRIME##_3 = {PRIME,PRIME,PRIME,PRIME}; \
   v32_4ui x##PRIME##_6 = {128%PRIME,128%PRIME,128%PRIME,128%PRIME};

   DO_FOR(CAv2_DECX, USED_PRIMES)


   for (i = 0; i < 8; i+=2) {
      m = 0;
#define CAv2_GET_OFFS(PRIME) \
      v32_4ui x##PRIME##_1 = *(v32_4ui *)&sctx->offsets_v2[m++*32 + i/2*4];

      DO_FOR(CAv2_GET_OFFS, USED_PRIMES);

      v32_4ui x8;


      v32_4ui g = { 0x0101010101010101ul, 0x0101010101010101ul, 0x0202020202020202ul, 0x0202020202020202ul };
      g <<= (v32_4ui){i,i,i,i};

      v32ui b;
      v32ui c;
      v32ui c1;
      v32ui m1;
      v32ui m2;

      out = (v32ui *)__builtin_assume_aligned(pcb->block, 32);

      /* There are xxx bytes per run */
      int j = 256;
      while (j--) {
         x8 = x7 << x67_1;
#define CAv2_OR_X(PRIME) x8 |= x7 << x##PRIME##_1;

         DO_FOR(CAv2_OR_X, USED_PRIMES)

         b = b0;
         int k;
         for (k = 0; k < 4; k++) {

            c = (v32ui)_mm256_shuffle_epi8((__m256i)x8, (__m256i)b);
            c = (((c&d) == d) & (v32ui)g);
            b += b1;

            c1 = (v32ui)_mm256_shuffle_epi8((__m256i)x8, (__m256i)b);
            c1 = (((c1&d) == d) & (v32ui)g);
            b += b1;

            m1 = (v32ui)_mm256_permute2f128_si256((__m256i)c, (__m256i)c1, 0x20);
            m2 = (v32ui)_mm256_permute2f128_si256((__m256i)c, (__m256i)c1, 0x31);
            *out++ |= m1 | m2;
         }

#define CAv2_INC_X(PRIME) \
         x##PRIME##_1 = (x##PRIME##_1 + x##PRIME##_3 + (x##PRIME##_3 & (x##PRIME##_1 < x##PRIME##_6))) & x4;

         DO_FOR(CAv2_INC_X, USED_PRIMES)
      }
      m = 0;
#define CAv2_PUT_OFFS(PRIME) \
      *(v32_4ui *)&sctx->offsets_v2[m++ * 32 + i/2*4] = x##PRIME##_1;

      DO_FOR(CAv2_PUT_OFFS, USED_PRIMES);
   }
}
#endif

#if 0
   v32ui *out = (v32ui *)__builtin_assume_aligned(pcb->block, 32);

   int c = CEIL_DIV(pcb->block_size, 4*32);
   int i;
   int offset[16];

   v32_4ui x1 = {(0*67)-(0*64),(1*67)-(1*64),(2*67)-(2*64),(3*67)-(3*64)};
   v32_4ui x3 = {7*67,7*67,7*67,7*67};
   v32_4ui x4 = {0xFF,0xFF,0xFF,0xFF};
   v32_4ui x5 = {67,67,67,67};
   v32_4ui x6 = {256 - (256/67), 43,43,43,43};
   v32_4ui x7 = {1,1,1,1};
   v32_4ui x8;

   uint32_t bitm;


   int i;
   for (i = 0; i < 100; i++) {
      x8 = x7 << x1;
      print_v32((v32ui*)&x8);
      x1 = (x1 + x3 + (x5 & (x1 < x6))) & x4;
   }



   v32ui a = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
         b, g, y;

   v32ui d =  {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80,
               0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};





   i = 0;
#define CA_SET_INITIAL_OFFSET(PRIME) offset[i++] = pcb->block_start_num / 30 % PRIME;
   DO_FOR( CA_SET_INITIAL_OFFSET, PRIMES_64_to_128)

   while (c--) {
      b = (v32ui) { 0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3 };
      g = (v32ui) { 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1 };

      y = __builtin_shuffle(a, b);


      i = 0;

   }
}
#endif

/*
 * The first time we don't know how many times to mark off the prime in the
 * block (ie don't really know 'n'). Anyway, keeping this separate to the
 * calc_block function allows the "calc_block" function to make more
 * assumptions and optimise better.
 */
static void
compute_block_first_time(struct prime_current_block *pcb, uint32_t sieve_prime, uint32_t *offsets)
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
compute_block(struct prime_current_block *pcb, uint32_t sieve_prime, uint32_t *offsets, uint32_t n)
{
   char *bms[8];
   int i;
   int n2 = n;
   const unsigned char *bits = a_x_b_bitmask[pp_to_bit(sieve_prime)];

   for (i = 0; i < 8; i++)
      bms[i] = (char *)pcb->block + offsets[i];

   while (n--) {
      for (i = 0; i < 8; i++) {
         *bms[i]   |= bits[i];
         bms[i]   += sieve_prime;
      }
   }

   for (i = 0; i < 8; i++)
      bms[i] = (char *)pcb->block + offsets[i];

   n = n2;

   while (n--) {
      for (i = 0; i < 8; i++) {
         *bms[i]   |= bits[i];
         bms[i]   += sieve_prime;
      }
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


static void
check_new_sieve_primes(struct lower_middle_ctx *sctx, struct prime_current_block *pcb)
{
   int bit;
   uint32_t sieve_prime;
   uint64_t multiplier_base;
   uint64_t multiplier_byte;

   uint64_t block_start_byte = num_to_bytes(pcb->block_start_num);

   sieve_prime = 67;

   multiplier_base = bytes_to_num(num_to_bytes(MAX(pcb->block_start_num / sieve_prime, sieve_prime)));

   for (bit = 0; bit < 8; bit++) {
      multiplier_byte = num_to_bytes((multiplier_base + ind_to_mod[bit]) * sieve_prime);
      while (multiplier_byte < block_start_byte)
         multiplier_byte += sieve_prime;

      while (multiplier_byte - block_start_byte >= UINT16_MAX)
         multiplier_byte -= sieve_prime;

      sctx->offsets_v3[bit] = multiplier_byte - block_start_byte;
   }
   compute_block_first_time(pcb, 67, &sctx->offsets_v3[0]);
}


int
lower_middle_calc_primes(struct prime_current_block *pcb, void *ctx)
{
   struct lower_middle_ctx *sctx = ctx;
   if (pcb->block_start_num == 0)
      check_new_sieve_primes(sctx, pcb);
   else
      compute_block(pcb, 67, sctx->offsets_v3, sctx->block_size/67);

   /*compute_all_v2(sctx, pcb);*/

   if (pcb->block_start_num == 0) {
      pcb->block[2] = (1<<4);
   }

   return 0;
}
