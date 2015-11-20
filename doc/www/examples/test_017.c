#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))

const char diff[] = {6,4,2,4,2,4,6,2};
const char bval[] = {1,7,11,13,17,19,23,29};

static inline int num2bit(uint64_t num) { return (num % 30) * 8 / 30; }


static void
calc_primes(char *bitmap, uint64_t max)
{
   char     *bmp,
            *bitmap_end;
   uint64_t  a,
             cur;
   uint32_t *as   = malloc(sizeof(uint32_t) * (sqrtl(max) / 4 + 1000)),
            *bmps = malloc(sizeof(uint32_t) * (sqrtl(max) / 4 + 1000)),
            *offs = malloc(sizeof(uint32_t) * 8 * (sqrtl(max) / 4 + 1000));
   int       i,
             j,
             a_i,
             cnt = 0;

   /*** Main algorithm start ***/

   /*
    * Find sieving primes
    */
   bitmap_end = bitmap + (int64_t)sqrtl(max) + 1;
   for (a = 7, a_i = 1; a <= sqrtl(max); a += diff[a_i++%8]) {
      if (bitmap[a / 30] & 1 << (a_i%8))
         continue;

      as[cnt] = a;
      for (i = 0; i < 8; i++)
         offs[cnt*8 + num2bit(a * bval[i])] = a * bval[i] / 30;

      for (bmp = bitmap + (uint64_t)a * (a / 30); bmp < bitmap_end; bmp += a)
         for (i = 0; i < 8; i++)
            *(bmp + offs[cnt*8 + i]) |= 1 << i;

      bmps[cnt++] = (bmp - bitmap);
   }

   /*
    * Calculate blocks
    */
   bitmap_end = bitmap + 32*1024;
   for (cur = 0; cur <= max; cur += 32*1024*30, bitmap_end += 32*1024) {

      for (i = 0; i < cnt; i++)
         for (; bmps[i] < cur/30 + 32*1024; bmps[i] += as[i])
            for (j = 0; j < 8; j++)
               *(bitmap + bmps[i] + offs[i*8 + j]) |= 1 << j;
   }

   /*** Main algorithm end ***/

   /* First byte gets clobbered for the first block */
   bitmap[0] = 0x1;

   free(as);
   free(bmps);
   free(offs);
}


static uint64_t
count_primes(char *bitmap, uint64_t max_bit) {
   uint64_t *p = (uint64_t *)bitmap,
             c = max_bit / 64,
         count = max_bit + 1;

   while (c--)
      count -= __builtin_popcountll(*p++);

   count -= __builtin_popcountll(*p & (~0ul >> (63 - (max_bit % 64))));

   return count;
}


int
main(int argc, char *argv[])
{
   uint64_t  max = argc > 1 ? atol(argv[1]) : 1000000000ull;
   char     *bitmap = calloc(max/30 + 2*32*1024 + sqrtl(max) + 1, 1);
   uint64_t  adj_max;
   int       a_i;

   /* adjust max to a multiple of 2,3,5 */
   for (adj_max = (max - 1) / 30 * 30 + 1, a_i = 0;
        adj_max + diff[a_i%8] <= max;
        adj_max += diff[a_i++%8])
      ;

   calc_primes(bitmap, adj_max);

   printf("%ld\n", count_primes(bitmap, adj_max / 30 * 8 + num2bit(adj_max)) + 3);

   free(bitmap);
   return EXIT_SUCCESS;
}
