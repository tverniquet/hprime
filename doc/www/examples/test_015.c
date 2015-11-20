#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

const char diff[] = {6,4,2,4,2,4,6,2};
const char bval[] = {1,7,11,13,17,19,23,29};

static inline int num2bit(uint64_t num) { return (num % 30) * 8 / 30; }


static void
calc_primes(char *bitmap, uint64_t max)
{
   char     *bmp,
            *bitmap_end = bitmap + max / 30 + 1,
             pattern[16 * 8];
   uint64_t  a;
   int       i,
             j,
             r,
             a_i,
             offsets[8];

   /*** Main algorithm start ***/

   for (i = 1; i < 4; i++) {
      bzero(pattern, sizeof pattern);
      for (r = bval[i], a_i = 0; r < bval[i] * 30 * 8; r += bval[i] * diff[a_i++ % 8])
         pattern[r / 30] |= 1 << num2bit(r);

      for (bmp = bitmap; bmp < bitmap_end; bmp += bval[i] * sizeof(uint64_t))
         for (j = 0; j < bval[i]; j++)
            *((uint64_t *)bmp + j) |= *((uint64_t *)pattern + j);
   }

   /*** Main algorithm end ***/

   for (a = 17, a_i = 4; a <= sqrtl(max); a += diff[a_i++%8]) {
      if (bitmap[a / 30] & 1 << (a_i%8))
         continue;

      for (i = 0; i < 8; i++)
         offsets[num2bit(a * bval[i])] = a * bval[i] / 30;

      for (bmp = bitmap + a * (a / 30); bmp < bitmap_end; bmp += a)
         for (i = 0; i < 8; i++)
            *(bmp + offsets[i]) |= 1 << i;
   }

   /* Fix first byte which is clobbered above */
   bitmap[0] = 0x1;
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
   char     *bitmap = calloc(max/30 + sqrtl(max) + 16*8, 1);
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
