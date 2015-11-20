#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

const char diff[] = {6,4,2,4,2,4,6,2};
const char bval[] = {1,7,11,13,17,19,23,29};

static inline int num2bit(uint64_t num) { return (num % 30) * 8 / 30; }


static void
calc_primes(char *bitmap, uint64_t max)
{
   char     *bmp,
            *bitmap_end = bitmap + max / 30 + 1;
   uint64_t  a;
   int       i,
             a_i,
             offsets[8];
   char      masks[8];

   /*** Main algorithm start ***/

   for (a = 7, a_i = 1; a <= sqrtl(max); a += diff[a_i++%8]) {
      if (bitmap[a / 30] & 1 << (a_i % 8))
         continue;

      for (i = 0; i < 8; i++) {
         offsets[i] = a * bval[i] / 30;
         masks[i] = 1 << num2bit(a * bval[i]);
      }

      for (bmp = bitmap + a * (a / 30);  bmp < bitmap_end; bmp += a)
         for (i = 0; i < 8; i++)
            *(bmp + offsets[i]) |= masks[i];
   }

   /*** Main algorithm end ***/

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
   char     *bitmap = calloc(1, max/30 + sqrtl(max) + 1);
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
