#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>


const char diff[] = {6,4,2,4,2,4,6,2};

static void
calc_primes(char *bitmap, uint64_t max)
{
   uint64_t  a,
             r;
   int       a_i,
             b_i;

   /*** Main algorithm start ***/

   for (a = 7, a_i = 1; a <= sqrtl(max); a += diff[a_i++%8]) {
      if (bitmap[a / 30] & 1 << ((a % 30) * 8 / 30))
         continue;

      for (r = a * a, b_i = a_i; r <= max; r += a * diff[b_i++%8])
         bitmap[r / 30] |= 1 << ((r % 30) * 8 / 30);
   }

   /*** Main algorithm end ***/

   /* mark 1 as not prime */
   bitmap[0] |= 0x1;
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
   char     *bitmap = calloc(max/30 + 1, 1);
   uint64_t  adj_max;
   int       a_i;

   /* adjust max to a multiple of 2,3,5 */
   for (adj_max = (max - 1) / 30 * 30 + 1, a_i = 0;
        adj_max + diff[a_i%8] <= max;
        adj_max += diff[a_i++%8])
      ;

   calc_primes(bitmap, adj_max);

   printf("%ld\n", count_primes(bitmap, adj_max / 30 * 8 + (adj_max % 30) * 8 / 30) + 3);

   free(bitmap);
   return EXIT_SUCCESS;
}
