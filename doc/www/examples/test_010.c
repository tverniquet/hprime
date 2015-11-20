#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>


static void
calc_primes(char *bitmap, uint64_t max)
{
   uint64_t a,
            r;

   /*** Main algorithm start ***/

   for (a = 3; a <= sqrtl(max); a += 2) {
      if (bitmap[a/2/8] & 1 << ((a/2) % 8))
         continue;

      for (r = a * a; r <= max; r += a * 2)
         bitmap[r/2/8] |= 1 << ((r/2) % 8);
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
   char     *bitmap = calloc(max/8/2 + 1, 1);

   calc_primes(bitmap, max);

   printf("%ld\n", count_primes(bitmap, (max-1)/2 + 1));

   free(bitmap);
   return EXIT_SUCCESS;
}
