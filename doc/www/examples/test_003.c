#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>


static void
calc_primes(char *bitmap, uint64_t max)
{
   char     *bmp,
            *bitmap_end = bitmap + max/8 + 1;
   uint64_t  a;
   int       i;

   /*** Main algorithm start ***/

   for (a = 2; a <= sqrtl(max); a++) {
      if (bitmap[a / 8] & 1 << (a % 8))
         continue;

      for (i = 0; i < 8; i++)
         for (bmp = bitmap + a * (a+i) / 8;  bmp < bitmap_end; bmp += a)
            *bmp |= 1 << ((a * (a+i)) % 8);
   }

   /*** Main algorithm end ***/

   /* mark bit 0 and 1 as not prime */
   bitmap[0] |= 0x3;
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
   char     *bitmap = calloc(max/8 + 1, 1);

   calc_primes(bitmap, max);

   printf("%ld\n", count_primes(bitmap, max));

   free(bitmap);
   return EXIT_SUCCESS;
}
