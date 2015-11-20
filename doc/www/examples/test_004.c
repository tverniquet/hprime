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

      for (bmp = bitmap + a * (a / 8);  bmp < bitmap_end; bmp += a)
         for (i = 0; i < 8; i++)
            *(bmp + (a * i / 8)) |= 1 << (a * i % 8);
   }

   /*** Main algorithm end ***/

   /* Reset the first byte as it is clobbered by the above algorithm*/
   bitmap[0] = 0x53;
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
   char     *bitmap = calloc(max/8 + sqrtl(max) + 1, 1);

   calc_primes(bitmap, max);

   printf("%ld\n", count_primes(bitmap, max));

   free(bitmap);
   return EXIT_SUCCESS;
}
