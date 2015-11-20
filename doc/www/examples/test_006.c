#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>


static void
calc_primes(char *bitmap, uint64_t max)
{
   uint64_t  a,
             r,
             i,
             j;
   int       low_a[] = {2,3,5,7,11,13};
   char      pattern[16*8],
            *bmp,
            *bitmap_end = bitmap + max/8 + 1;

   /*** Main algorithm start ***/

   for (i = 0; i < 6; i++) {
      bzero(pattern, sizeof pattern);
      a = low_a[i];

      for (r = 0; r < a*64; r += a)
         pattern[r / 8] |= 1 << r % 8;

      for (bmp = bitmap; bmp < bitmap_end; bmp += a * sizeof(uint64_t))
         for (j = 0; j < a; j++)
            *((uint64_t *)bmp + j) |= *((uint64_t *)pattern + j);
   }

   /*** Main algorithm end ***/

   for (a = 17; a < sqrtl(max); a++) {
      if (bitmap[a / 8] & 1 << (a % 8))
         continue;

      for (r = a * a; r <= max; r += a)
         bitmap[r / 8] |= 1 << (r % 8);
   }


   /* Reset the first bytes as they are clobbered by the above algorithm*/
   bitmap[0] = 0x53;
   bitmap[1] = 0xD7;
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
   char     *bitmap = calloc(max/8 + 16*8 + 1, 1);

   calc_primes(bitmap, max);

   printf("%ld\n", count_primes(bitmap, max));

   free(bitmap);
   return EXIT_SUCCESS;
}
