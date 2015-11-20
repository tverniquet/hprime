#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>


/*** Main algorithm start ***/
static inline void __attribute__((always_inline))
mark_multiples(char *bitmap, uint64_t max, uint64_t a, int j)
{
   char *bmp;
   int   i;

   for (bmp = bitmap + a * (a / 8); bmp <= bitmap + max/8; )
      for (i = 0; i < 8; i++, bmp += a/8 + ((i * j / 8) - (i-1) * j / 8))
         *bmp |= 1 << ((i * j) % 8);
}
/*** Main algorithm end ***/


static void
calc_primes(char *bitmap, uint64_t max)
{
   uint64_t a;

   /*** Main algorithm start ***/

   for (a = 2; a <= sqrtl(max); a++) {
      if (bitmap[a / 8] & 1 << (a % 8))
         continue;

      switch (a % 8) {
         case 0: mark_multiples(bitmap, max, a, 0); break;
         case 1: mark_multiples(bitmap, max, a, 1); break;
         case 2: mark_multiples(bitmap, max, a, 2); break;
         case 3: mark_multiples(bitmap, max, a, 3); break;
         case 4: mark_multiples(bitmap, max, a, 4); break;
         case 5: mark_multiples(bitmap, max, a, 5); break;
         case 6: mark_multiples(bitmap, max, a, 6); break;
         case 7: mark_multiples(bitmap, max, a, 7); break;
      }
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
