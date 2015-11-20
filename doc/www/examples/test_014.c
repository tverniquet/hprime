#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

const char diff[] = {6,4,2,4,2,4,6,2};

static inline int num2bit(uint64_t num) { return (num % 30) * 8 / 30; }

/*** Main algorithm start ***/

static inline void __attribute__((always_inline))
mark_multiples (char *bitmap, char *bitmap_end, uint64_t a, int jv)
{
   char *bmp;
   int   i,
         iv;

   bmp = bitmap + a * (a/30) + a/30;
   while (bmp < bitmap_end) {
      for (iv = 1, i = 0; i < 8; iv += diff[i++]) {
         *bmp |= 1 << num2bit(iv * jv);
          bmp += a/30 * diff[i]  +  (((iv + diff[i]) * jv) / 30 - (iv * jv) / 30);
      }
   }
}

/*** Main algorithm end ***/

static void
calc_primes(char *bitmap, uint64_t max)
{
   char     *bitmap_end = bitmap + max/30 + 1;
   uint64_t  a;
   int       a_i;

   /*** Main algorithm start ***/

   for (a = 7, a_i = 1; a <= sqrtl(max); a += diff[a_i++%8]) {
      if (bitmap[a / 30] & 1 << (a_i%8))
         continue;

      switch (a_i % 8) {
         case 0 : mark_multiples(bitmap, bitmap_end, a,  1); break;
         case 1 : mark_multiples(bitmap, bitmap_end, a,  7); break;
         case 2 : mark_multiples(bitmap, bitmap_end, a, 11); break;
         case 3 : mark_multiples(bitmap, bitmap_end, a, 13); break;
         case 4 : mark_multiples(bitmap, bitmap_end, a, 17); break;
         case 5 : mark_multiples(bitmap, bitmap_end, a, 19); break;
         case 6 : mark_multiples(bitmap, bitmap_end, a, 23); break;
         case 7 : mark_multiples(bitmap, bitmap_end, a, 29); break;
      }
   }

   /*** Main algorithm end ***/

   /* Reset first byte clobbered above */
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
   char     *bitmap = calloc(max/30 + sqrtl(max) + 1, 1);
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
