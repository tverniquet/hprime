#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>

#define MIN(A,B) ((A) < (B) ? (A) : (B))

static void
calc_primes(char *bitmap, uint64_t max)
{
   uint32_t *as = malloc(sizeof (uint32_t) * (sqrtl(max) / 4 + 1000));
   uint64_t *rs = malloc(sizeof (uint64_t) * (sqrtl(max) / 4 + 1000));
   int       cnt = 0,
             i;
   uint64_t  a,
             r,
             cur;

   /*** Main algorithm start ***/

   for (a = 2; a <= sqrtl(max); a++) {
      if (bitmap[a / 8] & 1 << (a % 8))
         continue;

      for (r = a * a; r < sqrtl(max); r += a)
         bitmap[r / 8] |= 1 << (r % 8);

      as[cnt] = a; rs[cnt] = r; cnt++;
   }


   for (cur = 0; cur <= max; cur += 32*1024*8) {

      for (i = 0; i < cnt; i++)
         for (; rs[i] < cur + 32*1024*8; rs[i] += as[i])
            bitmap[rs[i]/8] |= 1 << (rs[i] % 8);

   }

   /*** Main algorithm end ***/

   /* mark bit 0 and 1 as not prime */
   bitmap[0] |= 0x3;
   free(as);
   free(rs);
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
   char     *bitmap = calloc(1, max/8 + 32*1024);

   calc_primes(bitmap, max);

   printf("%ld\n", count_primes(bitmap, max));

   free(bitmap);
   return EXIT_SUCCESS;
}
