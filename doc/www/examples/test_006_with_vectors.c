#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <inttypes.h>

/*** Main algorithm start ***/

typedef uint64_t v32_4ui __attribute__((vector_size(32)));

static inline void
mark_low_prime (char *bitmap, uint64_t max, int a)
{
   int r, n;
   char     *pattern = aligned_alloc(32, a*32);
   register v32_4ui v0, v1, v2, v3, v4, v5, v6, v7, v8, v9, va, vb, vc;
   v32_4ui *out = (v32_4ui *)bitmap;

   bzero(pattern, a*32);

   for (r = 0; r < a*32*8; r += a)
      pattern[r/8] |= 1ul << r%8;

   switch (a) {
      case 13:
         vc = *(v32_4ui *)(pattern + 32*12);
         vb = *(v32_4ui *)(pattern + 32*11);
      case 11:
         va = *(v32_4ui *)(pattern + 32*10);
         v9 = *(v32_4ui *)(pattern + 32*9);
         v8 = *(v32_4ui *)(pattern + 32*8);
         v7 = *(v32_4ui *)(pattern + 32*7);
      case 7:
         v6 = *(v32_4ui *)(pattern + 32*6);
         v5 = *(v32_4ui *)(pattern + 32*5);
      case 5:
         v4 = *(v32_4ui *)(pattern + 32*4);
         v3 = *(v32_4ui *)(pattern + 32*3);
      case 3:
         v2 = *(v32_4ui *)(pattern + 32*2);
      case 2:
         v1 = *(v32_4ui *)(pattern + 32*1);
         v0 = *(v32_4ui *)(pattern + 32*0);
   }

   n = max/8/32/a + 1;
   while (n--) {
      *out++ |= v0;
      *out++ |= v1;
      if (a == 2) continue;
      *out++ |= v2;
      if (a == 3) continue;
      *out++ |= v3;
      *out++ |= v4;
      if (a == 5) continue;
      *out++ |= v5;
      *out++ |= v6;
      if (a == 7) continue;
      *out++ |= v7;
      *out++ |= v8;
      *out++ |= v9;
      *out++ |= va;
      if (a == 11) continue;
      *out++ |= vb;
      *out++ |= vc;
   }
}

/*** Main algorithm end ***/


void
calc_primes(char *bitmap, uint64_t max)
{
   uint64_t  a,
             r;

   /*** Main algorithm start ***/

   mark_low_prime (bitmap, max, 2);
   mark_low_prime (bitmap, max, 3);
   mark_low_prime (bitmap, max, 5);
   mark_low_prime (bitmap, max, 7);
   mark_low_prime (bitmap, max, 11);
   mark_low_prime (bitmap, max, 13);

   /*** Main algorithm end ***/

   for (a = 17; a < sqrtl(max); a++) {
      if (bitmap[a/8] & 1 << (a%8))
         continue;

      for (r = a * a; r <= max; r += a)
         bitmap[r/8] |= 1 << (r%8);
   }


   /* The first 16 bytes is ruined by the pattern setting */
   bitmap[0] = 0x53;
   bitmap[1] = 0xD7;
}


uint64_t
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
   char     *bitmap = aligned_alloc(32, max/8 + 32*16*2 + 1);

   bzero(bitmap, max/8 + 32*16*2 + 1);

   calc_primes(bitmap, max);

   printf("%ld\n", count_primes(bitmap, max));

   free(bitmap);
   return EXIT_SUCCESS;
}
