#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>


static uint64_t
calc_primes(char *bitmap, uint64_t max)
{
   uint64_t a,
            r,
            num_primes = 0;

   /*** Main algorithm start ***/

   for (a = 2; a <= max; a++) {
      if (bitmap[a / 8] & 1 << (a % 8))
         continue;

      num_primes++;

      for (r = a * a; r <= max; r += a)
         bitmap[r / 8] |= 1 << (r % 8);
   }

   /*** Main algorithm end ***/

   return num_primes;
}


int
main(int argc, char *argv[])
{
   uint64_t  max = argc > 1 ? atol(argv[1]) : 1000000000ull;
   char     *bitmap = calloc(max/8 + 1, 1);

   printf("%ld\n", calc_primes(bitmap, max));

   free(bitmap);
   return EXIT_SUCCESS;
}
