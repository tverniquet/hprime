#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include <math.h>

#define MAX(A,B) ((A) > (B) ? (A) : (B))
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#define CEIL_TO(A, B) (((A) + (B) - 1) / (B) * (B))

struct ctx {
   char     *bitmap;
   uint64_t  max;
   uint64_t  cur;
   uint64_t  next;

   uint32_t *as;
   uint64_t *rs;
   int       cnt;
};


static void
init_ctx (struct ctx *ctx, uint64_t max)
{
   uint32_t a;
   uint32_t r;

   ctx->bitmap = calloc(1, 32*1024 );
   ctx->as = malloc(sizeof (uint32_t) * (sqrtl(max) / 4 + 1000));
   ctx->rs = malloc(sizeof (uint64_t) * (sqrtl(max) / 4 + 1000));
   ctx->max = max;
   ctx->cur = 0;
   ctx->next = 0;
   ctx->cnt = 0;

   /* Bootstrap primes for the first block */
   for (a = 2; a <= sqrtl(32*1024*8); a++) {
      if (ctx->bitmap[a / 8] & 1 << (a % 8))
         continue;

      for (r = a * a; r < sqrtl(32*1024*8); r += a)
         ctx->bitmap[r / 8] |= 1 << (r % 8);

      ctx->as[ctx->cnt] = a; ctx->rs[ctx->cnt] = r; ctx->cnt++;
   }

   ctx->bitmap[0] |= 0x3;
}



static int
calc_block(struct ctx *ctx)
{
   int       i;
   char     *bitmap;
   uint32_t  a;

   ctx->cur = ctx->next;
   ctx->next += 32*1024*8;

   if (ctx->cur > ctx->max)
      return 0;

   if (ctx->cur != 0)
      bzero(ctx->bitmap, 32*1024);

   bitmap = ctx->bitmap - ctx->cur/8;

   /*** Main algorithm start ***/

   for (i = 0; i < ctx->cnt; i++)
      for (; ctx->rs[i] < ctx->next; ctx->rs[i] += ctx->as[i])
         bitmap[ctx->rs[i] / 8] |= 1 << (ctx->rs[i] % 8);

   for (a = MAX(ctx->cur,ctx->as[ctx->cnt-1]); a < MIN(sqrtl(ctx->max), ctx->cur + 32*1024*8); a++) {
      if (bitmap[a / 8] & 1 << (a % 8))
         continue;

      ctx->as[ctx->cnt] = a; ctx->rs[ctx->cnt] = a * a; ctx->cnt++;
   }

   /*** Main algorithm end ***/

   return 1;
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
   uint64_t  max = argc > 1 ? atol(argv[1]) : 1000000000ul;
   uint64_t  num_primes = 0ul;
   struct ctx ctx;

   /*** Main algorithm start ***/

   init_ctx(&ctx, max);

   while (calc_block(&ctx))
      num_primes += count_primes(ctx.bitmap, MIN(32*1024*8 - 1, max - ctx.cur));

   /*** Main algorithm end ***/

   printf("%ld\n", num_primes);

   return EXIT_SUCCESS;
}
