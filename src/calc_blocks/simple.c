#include <stdlib.h>
#include <string.h>

#include "misc.h"
#include "wheel.h"
#include "ctx.h"

struct simple_ctx
{
   uint32_t  start_prime;
   uint32_t  end_prime;
   uint32_t *primelist;
   uint32_t  primelist_count;
};


int
simple_init(struct prime_ctx *pctx, uint32_t start_prime, uint32_t end_prime, void **ctx)
{
   struct simple_ctx *sctx = malloc(sizeof (struct simple_ctx));
   *ctx = sctx;

   sctx->start_prime = start_prime;
   sctx->end_prime = MIN(pctx->run_info.max_sieve_prime, end_prime);

   sctx->primelist = malloc(sizeof(uint32_t) * (sctx->end_prime - start_prime) / 4 + 1000);
   sctx->primelist_count = 0;

   return 0;
}


int
simple_free(void *ctx)
{
   struct simple_ctx *sctx = ctx;
   FREE(sctx->primelist);
   FREE(ctx);
   return 0;
}


int
simple_add_sieving_primes(uint32_t *primelist, uint32_t *ind, uint32_t size, void *ctx)
{
   struct simple_ctx *sctx = ctx;
   for (; *ind < size; (*ind)++) {
      if (primelist[*ind] < sctx->start_prime)
         continue;
      if (primelist[*ind] > sctx->end_prime)
         return 0;

      sctx->primelist[sctx->primelist_count++] = primelist[*ind];
   }
   return 0;
}


int
simple_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx)
{
   (void)pctx;
   (void)target_num;
   (void)ctx;

   return 0;
}


/*
 * This does not overrun the buffer so should be safe for larger primes.
 * However, it will not be efficient for large primes due to all of the checking.
 *
 * No testing has been done for really large numbers which may affect the
 * pointer arithmetic.
 *
 * No state is stored which means it is inherently safe for threaded
 * calculations.
 */
static void
mark_off_prime(struct prime_current_block *pcb, uint32_t sieve_prime)
{
   char     *bmp;
   int       i;
   int       offsets[8];
   const unsigned char bits[]  = {0x01,0x02,0x04,0x08,0x10,0x20,0x40,0x80};

   for (i = 0; i < 8; i++)
      offsets[pp_to_bit(ind_to_mod[i] * sieve_prime)] = num_to_bytes(ind_to_mod[i] * sieve_prime);

   bmp = pcb_initial_offset(pcb, sieve_prime);

   /*
    * Only perform range checks for the first and last set of 8
    */

   for (i = 0; i < 8; i++)
      if (pcb_inrange(pcb, bmp + offsets[i]))
         *(bmp + offsets[i]) |= bits[i];

   if ((bmp += sieve_prime) > pcb_end(pcb))
      return;

   for ( ; bmp < pcb_end(pcb) - sieve_prime; bmp += sieve_prime)
      for (i = 0; i < 8; i++)
         *(bmp + offsets[i]) |= bits[i];

   for (i = 0; i < 8; i++)
      if (pcb_inrange(pcb, bmp + offsets[i]))
         *(bmp + offsets[i]) |= bits[i];
}


int
simple_calc_primes(struct prime_thread_ctx *ptx, void *ctx)
{
   struct simple_ctx *sctx = ctx;
   uint32_t i;

   if (sctx->primelist_count == 0)
      return 0;

   if (sctx->primelist[0] == 7)
      memset(ptx->current_block.block, 0, ptx->current_block.block_size);

   for (i = 0; i < sctx->primelist_count; i++) {
      if (sctx->primelist[i] > ptx->current_block.sqrt_end_num)
         break;
      mark_off_prime(&ptx->current_block, sctx->primelist[i]);
   }
   return 0;
}
