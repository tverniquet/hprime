#ifndef _HARU_PLAN_REGISTER_H
#define _HARU_PLAN_REGISTER_H

#include "plans.h"


/**
 * @FILE - Register the plans to use to calculate the blocks of primes
 */


#define DECLARE_PLAN_ENTRY_FUNCTIONS(name) \
   int name##_init(struct prime_ctx *pctx, uint32_t start_sieve_prime, uint32_t end_sieve_prime, void **ctx); \
   int name##_free(void *ctx); \
   int name##_skip_to(struct prime_thread_ctx *pctx, uint64_t target_num, void *ctx); \
   int name##_add_sieving_primes(uint32_t *primes, uint32_t *ind, uint32_t max_ind, void *ctx); \
   int name##_calc_primes(struct prime_thread_ctx *, void *ctx)


#define USE_PLAN_ENTRY_FUNCTIONS(name) \
   name##_init, \
   name##_free, \
   name##_skip_to, \
   name##_add_sieving_primes, \
   name##_calc_primes

/**
 * slow - The slow method:
 *
 * This performs "a x b" and marks the resulting multiple as "not-prime"
 * within the block.
 *
 * 'b' is incremented according to the "possible primes" in the 8/30 wheel. That
 * is, b is never a multiple of 2,3,5.
 *
 * psudo code:
 *
 * For each block:
 *   For each 'a' (ie each prime to sqrt(block_end)):
 *      b = max( block_start / a, a);
 *      For each b where 'a x b' is in the block and b is a wheel number
 *         Mark a x b as not-prime
 */
DECLARE_PLAN_ENTRY_FUNCTIONS(slow);


/**
 * simple - Direct ("simple") implementation of some major optimisations
 *
 * The main optimisation in this implemention compared to the "slow"
 * implentation is in avoiding:
 *   - The multiply needed to calculate 'a * b'
 *   - The calculations required each time to find the resulting byte/bit to
 *     mark "not-prime"
 *
 * This is done by noting:
 *   a * b = (a_byte + a_bit) * (b_byte + b_bit)
 *
 *         =   a_byte * b_byte
 *           + a_byte * b_bit
 *           + a_bit  * b_byte
 *           + a_bit  * b_bit  *** influences bit position of result
 *
 * From here, only the last term influences which 'bit' to set. Ie, all of the
 * first terms are multiples of "byte". (In the wheel a "byte" is just 30, but
 * also corresponds directly to a byte in the wheel storage)
 *
 * The optimisation is to choose 'b_bit', but then only increment b_byte while
 * marking off possible 'a * b' within the block. This means that the resulting
 * bit mask does not change for that run as a_bit and b_bit are constant.  This
 * step is then done for each of the 8 possible 'b_bit'.
 *
 * We can continue changing the above function to note:
 *
 *   a * b = a * (b_byte + b_bit) = a *  b_byte + a * b_bit
 *
 * So if 'a' is constant and 'b_bit' is constant, then incrementing
 * b_byte means that the next result can be achieved by adding 'a' bytes
 * to the result.
 *
 *
 * So finally, the psudocode is:
 *
 *   For each block:
 *      For each 'a' (ie each prime to sqrt(block_end))
 *         b_byte = max (block_start / a, a) / 30
 *         for each bit (0,1,2,3,4,5,6,7)
 *            b = b_byte * 30 + wheel_value_of(bit);
 *            calculate the res_byte and res_bit from a * b
 *            while res_byte is in the block:
 *               apply res_bit mask to res_byte
 *               increment res_byte by 'a'
 *
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(simple);


/**
 * simple_middle - Small advancement on the general 'simple' plan
 *
 * This only works for primes that are < block_size, and to a maximum of 2^16.
 *
 * The main optimisations from the 'simple' plan are:
 *  - calculates all 8 'b_bits' for a given prime at the same time
 *  - store the offset for these primes between blocks to save re-calculating
 *    each time
 *
 * NOTE: this requires a fair bit of space (20 bytes per prime)
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(simple_middle);


/**
 * @file load_unaligned implementation for lower primes
 *
 * Instead of shuffling the 32byte registers, an initial set
 * of buffers is first calculated for each prime.
 *
 * These are then loaded onto the stack for quick reference (did seem to make a
 * difference) ????
 *
 * Then an intrinsic is used to perform an unaligned load according to an offset
 * and the offset is shifted according to the prime. Currently this calculation
 * of the offset seems a little clunky.
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(load_unaligned);

/**
 * read_offs - Optimisation on the simple_middle plan
 *
 * This only works for primes that are < block_size, and to a maximum of 2^15.
 *
 * The main optimisations from the 'simple_middle' plan are:
 *   - Offsets are calculated initially but are used read_only. A separate
 *     offset is kept which is then applied to all bit offsets.
 *
 *     The idea was that there doesn't need to be writing back to the L2 cache
 *     (though I don't know if this even happened or if it was a bottle neck)
 *
 *     NOTE: a main advantage of this is that it will work between threads
 *
 *   - the primes overrun the buffer instead of performing > buflen check
 *   - The initial < bufstart check became a conditional += sieve_prime
 *
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(read_offs);


/**
 * calc_offs - Calculate the offsets
 *
 * This only works for primes that are < block_size, and to a maximum of 2^15.
 *
 * It calculates the 8 relative offsets for 16 primes at a time in ymm
 * registers.
 *
 * Only a single 16bit int per id is stored each time.
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(calc_offs);


/**
 * lu_calc_offs - Calculate the offsets
 *
 * This only works for primes that are > 2^15. It is not optimal
 * once no bits are set in any one 32*1024 block.
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(lu_calc_offs);


/**
 * keep_byte - instead of all offsets, keep some information and calculate 'b'
 *
 * This seems to be slower...
 */

DECLARE_PLAN_ENTRY_FUNCTIONS(keep_byte);


const struct prime_plan plans[] = {
   {
      "calc middle primes", 32*1024, 4,
         {{"unaligned",     1,          0,          96,  USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          {"calc_offs",     1,         96,     (1<<15),  USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          {"lu_calc_offs",  1,    (1<<15),      400000,  USE_PLAN_ENTRY_FUNCTIONS(lu_calc_offs)},
          {"simple_sieve",  1,     400000,  UINT32_MAX,  USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "calc lower upper primes", 32*1024, 4,
         {{"unaligned",     1,          0,          96,  USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          {"read_offs",     1,         96,     (1<<15),  USE_PLAN_ENTRY_FUNCTIONS(read_offs)},
          {"lu_calc_offs",  1,    (1<<15),      400000,  USE_PLAN_ENTRY_FUNCTIONS(lu_calc_offs)},
          {"simple_sieve",  1,     400000,  UINT32_MAX,  USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "read offset to 10^15", 32*1024, 4,
         {{"unaligned",     1,          0,          96,  USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          {"read_offs",     1,         96,     (1<<15),  USE_PLAN_ENTRY_FUNCTIONS(read_offs)},
          {"simple_middle", 1,    (1<<15),  (1<<16)-30,  USE_PLAN_ENTRY_FUNCTIONS(simple_middle)},
          {"simple_sieve",  1, (1<<16)-30,  UINT32_MAX,  USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "simple middle", 32*1024, 3,
         {{"unaligned",     1,          0,          96,  USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          {"simple_middle", 1,         96,  (1<<16)-30,  USE_PLAN_ENTRY_FUNCTIONS(simple_middle)},
          {"simple_sieve",  1,     400000,  UINT32_MAX,  USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "load unaligned", 32*1024, 2,
         {{"unaligned",     1,          0,          96,  USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          {"simple_sieve",  1,         96,  UINT32_MAX,  USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "breakdown - best", 32*1024, 14,
         {{ "to   16", 1,       0,  (1<<4), USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          { "to   32", 1,  (1<<4),  (1<<5), USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          { "to   64", 1,  (1<<5),  (1<<6), USE_PLAN_ENTRY_FUNCTIONS(load_unaligned)},
          { "to  128", 1,  (1<<6),  (1<<7), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to  256", 1,  (1<<7),  (1<<8), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to  512", 1,  (1<<8),  (1<<9), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to 1024", 1,  (1<<9), (1<<10), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to 2048", 1, (1<<10), (1<<11), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to 4096", 1, (1<<11), (1<<12), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to   8k", 1, (1<<12), (1<<13), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to  16k", 1, (1<<13), (1<<14), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to  32k", 1, (1<<14), (1<<15), USE_PLAN_ENTRY_FUNCTIONS(calc_offs)},
          { "to  64k", 1, (1<<15), (1<<16), USE_PLAN_ENTRY_FUNCTIONS(lu_calc_offs)},
          { "rest",    1, (1<<16), UINT32_MAX, USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "breakdown simple", 32*1024, 14,
         {{ "to   16", 1,       0,  (1<<4), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to   32", 1,  (1<<4),  (1<<5), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to   64", 1,  (1<<5),  (1<<6), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to  128", 1,  (1<<6),  (1<<7), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to  256", 1,  (1<<7),  (1<<8), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to  512", 1,  (1<<8),  (1<<9), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to 1024", 1,  (1<<9), (1<<10), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to 2048", 1, (1<<10), (1<<11), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to 4096", 1, (1<<11), (1<<12), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to   8k", 1, (1<<12), (1<<13), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to  16k", 1, (1<<13), (1<<14), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to  32k", 1, (1<<14), (1<<15), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "to  64k", 1, (1<<15), (1<<16), USE_PLAN_ENTRY_FUNCTIONS(simple)},
          { "rest",    1, (1<<16), UINT32_MAX, USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "simple", 32*1024, 1,
         {{ "simple_sieve", 1, 0, UINT32_MAX, USE_PLAN_ENTRY_FUNCTIONS(simple)}}
   },
   {
      "slow", 32*1024, 1,
         {{ "slow_sieve", 1, 0, UINT32_MAX, USE_PLAN_ENTRY_FUNCTIONS(slow)}}
   }
};
#endif
