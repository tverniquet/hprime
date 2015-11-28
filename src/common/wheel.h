#ifndef _HARU_WHEEL_H
#define _HARU_WHEEL_H

#include <inttypes.h>

/**
 * @FILE This file contains functions for interacting with the 8/30 wheel
 * storage block.
 *
 * The primes are discovered by marking multiples of primes as "not prime" on a
 * 8/30 wheel, one block at a time.
 */


struct prime_current_block
{
   char     *block;
   uint32_t  block_size;
   uint64_t  block_start_num;
   uint64_t  block_end_num;
   uint64_t  sqrt_end_num;

   uint64_t  block_start_byte;
   uint64_t  block_num;
};


/* Note: num is exclusive (ie 30 is 1 byte) */
uint64_t num_to_bytes(uint64_t num);

/* The bit within a byte of 'num' - only for possible primes */
uint64_t pp_to_bit(uint64_t num);

/*
 * The bit within a byte of 'num' - for all numbers.
 *
 * Note non-pp are rounded up to the nearest possible prime
 */
uint64_t num_to_bit(uint64_t num);

/* Note: This is the number 1 more than that covered by bytes (ie 1 byte covers 30 nums) */
uint64_t bytes_to_num(uint64_t bytes);


/**
 * The byte just past the end of the block
 */
char * pcb_end(struct prime_current_block *pcb);


/*
 * Checks if the offset is within the current block
 */
int pcb_inrange(struct prime_current_block *pcb, char *p);


/*
 * This calculation was done in many places and is a little messy.
 *
 * Basically all that is needed is -(block_start_byte % prime);
 *
 *            block_start_byte
 *            |-------
 * .   .   .xx
 *
 *
 *  .  - dots are multiples of a
 *  xx - block_start_byte % prime (ie the remainder)
 *
 * This returns a byte position that is aligned to 'a'. The 8 multiples within
 * the 'a' byte repeat pattern need to be checked to see if they are inside the
 * block. That is it is the initial offset for the block for that prime (if
 * using this technique)
 *
 * The only real complexity is that it is desirable to skip until a * a, so
 * this choice needs to be made (otherwise for the first block primes mark
 * themselves as multiples and not prime)
 *
 */
char * pcb_initial_offset(struct prime_current_block *pcb, uint32_t prime);


/*
 * Return the first prime from the given byte/bit offset in the current
 * block.
 *
 * The byte/bit are incremented so this this can be called in a loop to return
 * all primes in the current block, however the byte/bit don't actually refer
 * to the returned prime.
 *
 * Returns 0 if there are no more primes in the block
 */
uint64_t get_next_prime(struct prime_current_block *pcb, uint32_t *byte, uint32_t *bit);


/*
 * Given a byte and a bit, return what the number that the possible prime represents
 */
uint64_t get_possible_prime(uint64_t *byte, uint32_t *bit);


/*
 * Returns the next 'possible prime' in the 8/30 wheel.
 *
 * The byte/bit will be set according to the returned possible prime
 */
uint64_t get_next_possible_prime(uint64_t *byte, uint32_t *bit);

/*
 * Returns the previous 'possible prime' in the 8/30 wheel.
 *
 * The byte/bit will be set according to the returned possible prime
 */
uint64_t get_previous_possible_prime(uint64_t *byte, uint32_t *bit);


/*
 * Returns the next "result set" in the 8/30 wheel.  That is, base_byte and
 * b_bit will point to the next values in the result set.
 */
void get_next_result_set(uint32_t a_byte, uint32_t a_bit, uint64_t *base_byte, uint32_t *b_bit);


const unsigned char a_x_b_bytes[8][8]   __attribute__((aligned(32)));
const unsigned char a_x_b_byte_diffs[8][8] __attribute__((aligned(32)));
const unsigned char a_x_b_bits[8][8]    __attribute__((aligned(32)));
const unsigned char a_x_b_bitmask[8][8] __attribute__((aligned(32)));
const unsigned char a_x_b_bitvals[8][8] __attribute__((aligned(32)));

const unsigned int pp_diffs[8] __attribute__((aligned(32)));

const char ind_to_mod[8] __attribute__((aligned(32)));
const char mod_to_ind[30] __attribute__((aligned(32)));

#endif
