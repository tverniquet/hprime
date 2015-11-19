#ifndef _HARU_PRIME_H
#define _HARU_PRIME_H

struct prime_ctx;
struct prime_thread_ctx;

int calc_next_block (struct prime_ctx *ctx);

int calc_blocks (struct prime_ctx *ctx, int (*fn)(struct prime_thread_ctx *, void *), void *);

/*
 * Test a block
 */
void slow_test_primes(struct prime_ctx *pctx);

#endif

