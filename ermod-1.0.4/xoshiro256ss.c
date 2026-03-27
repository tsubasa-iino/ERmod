/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)
    Modified in 2024 by Shun Sakuraba

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide. This software is distributed without any warranty.

See <http://creativecommons.org/publicdomain/zero/1.0/>. */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>

/* for time and gettimeofday*/
#include <time.h>
#include <sys/time.h>

/* for getpid */
#include <sys/types.h>
#include <unistd.h>

/* This is xoshiro256** 1.0, one of our all-purpose, rock-solid
   generators. It has excellent (sub-ns) speed, a state (256 bits) that is
   large enough for any parallel application, and it passes all tests we
   are aware of.

   For generating just floating-point numbers, xoshiro256+ is even faster.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. */

static inline uint64_t rotl(const uint64_t x, int k) {
	return (x << k) | (x >> (64 - k));
}

typedef struct {
    uint64_t s[4];
} xoshiro256ss_state;

uint64_t xoshiro256ss_next(xoshiro256ss_state *st);

static uint64_t splitmix64_next(uint64_t *x) {
	uint64_t z = (*x += 0x9e3779b97f4a7c15);
	z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
	z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
	return z ^ (z >> 31);
}

xoshiro256ss_state* xoshiro256ss_init_state_with_seed(uint64_t seed)
{
    xoshiro256ss_state *st;
    st = malloc(sizeof(xoshiro256ss_state));
    st->s[0] = splitmix64_next(&seed);
    st->s[1] = splitmix64_next(&seed);
    st->s[2] = splitmix64_next(&seed);
    st->s[3] = splitmix64_next(&seed);

    /* warmup */
    xoshiro256ss_next(st);
    xoshiro256ss_next(st);
    xoshiro256ss_next(st);
    xoshiro256ss_next(st);
    return st;
}

xoshiro256ss_state* xoshiro256ss_init_state_with_seeds(uint64_t seed[4])
{
    xoshiro256ss_state *st;
    st = malloc(sizeof(xoshiro256ss_state));
    st->s[0] = seed[0];
    st->s[1] = seed[1];
    st->s[2] = seed[2];
    st->s[3] = seed[3];
    return st;
}

void xoshiro256ss_get_internal_seeds(xoshiro256ss_state* st, uint64_t internal[4])
{
    internal[0] = st->s[0];
    internal[1] = st->s[1];
    internal[2] = st->s[2];
    internal[3] = st->s[3];
}

/* Unused in ERmod */
xoshiro256ss_state* xoshiro256ss_init_state(void)
{
    xoshiro256ss_state *ret;
    uint64_t seed = 0;
    uint64_t seeds[4];
    int t;
    FILE *fp;
    size_t rsize;
    struct timeval tv;

    /* exists on Linux, MacOSX, Cygwin. Cygwin /dev/urandom falls back to LCG randgen though */
    fp = fopen("/dev/urandom", "rb");
    rsize = fread(seeds, sizeof(*seeds), sizeof(seeds) / sizeof(*seeds), fp);
    /* if fread should fail, it'll be OK because of additional entropy below (plus uninitialized seeds)*/
    (void) rsize;
    fclose(fp);

    t = time(NULL);
    gettimeofday(&tv, NULL);

    /*  values below may be equal if there are subsequent calls to this function.
        by taking xor with the pointer the possibility is minimized (unless free'ed asap) */
    ret = malloc(sizeof(xoshiro256ss_state));
    seeds[0] ^= (uint64_t) ret;
    ret->s[0] = seeds[0];
    seeds[1] ^= getpid();
    ret->s[1] = seeds[1];
    seeds[2] ^= (uint64_t) t;
    ret->s[2] = seeds[2];
    seeds[3] ^= (uint64_t) tv.tv_usec;
    ret->s[3] = seeds[3];

    /* warmup */
    xoshiro256ss_next(ret);
    xoshiro256ss_next(ret);
    xoshiro256ss_next(ret);
    xoshiro256ss_next(ret);

    return ret;
}

void xoshiro256ss_fini(xoshiro256ss_state* st)
{
    if(st) free(st);
    return;
}

uint64_t xoshiro256ss_next(xoshiro256ss_state *st)
{
    uint64_t *s = st->s;
	const uint64_t result = rotl(s[1] * 5, 7) * 9;

	const uint64_t t = s[1] << 17;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 45);

	return result;
}

/* This is the jump function for the generator. It is equivalent
   to 2^128 calls to next(); it can be used to generate 2^128
   non-overlapping subsequences for parallel computations. */

void xoshiro256ss_jump(xoshiro256ss_state *st)
{
    uint64_t *s = st->s;
	static const uint64_t JUMP[] = { 0x180ec6d33cfd0aba, 0xd5a61266f0c9392c, 0xa9582618e03fc9aa, 0x39abdc4529b1661c };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof JUMP / sizeof *JUMP; i++) {
		for(int b = 0; b < 64; b++) {
			if (JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			xoshiro256ss_next(st);
		}
    }
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

/* This is the long-jump function for the generator. It is equivalent to
   2^192 calls to next(); it can be used to generate 2^64 starting points,
   from each of which jump() will generate 2^64 non-overlapping
   subsequences for parallel distributed computations. */

void xoshiro256ss_long_jump(xoshiro256ss_state *st)
{
    uint64_t *s = st->s;
	static const uint64_t LONG_JUMP[] = { 0x76e15d3efefdcbbf, 0xc5004e441c522fb3, 0x77710069854ee241, 0x39109bb02acbe635 };

	uint64_t s0 = 0;
	uint64_t s1 = 0;
	uint64_t s2 = 0;
	uint64_t s3 = 0;
	for(int i = 0; i < sizeof LONG_JUMP / sizeof *LONG_JUMP; i++) {
		for(int b = 0; b < 64; b++) {
			if (LONG_JUMP[i] & UINT64_C(1) << b) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			xoshiro256ss_next(st);
		}
    }
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

double xoshiro256ss_next_double(xoshiro256ss_state *st)
{
    uint64_t v;
    double ret;
    v = xoshiro256ss_next(st);
    ret = (v >> 11) * 0x1.0p-53;
    return ret;
}

/* returns 0 <= x < 2^31 */
int32_t xoshiro256ss_next_int31(xoshiro256ss_state *st)
{
    uint64_t v;
    v = xoshiro256ss_next(st);
    return (int32_t)(v & 0x7fffffffULL);
}



