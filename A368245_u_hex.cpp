#include <cstdio>
#include <cstdlib>
#include <cmath>
//#include <immintrin.h>
//#include <smmintrin.h>
#include <cstdint>
#include <bit>
#include <gmp.h>
#include <thread>
#include <vector>

// g++ A368245_u_simd_hex.cpp -o sequence_SIMD -Wall -Wextra -fuse-ld=lld -Wshadow -g -O3 -std=c++20  -march=native -ffast-math -lgmp

// The A368245-like sequence but in hexadecimal base

// With test, we can see that the interesting range is +-100 around the value
// and it forms a linear function y = 90/1000 x
// The code is optimized to skip the unuseful zones
#define COEFF 13.85
#define WIDTH 50 // Width of the band

/*typedef union alignas(32)
{
	__m256i m256;

	// In little-endian

	// To get uint4_t:
	// uint8_t high = u.u8[i] >> 4;
	// uint8_t low  = u.u8[i] & 0xF;
	uint8_t u8[32];
	uint16_t u16[16];
	uint32_t u32[8];
	uint64_t u64[4];

} u256;*/


unsigned long long sum_hex_digits_power(
	unsigned long long base,
	unsigned short exp)
{
	mpz_t n;
	mpz_init(n);
	mpz_ui_pow_ui(n, base, exp);

	unsigned long long total = 0;
	size_t count = mpz_sizeinbase(n, 16); // Number of hex digits

	// Mask-and-shift for fastest hex digit extraction
	for (size_t i = 0; i < count; i++) {
		unsigned int digit = mpz_tdiv_ui(n, 16); // Get last digit mod 16
		total += digit;
		mpz_fdiv_q_2exp(n, n, 4); // Divide by 16
	}

	mpz_clear(n);
	return total;
}

/**
 * FASTEST 256-bit × 256-bit → 256-bit multiply (mod 2^256)
 *
 * Algorithm: Schoolbook multiplication with 4x4 partial products
 *
 * Performance characteristics:
 * - 10 MULX instructions (1 per term in lower triangle)
 * - ~13-16 cycles latency on Zen 3/4, Ice Lake (MULX = 3 cycles, can pipeline 4-5)
 * - Zero dependencies within each output limb computation
 * - All 4 MULX units can operate in parallel
 */
/*u256 mul256(u256 a, u256 b)
{
	u256 r;
	r.m256 = _mm256_set1_epi64x(0);
	unsigned long long hi, lo;

	// Use unsigned __int128 for clean, compiler-optimizable carry handling
	unsigned __int128 acc;

	// ========== OUTPUT LIMB 0 ==========
	// r[0] = (a[0] * b[0]) & 0xFFFFFFFFFFFFFFFF
	lo = _mulx_u64(a.u64[0], b.u64[0], &hi);
	acc = ((unsigned __int128)hi << 64) | lo;
	r.u64[0] = (uint64_t)acc;
	acc >>= 64; // carry to next limb

	// ========== OUTPUT LIMB 1 ==========
	// r[1] = (a[0]*b[1] + a[1]*b[0] + carry) & 0xFFFFFFFFFFFFFFFF
	lo = _mulx_u64(a.u64[0], b.u64[1], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	lo = _mulx_u64(a.u64[1], b.u64[0], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	r.u64[1] = (uint64_t)acc;
	acc >>= 64;

	// ========== OUTPUT LIMB 2 ==========
	// r[2] = (a[0]*b[2] + a[1]*b[1] + a[2]*b[0] + carry) & 0xFFFFFFFFFFFFFFFF
	lo = _mulx_u64(a.u64[0], b.u64[2], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	lo = _mulx_u64(a.u64[1], b.u64[1], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	lo = _mulx_u64(a.u64[2], b.u64[0], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	r.u64[2] = (uint64_t)acc;
	acc >>= 64;

	// ========== OUTPUT LIMB 3 ==========
	// r[3] = (a[0]*b[3] + a[1]*b[2] + a[2]*b[1] + a[3]*b[0] + carry) & 0xFFFFFFFFFFFFFFFF
	lo = _mulx_u64(a.u64[0], b.u64[3], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	lo = _mulx_u64(a.u64[1], b.u64[2], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	lo = _mulx_u64(a.u64[2], b.u64[1], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	lo = _mulx_u64(a.u64[3], b.u64[0], &hi);
	acc += ((unsigned __int128)hi << 64) | lo;

	r.u64[3] = (uint64_t)acc;

	// We discard the carry beyond bit 256 (the high bits of acc)

	return r;
}

// Same here but its a pow
static inline u256 pow_ull_to_m256(uint64_t a, uint16_t s)
{
	u256 base;
	base.m256 = _mm256_set_epi64x(0, 0, 0, a);
	u256 result;
	result.m256 = _mm256_set_epi64x(0, 0, 0, 1);

	while (s)
	{
		if (s & 1)
			result = mul256(result, base);
		s >>= 1;
		if (!s)
			break;
		base = mul256(base, base);
	}

	return result;
}

// max return = 15*64 = 960 (64 hex digits, 15 (F) is the max for each one)
// Fits on a u16
// Change this to get an other base
#define ABS_MAX_SUM 960ull
// Faster method than the one below
static inline int16_t sum_digits_hex(u256 number)
{
	const __m256i lut = _mm256_setr_epi8(
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,
		0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15);

	__m256i lo = _mm256_and_si256(number.m256, _mm256_set1_epi8(0x0F));
	__m256i hi = _mm256_and_si256(_mm256_srli_epi16(number.m256, 4), _mm256_set1_epi8(0x0F));

	lo = _mm256_shuffle_epi8(lut, lo);
	hi = _mm256_shuffle_epi8(lut, hi);

	__m256i sum = _mm256_add_epi8(lo, hi);

	__m256i sad = _mm256_sad_epu8(sum, _mm256_setzero_si256());

	__m128i low = _mm256_castsi256_si128(sad);
	__m128i high = _mm256_extracti128_si256(sad, 1);

	__m128i total = _mm_add_epi64(low, high);

	return (int16_t)_mm_cvtsi128_si64(total) +
		   (int16_t)_mm_extract_epi64(total, 1);
}*/

/*
static inline int16_t sum_digits_16(u256 number)
{
	// Goal: sum all u4 (hex) digits.
	// Each group of 4 bits represent a hex digit, there is 64 digits
	// ABCD EFGH abcd efgh ... (32 bytes) ... stuv wxyz
	//   \   /     \   /                        \   /
	//  (1)+(2)   (3)+(4)  ... (32  sums) ... (63)+(64)
	//        \   /                               |
	//      1+2+...+16     ... (4   sums) ... 49+...+64
	//             \          |          |          /
	//            ------------ Final  sum ------------

	// To do it:
	// abcdefgh
	// first calc the mask for the first u4 of each byte then the second half
	// 0000efgh and abcd0000
	// then sum them

	__m256i maskLow  = _mm256_set1_epi8(0x0F);

	__m256i digitsSummed = _mm256_adds_epu8( // Sum the 2 hex digits
		_mm256_and_si256(number.m256, maskLow), // the 0000efgh part
		_mm256_srli_epi16(number.m256, 4)); // The 0000abcd part

	// Now we have 32 sums of the 64 digits, so lets add the 32 bytes together
	// This instruction creates 4 u16 (stored ion u64) with the sums
	__m256i bytesSummed = _mm256_sad_epu8(digitsSummed, _mm256_setzero_si256());

	__m128i low  = _mm256_castsi256_si128(bytesSummed);
	__m128i high = _mm256_extracti128_si256(bytesSummed, 1);

	// Contains 2 numbers, so lets sum them
	__m128i sum = _mm_add_epi64(low, high);

	int16_t sLow = (int16_t)_mm_cvtsi128_si64(sum);
	int16_t sHigh = (int16_t)_mm_extract_epi64(sum, 1);

	return (int16_t)(sLow + sHigh);
}*/

// Single CPU instruction to approximate log2
inline int fast_log2(unsigned long long x)
{
	return std::bit_width(x) - 1;
}

inline long long fast_abs(long long x)
{
	long long mask = x >> 63;
	return (x + mask) ^ mask;
}

void process_line(
	u_int8_t *line,
	const unsigned long long max,
	const unsigned short n)
{
	unsigned long long a;

	// approximate maximum digit sum difference
	//long double max_ax = 9.0L * ceill(n * log10l(max));

#ifdef SKIP_USELESS
	// Skips useless zones
	unsigned long long range_start = std::max(1, (unsigned long long)((n_ld - WIDTH) * COEFF));
	unsigned long long range_end = std::min(max, (unsigned long long)((n_ld + WIDTH) * COEFF));

	// Fill pixels before the optimized range with neutral gray
	for (a = 1; a < range_start - 1; a++)
	{
		unsigned char byte_val = 255;
		fwrite(&byte_val, 1, 1, stdout);
	}

	// Process the optimized range
	for (a = range_start;
		 a < range_end;
		 a++)
#else
	// Process all the range
	for (a = 1; a < max; a++)
#endif
	{
		// Max of a is 2^64. Max of C is 2^256. so n must be <= 4
		// Update: no more issue with GMP
		//if (fast_log2(a) * n > 256)
		//	goto early_break;

		// Check if the result is too big.
		// In that case, it will never be able to work so there is no point
		// on going any farther because every sum_digits(a^n) - n will be
		// less than a
		// TODO it seems that its not working
		/*if (a > ABS_MAX_SUM - (uint64_t)n)
			goto early_break;*/


		unsigned long long result = sum_hex_digits_power(a, n);
		long long difference = result - n - a;

		if (!difference) // We found a working case
		{
			mpz_t r;
			mpz_init(r);
			mpz_ui_pow_ui(r, a, n);
			// Print the result (decimal)
			gmp_fprintf(stderr, "%Zd - %hu\n", r, n);

			//fwrite(&zero, 1, 1, stdout);
			*line++ = (uint8_t)0;
		}
		else
		{
			// Distribution to help visualise the range
			// It's kinna ugly (only 63 levels) but its fast
			// The *16 limits greatly the range but hey,
			// its there for visualisation bot accuracy
			uint16_t y = fast_log2(fast_abs(difference)) << 4;

			// Clamp y to [0, 200] first, then add 35 to get [55, 255]
			if (y > 200) y = 200;

			uint8_t byte_val = (uint8_t)(y + 55);
			//fwrite(&byte_val, 1, 1, stdout);
			*line++ = byte_val;
		}
	}

#ifdef SKIP_USELESS
	// Fill pixels after the optimized range with neutral gray
	for (a = range_end + 1; a < max; a +)
	{
		unsigned char byte_val = 255;
		fwrite(&byte_val, 1, 1, stdout);
	}
#endif

	return;
/*
early_break:
	// Fill the rest of the line
	for (; a < max; a++)
	{
		unsigned char byte_val = 245;
		fwrite(&byte_val, 1, 1, stdout);
	}
	return;*/
}

// a.out 1000000 2
// Returns PGM in cout, and infos in cerr
int main(int argc, char **argv)
{
	if (argc != 3)
		return -1;
	const unsigned n_threads = std::max(std::thread::hardware_concurrency(),1u);
	// Make n divisible by n_threads
	const unsigned short n_max = (atoi(argv[2]) / n_threads) * n_threads;
	const unsigned long long max = (uint64_t)atoll(argv[1]);

	printf("P5\n%llu %hu\n255\n", max - 1ULL, n_max);

	// Start at pow 2: any n^0 is 1 so only 10^x flags
	// Any n^1 is n so only 1-9 digits flags
	for (unsigned short n = 2; n < n_max+2; n += n_threads)
	{
		u_int8_t line[n_threads][max]; 
		std::vector<std::thread> threads;

		for (u_int8_t i = 0; i < n_threads; i++)
		{
			// WARN: dont put max too high otherwise
			// you will need To of RAM (and RAM is expensive
			// right now, 2026-02) (yes today we muggles
			// put 16gig of RAM in our computers)
		
			threads.emplace_back(process_line, (u_int8_t*)line[i], max, i+n);
		}

		for (u_int8_t i = 0; i < n_threads; i++)
		{
			threads[i].join();
			fwrite(line[i], 1, max - 1ULL, stdout);
		}
	}

	fprintf(stderr, "\n");
	return 0;
}
