#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>
#include <smmintrin.h>
#include <cstdint>
#include <bit>
#include <gmp.h>
#include <thread>
#include <vector>
#include <cstring>

// g++ A368245_u_simd_hex.cpp -o sequence_SIMD -Wall -Wextra -fuse-ld=lld -Wshadow -g -O3 -std=c++20  -march=native -ffast-math -lgmp

// The A368245-like sequence but in hexadecimal base and decimal base - 
// just adjust the compilation options
// Just adjust this flag
// "sum_dec_digits_power" or "sum_hex_digits_power" or "sum_oct_digits_power"
//#define sum_digits_power(b,e) sum_dec_digits_power(b,e)

// Dont touch below
#ifndef sum_digits_power
#	ifdef DEC
#		define sum_digits_power(b,e) sum_dec_digits_power(b,e)
#	elif defined(OCT)
#		define sum_digits_power(b,e) sum_oct_digits_power(b,e)
#	else // ifdef HEX
#		define sum_digits_power(b,e) sum_hex_digits_power(b,e)
#	endif
#endif

// With test, we can see that the interesting range is +-100 around the value
// and it forms a linear function y = 90/1000 x
// The code is optimized to skip the unuseful zones
#define COEFF 13.85
#define WIDTH 50 // Width of the band



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
	for (size_t i = 0; i < count; i++)
	{
		unsigned int digit = mpz_get_ui(n) & 0xF; // Get last digit mod 16
		total += digit;
		mpz_fdiv_q_2exp(n, n, 4); // Divide by 16
	}

	mpz_clear(n);
	return total;
}

unsigned long long sum_oct_digits_power(
	unsigned long long base,
	unsigned short exp)
{
	mpz_t n;
	mpz_init(n);
	mpz_ui_pow_ui(n, base, exp);

	unsigned long long total = 0;
	size_t count = mpz_sizeinbase(n, 8);  // Number of octal digits

	// Pure bit operations - extract last 3 bits, shift right 3
	for (size_t i = 0; i < count; i++)
	{
		unsigned int digit = mpz_get_ui(n) & 07;  // Last octal digit (0-7)
		total += digit;
		mpz_fdiv_q_2exp(n, n, 3);  // Shift right 3 bits (divide by 8)
	}

	mpz_clear(n);
	return total;
}

// Summing with this method seems WAY faster (*3)
unsigned long long sum_dec_digits_power(
	unsigned long long base,
	unsigned short exp)
{
	mpz_t num;
	mpz_init(num);
	mpz_ui_pow_ui(num, base, exp);

	char* str = mpz_get_str(nullptr, 10, num);
	size_t len = strlen(str);

	__m256i sum_vec = _mm256_setzero_si256();  // 16-bit lanes
	uint64_t sum = 0;

	size_t i = 0;
	for (; i + 16 <= len; i += 16) {
		__m128i chunk = _mm_loadu_si128((__m128i*)&str[i]);
		__m128i digits = _mm_sub_epi8(chunk, _mm_set1_epi8('0'));

		// Unpack to 16-bit and accumulate
		__m256i wide = _mm256_cvtepu8_epi16(digits);
		sum_vec = _mm256_add_epi16(sum_vec, wide);
	}

	// Extract from 16-bit vector
	uint16_t temp[16];
	_mm256_storeu_si256((__m256i*)temp, sum_vec);
	for (int j = 0; j < 16; j++)
		sum += temp[j];

	// Remainder
	for (; i < len; i++)
		sum += str[i] - '0';

	free(str);
	return sum;
}

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


		unsigned long long result = sum_digits_power(a, n);
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
