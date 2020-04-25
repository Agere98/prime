#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <omp.h>

// ReSharper disable CommentTypo
// ReSharper disable StringLiteralTypo
// ReSharper disable CppClangTidyCppcoreguidelinesMacroUsage
// ReSharper disable CppCStyleCast

#define MAX_SIZE (1<<28) // 268 435 456
#define MAX_PRIME (long long)(2*MAX_SIZE+3) // 536 870 915
#define MAX_UB (MAX_PRIME*MAX_PRIME-1LL) // 288 230 379 372 937 224
#define MAX_WIDTH (2*MAX_SIZE) // 536 870 912
#define INDEX(x) ((int)((x)-3)/2)
#define INDEX_B(x, lb) ((int)((x)-(lb))/2)
char sieve[MAX_SIZE + 1];
long long primes[1 << 25];

int find_primes(long long lb, const long long ub)
{
	auto count = 0;
	const auto sq = (int)sqrt(ub);
	char mark = 1;
#pragma omp parallel
	{
		const auto id = omp_get_thread_num();
		const auto total = omp_get_num_threads();
		auto timer = omp_get_wtime();
		// Wyznaczanie liczb pierwszych mniejszych lub równych sqrt(ub)
#pragma omp for schedule(dynamic)
		for (auto i = 0; i < (sq - 1) / 2; i++)
		{
			if (sieve[i] == mark)continue;
			const auto prime = 2 * i + 3;
			if ((long long)prime * (long long)prime > sq)break;
			for (auto j = INDEX(prime * prime); j <= INDEX(sq); j += prime)
			{
				sieve[j] = mark;
			}
		}
#pragma omp single
		for (auto i = 0; i < (sq - 1) / 2; i++)
		{
			if (sieve[i] != mark)
			{
				primes[count++] = 2LL * i + 3LL;
			}
		}
		timer = omp_get_wtime() - timer;
		printf("Faza 1: %fs\n", timer);
		timer = omp_get_wtime();
		// Wyznaczanie liczb pierwszych w przedziale [lb, ub]
		mark = 2;
		if (lb < 2LL)lb = 2LL;
		const auto interval = (ub - lb) / total;
		const auto _lb = lb + id * interval;
		const auto _ub = (id + 1 == total) ? ub : _lb + interval - 1LL;
		for (auto i = 0; i < count; i++)
		{
			const auto prime = primes[i];
			auto start = _lb - _lb % prime;
			if (start % 2 == 0)start += prime;
			if (start < _lb || start == prime)start += 2 * prime;
			for (auto j = INDEX_B(start, lb); j <= INDEX_B(_ub, lb); j += prime)
			{
				sieve[j] = mark;
			}
		}
		timer = omp_get_wtime() - timer;
		printf("Faza 2: %fs\n", timer);
	}
	// Odczytywanie wyników
	count = 0;
	if (lb <= 2LL && ub >= 2LL)primes[count++] = 2LL;
	lb += 1LL - (lb % 2LL);
	for (auto i = 0; lb + i + i <= ub; i++)
	{
		if (sieve[i] != mark)
		{
			primes[count++] = lb + i + i;
		}
	}
	return count;
}

int main(const int argc, char* argv[])
{
	if (argc < 4)
	{
		printf("Parametry: M, N, num_threads");
		return 0;
	}
	const auto num_threads = strtol(argv[3], nullptr, 0);
	const auto m = strtoll(argv[1], nullptr, 0);
	const auto n = strtoll(argv[2], nullptr, 0);
	// Parametry skrajne: 288230378836066313 288230379372937224
	if (m < 1LL || n < 1LL || m > n || n > MAX_UB || (n - m + 1LL) > MAX_WIDTH)
	{
		printf("Nieprawidłowy przedział.\n");
		return 0;
	}
	omp_set_num_threads(num_threads);
	const auto start = omp_get_wtime();
	const auto p = find_primes(m, n);
	const auto stop = omp_get_wtime();
	if (argc > 4)
	{
		for (auto i = 0; i < p; i++)
		{
			printf("%lld ", primes[i]);
			if (i % 10 == 9 || i == p - 1)printf("\n");
		}
	}
	printf("W zakresie od %lld do %lld znaleziono %d liczb pierwszych.\n", m, n, p);
	printf("Czas przetwarzania: %fs.\n", stop - start);
	return 0;
}
