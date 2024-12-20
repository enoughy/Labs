#include <stdio.h>
#include <math.h>
#include <time.h>

int main() {
	struct timespec start, end;
	unsigned long long N;
	long double res = 0;
	scanf("%llu", &N);
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);
	unsigned long long n = 0;
	while (n < N) {
        if (n%2==0) {
        	res += 1 / (2 * n + 1);
        }
        else {
            res += 1 / (2 * n + 1);
        }
		n++;
	}
	res *= 4;
	printf("%Lf \n", res);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	printf("Time taken: %lf sec.\n", end.tv_sec-start.tv_sec + 0.000000001*(end.tv_nsec-start.tv_nsec));
}
