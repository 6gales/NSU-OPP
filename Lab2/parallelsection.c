#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#define N 12000
#define t 10e-6
#define eps 10e-9

struct timeval tv1, tv2, dtv;

void time_start() { gettimeofday(&tv1, NULL); } 

size_t time_stop()
{
	gettimeofday(&tv2, NULL); 
	dtv.tv_sec = tv2.tv_sec - tv1.tv_sec; 
	dtv.tv_usec = tv2.tv_usec - tv1.tv_usec; 
	if (dtv.tv_usec < 0) { dtv.tv_sec--; dtv.tv_usec += 1000000; } 
	return dtv.tv_sec * 1000 + dtv.tv_usec / 1000; 
}

void fillA(double *A)
{
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			A[i * N + j] = 1.0;
			if (i == j)
				A[i * N + j]++;
		}
	}
}

void fillB(double *b)
{
	for (size_t i = 0; i < N; i++)
	{
		b[i] = N + 1;
	}
}

void fillX(double *x)
{
	for (size_t i = 0; i < N; i++)
	{
		x[i] = 0;
	}
}

int main(int argc, char **argv)
{
	omp_set_num_threads((argc > 1 ? atoi(argv[1]) : 8));

	int n = omp_get_max_threads(),
		coef = N / n,
		rest = N % n;

	double *A = malloc(sizeof(double) * (N * N + N + N + N + n)),
		*b = A + N * N,
		*x = b + N,
		*tmp = x + N,
		*data = tmp + N;

	fillA(A);
	fillB(b);
	fillX(x);
	
	int *display = malloc(sizeof(int) * (n + 1));
	double bLen, tmpLen;
	
	display[0] = 0;
	for (size_t i = 0; i < n; i++)
	{
		display[i + 1] = display[i] + coef + (i < rest);
	}

	time_start();
	#pragma omp parallel
	{
		int rank = omp_get_thread_num();

		data[rank] = 0;
		for (int i = display[rank]; i < display[rank + 1]; i++)
		{
			data[rank] += b[i] * b[i];
		}
		#pragma omp barrier

		#pragma omp single
		{
			double sum = 0.0;
			for (int i = 0; i < n; i++)
			{
				sum += data[i];
			}
			bLen = sqrt(sum);
		}

		while (1)
		{
			data[rank] = 0.0;
			for (int i = display[rank]; i < display[rank + 1]; i++)
			{
				tmp[i] = 0.0;
				for (size_t j = 0; j < N; j++)
				{
					tmp[i] += A[i * N + j] * x[j];
				}
				tmp[i] -= b[i];
				data[rank] += tmp[i] * tmp[i];
				tmp[i] = x[i] - t * tmp[i];
			}
			#pragma omp barrier
	
			#pragma omp single
			{
				double sum = 0.0;
				for (int i = 0; i < n; i++)
				{
					sum += data[i];
				}
				tmpLen = sqrt(sum);
			}
			
			if (tmpLen / bLen < eps)
				break;

			#pragma omp single
			{
				double *ptr = x;
				x = tmp;
				tmp = ptr;
			}
		}
	}
	
	printf("Elapsed in %u\n", time_stop());
	
	double accur = 0.0;
	for (size_t i = 0; i < N; i++)
	{
		if (1.0 - x[i] > accur)
			accur = 1 - x[i];
	}
	printf("Accurasity: %f\n", accur / N);

	free(A);
	free(display);
	return 0;
}