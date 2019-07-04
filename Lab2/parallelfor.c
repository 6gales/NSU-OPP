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
	#pragma omp parallel for
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

void sub(double *rv, double *lv)
{
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++)
	{
		rv[i] -= lv[i];
	}
}

void subScalar(double *rv, double *lv, double scalar)
{
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++)
	{
		rv[i] -= scalar * lv[i];
	}
}

double *multiply(double *A, double *b)
{
	double *res = malloc(sizeof(double) * N);
	
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++)
	{
		res[i] = 0.0;
	}

	#pragma omp parallel for
	for (size_t i = 0; i < N; i++)
	{
		for (size_t j = 0; j < N; j++)
		{
			res[i] += A[i * N + j] * b[j];
		}
	}
	return res;
}

double length2(double *matrix)
{
	double sum = 0.0;
	#pragma omp parallel for reduction(+:sum)
	for (size_t i = 0; i < N; i++)
	{
		sum += matrix[i] * matrix[i];
	}
	return sqrt(sum);
}

double *intermediateCalculations(double *A, double *x, double *b)
{
	double *tmp = multiply(A, x);
	sub(tmp, b);
	return tmp;
}

int main(int argc, char **argv)
{
	omp_set_num_threads((argc > 1 ? atoi(argv[1]) : 8));
	double *A = malloc(sizeof(double) * (N * N + N + N + N)),
		*b = A + N * N,
		*x = b + N,
		*tmp = x + N;

	fillA(A);
	
	#pragma omp parallel for
	for (size_t i = 0; i < N; i++)
	{
		b[i] = N + 1;
		x[i] = 0;
	}

	time_start();	
	const double bLen = length2(b);
	
	while (1)
	{
		double sum = 0.0;
		#pragma omp parallel for reduction(+:sum)
		for (size_t i = 0; i < N; i++)
		{
			tmp[i] = 0.0;
			double nestedSum = 0.0;
			#pragma omp parallel for reduction(+:sum)
			for (size_t j = 0; j < N; j++)
			{
				nestedSum += A[i * N + j] * x[j];
			}
			tmp[i] = nestedSum - b[i];
			sum += tmp[i] * tmp[i];
			tmp[i] = x[i] - t * tmp[i];
		}

		if (sqrt(sum) / bLen < eps)
			break;

		double *ptr = x;
		x = tmp;
		tmp = ptr;
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
	return 0;
}