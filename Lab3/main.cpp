#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DIMS 2
#define SIZE_A 100
#define SIZE_B 200
#define N 300

void fill(double *matrix, size_t sizeX, size_t sizeY, double value)
{
	for (size_t i = 0; i < sizeX; i++)
	{
		for (size_t j = 0; j < sizeY; j++)
		{
			matrix[sizeY * i + j] = value;
		}
	}
}

void transpose(double **matrix, size_t sizeX, size_t sizeY)
{
	double *tmp = new double[sizeX * sizeY];
	
	for (size_t i = 0; i < sizeX; i++)
	{
		for (size_t j = 0; j < sizeY; j++)
		{
			tmp[j * sizeX + i] = (*matrix)[sizeY * i + j];
		}
	}
	delete[] *matrix;
	*matrix = tmp;
}

void fillDims(int *dims, int size)
{
	int a = (int)sqrt(size);
	
	if (size % a == 0)
	{
		dims[0] = size / a;
		dims[1] = a;
	}
	else
	{
		dims[0] = size;
		dims[1] = 1;
	}
}

void initComms(MPI_Comm cartesianComm, MPI_Comm *commLines, MPI_Comm *commColumns)
{
	int rdimsLines[] = { 0, 1 },
		rdimsColumns[] = { 1, 0 };

	MPI_Cart_sub(cartesianComm, rdimsLines, commLines);
	MPI_Cart_sub(cartesianComm, rdimsColumns, commColumns);
}

void initCommRankSize(MPI_Comm comm, int *rank, int *size)
{
	MPI_Comm_size(comm, size);
	MPI_Comm_rank(comm, rank);
}

int main(int argc, char **argv)
{
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm cartesianComm;

	int dims[DIMS],
		periods[] = { 0, 0 };

	fillDims(dims, size);

	MPI_Cart_create(MPI_COMM_WORLD, DIMS, dims, periods, 0, &cartesianComm);

	int coords[DIMS];
	MPI_Cart_coords(cartesianComm, rank, DIMS, coords);
	std::cout << "#" << rank << " coords: " << coords[0] << " " << coords[1] << std::endl;

	double *A = nullptr,
		*B = nullptr;
	if (rank == 0)
	{
		A = new double[N * SIZE_A];
		B = new double[N * SIZE_B];

		fill(A, SIZE_A, N, 1.0);
		fill(B, N, SIZE_B, 1.0);
		transpose(&B, N, SIZE_B);
	}

	MPI_Comm commLines,
		commColumns;
	initComms(cartesianComm, &commLines, &commColumns);

	int lineRank, lineSize,
		columnRank, columnSize;
	initCommRankSize(commLines, &lineRank, &lineSize);
	initCommRankSize(commColumns, &columnRank, &columnSize);

	std::cout << "#" << rank << ": " << lineRank << " " << columnRank << std::endl;

	int lineCoef = SIZE_B / lineSize,
		columnCoef = SIZE_A / columnSize,
		columnNum = N * lineCoef,
		lineNum = N * columnCoef;

	std::cout << "#" << rank << " sizes: " << lineSize << " " << columnSize << std::endl;
	std::cout << "#" << rank << " coefs: " << lineCoef << " " << columnCoef << std::endl;

	double *linesA = new double[lineNum + columnNum],
		*columnsB = linesA + lineNum;

	if (columnRank == 0)
	{
		MPI_Scatter(B, columnNum, MPI_DOUBLE, columnsB, columnNum, MPI_DOUBLE, 0, commLines);
	}

	if (lineRank == 0)
	{
		MPI_Scatter(A, lineNum, MPI_DOUBLE, linesA, lineNum, MPI_DOUBLE, 0, commColumns);
	}

	MPI_Bcast(columnsB, columnNum, MPI_DOUBLE, 0, commColumns);
	MPI_Bcast(linesA, lineNum, MPI_DOUBLE, 0, commLines);

	for (int i = 0; i < columnNum; i++)
	{
		std::cout << "b#" << rank << ": " << columnsB[i] << " ";
	}
	std::cout << std::endl;

	for (int i = 0; i < lineNum; i++)
	{
		std::cout << "a#" << rank << ": " << linesA[i] << " ";
	}
	std::cout << std::endl;

	double *result = new double[lineCoef * columnCoef];
	for (int i = 0; i < columnCoef; i++)
	{
		for (int j = 0; j < lineCoef; j++)
		{
			result[i * lineCoef + j] = 0.0;
			for (int k = 0; k < N; k++)
			{
				result[i * lineCoef + j] += linesA[i * N + k] * columnsB[j * N + k];
			}
		}
	}

	for (int i = 0; i < columnCoef; i++)
	{
		for (int j = 0; j < lineCoef; j++)
		{
			std::cout << "res#" << rank << ": " << result[lineCoef * i + j] << " ";
		}
		std::cout << std::endl;
	}

	double *partC = nullptr;
	if (lineRank == 0)
		partC = new double[columnCoef * SIZE_B];

	for (int i = 0; i < columnCoef; i++)
	{
		MPI_Gather(result + i * lineCoef, lineCoef, MPI_DOUBLE,
			partC + i * SIZE_B, lineCoef, MPI_DOUBLE, 0, commLines);
	}

	if (lineRank == 0)
	{
		for (int i = 0; i < columnCoef; i++)
		{
			for (int j = 0; j < SIZE_B; j++)
			{
				std::cout << "#" << rank << ": " << partC[SIZE_B * i + j] << " ";
			}
			std::cout << std::endl;
		}
	}

	double *C = nullptr;
	if (rank == 0)
		C = new double[SIZE_A * SIZE_B];

	if (lineRank == 0)
	{
		MPI_Gather(partC, columnCoef * SIZE_B, MPI_DOUBLE,
			C, columnCoef * SIZE_B, MPI_DOUBLE, 0, commColumns);
		delete[] partC;
	}

	if (rank == 0)
	{
		for (int i = 0; i < SIZE_A; i++)
		{
			for (int j = 0; j < SIZE_B; j++)
			{
				std::cout << C[SIZE_B * i + j] << " ";
			}
			std::cout << std::endl;
		}
		delete[] A;
		delete[] B;
		delete[] C;
	}

	delete[] linesA;
	delete[] result;

	MPI_Finalize();
	return 0;
}