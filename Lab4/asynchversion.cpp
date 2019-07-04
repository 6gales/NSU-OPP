#include <iostream>
#include <functional>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define DIMS 3
#define Z 2 //plane
#define Y 0 //column
#define X 1 //line

#define jn 64
#define in 64
#define kn 64

//#define F(k, z, y, x) partF[((k * pinfo.zCoef * pinfo.yCoef * pinfo.xCoef) + (z * pinfo.yCoef * pinfo.xCoef) + (y * pinfo.xCoef) + x)]

#define index(k, z, y, x, zSize, ySize, xSize) \
(((k) * (zSize) * (ySize) * (xSize)) + ((z) * (ySize) * (xSize)) + ((y) * (xSize)) + (x))

#define defInd(k, z, y, x) index((k), (z), (y), (x), pinfo.zCoef, pinfo.yCoef, pinfo.xCoef)
#define defIndIn(k, z, y, x) index((k), (z), (y), (x), zCoef, yCoef, xCoef)
#define tableInd(z, y, x) index(0, (z), (y), (x), planeSize, columnSize, lineSize)

//inline int index(int k, int z, int y, int x, int zSize, int ySize, int xSize)
//{
//	return ((k * zSize * ySize * xSize) + (z * ySize * xSize) + (y * xSize) + x);
//}

#define a 1

#define Fres(x, y, z) ((x) + (y) + (z))
#define Ro(x, y, z) (-a * ((x) + (y) + (z)))

#define LOOP_BODY pinfo.refreshBorders(i, j, k, partF, prevIt); \
					partF[defInd(currIt, i, j, k)] = \
						((partF[defInd(prevIt, i + 1, j, k)] \
							+ partF[defInd(prevIt, i - 1, j, k)]) / owz \
						+ (partF[defInd(prevIt, i, j + 1, k)] \
							+ partF[defInd(prevIt, i, j - 1, k)]) / owy \
						+ (partF[defInd(prevIt, i, j, k + 1)] \
							+ partF[defInd(prevIt, i, j, k - 1)]) / owx \
						- Ro((pinfo.lineRank * pinfo.xmul + k) * pinfo.hx, \
							(pinfo.columnRank * pinfo.ymul + j) * pinfo.hy, \
							(pinfo.planeRank * pinfo.zmul + i) * pinfo.hz)) / c; \
					pinfo.fillSendBuff(i, j, k, partF[defInd(currIt, i, j, k)]); \
					if (fabs(partF[defInd(currIt, i, j, k)] \
						- partF[defInd(prevIt, i, j, k)]) > e) \
					{ flag = 1; } \
					else if ((differ = fabs(partF[defInd(currIt, i, j, k)] \
						- Fres((pinfo.lineRank * pinfo.xmul + k) * pinfo.hx, \
							(pinfo.columnRank * pinfo.ymul + j) * pinfo.hy, \
							(pinfo.planeRank * pinfo.zmul + i) * pinfo.hz))) > max) \
					{ max = differ; }

class ProcessInfo
{
	void fillDims(int *dims, int size, int axis = 0)
	{
		int _sqrt = (int)sqrt(size);
		int _qbrt = (int)pow(size, 1.0 / 3.0);
		if (_qbrt > 1 && size % _qbrt == 0)
		{
			int subsize = size / _qbrt,
				_subsqrt = (int)sqrt(subsize);
			if (subsize % _subsqrt == 0)
			{
				dims[axis] = subsize / _subsqrt;
				dims[(axis + 1) % DIMS] = _subsqrt;
				dims[(axis + 2) % DIMS] = _qbrt;
			}
			else
			{
				dims[axis] = size / _qbrt;
				dims[(axis + 1) % DIMS] = _qbrt;
				dims[(axis + 2) % DIMS] = 1;
			}
		}
		else if (size % _sqrt == 0)
		{
			dims[axis] = size / _sqrt;
			dims[(axis + 1) % DIMS] = _sqrt;
			dims[(axis + 2) % DIMS] = 1;
		}
		else
		{
			dims[axis] = size;
			dims[(axis + 1) % DIMS] = 1;
			dims[(axis + 2) % DIMS] = 1;
		}
	}

	void initSubComms()
	{
		int rdimsLines[] = { 0, 1, 0 },
			rdimsColumns[] = { 1, 0, 0 },
			rdimsPlanes[] = { 0, 0, 1 };

		MPI_Cart_sub(commCube, rdimsLines, &commLines);
		MPI_Cart_sub(commCube, rdimsColumns, &commColumns);
		MPI_Cart_sub(commCube, rdimsPlanes, &commPlanes);
	}

	void initCommRankSize(MPI_Comm comm, int *rank, int *size)
	{
		MPI_Comm_size(comm, size);
		MPI_Comm_rank(comm, rank);
	}

	void initComms()
	{//0, 0 //1, 2 //2, 1
		int dims[DIMS],
			periods[] = { 0, 0, 0 };

		fillDims(dims, size, 0);

		MPI_Cart_create(MPI_COMM_WORLD, DIMS, dims, periods, 0, &commCube);

		initSubComms();

		initCommRankSize(commLines, &lineRank, &lineSize);
		initCommRankSize(commColumns, &columnRank, &columnSize);
		initCommRankSize(commPlanes, &planeRank, &planeSize);

		std::cout << "#" << rank << " ranks, line: " << lineRank << ", column: " << columnRank << ", plane: " << planeRank << std::endl;
		std::cout << "#" << rank << " sizes, line: " << lineSize << ", column: " << columnSize << ", plane: " << planeSize << std::endl;
	}

	void putBorders()
	{
		xCoef = kn / lineSize + 2;
		yCoef = jn / columnSize + 2;
		zCoef = in / planeSize + 2;
		std::cout << "#" << rank << " coefs, z: " << zCoef << ", y: " << yCoef << ", x: " << xCoef << std::endl;
		xtill = xCoef - 1;
		ytill = yCoef - 1;
		ztill = zCoef - 1;

		xmul = xCoef - 2;
		ymul = yCoef - 2;
		zmul = zCoef - 2;

		std::cout << "#" << rank << " muls: " << xmul << " " << ymul << " " << zmul << std::endl;

		effectiveX = xtill - 1;
		effectiveY = ytill - 1;
		effectiveZ = ztill - 1;
	}

	void fillRankTable()
	{
		MPI_Cart_coords(commCube, rank, DIMS, coords);
		std::cout << "#" << rank << " coords: " << coords[0] << " " << coords[1] << " " << coords[2] << std::endl;

		rankTable = new int[lineSize * columnSize * planeSize];
		
		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0)
			rankTable[tableInd(coords[Z], coords[Y], coords[X])] = 0;
		int otherCoords[DIMS];
		
		for (int i = 1; i < size; i++)
		{
			if (rank == i)
				MPI_Send(coords, DIMS, MPI_INT, 0, i, MPI_COMM_WORLD);
			if (rank == 0)
			{
				MPI_Recv(otherCoords, DIMS, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				rankTable[tableInd(otherCoords[Z], otherCoords[Y], otherCoords[X])] = i;
			}
		}
		
		MPI_Bcast(rankTable, planeSize * columnSize * lineSize, MPI_INT, 0, MPI_COMM_WORLD);
	}

	void putFlags()
	{
		firstZ = planeRank == 0;
		lastZ = planeRank == planeSize - 1;

		firstY = columnRank == 0;
		lastY = columnRank == columnSize - 1;

		firstX = lineRank == 0;
		lastX = lineRank == lineSize - 1;
	}

	void setVectorsValue(double *va, double *vb, double *vc, double value, int size)
	{
		for (int i = 0; i < size; i++)
		{
			va[i] = value;
			vb[i] = value;
			vc[i] = value;
		}
	}

	void nullBuffs(int position)
	{
		sendBuf[position] = nullptr;
		recvBuf[position] = nullptr;
	}
	
	void allocateBuffs(int position, int size, int neighborRank, int i)
	{
		sideSizes[i] = size;

		sendBuf[position] = new double[size];
		recvBuf[position] = new double[size];

		setVectorsValue(sendBuf[position], sendBuf[buffSize + position],
			recvBuf[position], 0.0, size);

		neighbors[i] = neighborRank;
		effectiveSides[i] = position;
		std::cout << "Allocated #" << rank << "'s buffer " << position
			<< ", neighbor is #" << neighborRank << std::endl;
	}

	void initBuffs()
	{
		effectiveSidesNum = !firstZ + !lastZ + !firstY + !lastY + !firstX + !lastX;
		std::cout << "#" << rank << "'s effective sides num " << effectiveSidesNum << std::endl;
		effectiveSides = new int[effectiveSidesNum];
		sideSizes = new int[effectiveSidesNum];
		neighbors = new int[effectiveSidesNum];
		requests = new MPI_Request[effectiveSidesNum * 2];
		statuses = new MPI_Status[effectiveSidesNum];
		int i = 0;

		if (firstZ)
		{
			nullBuffs(POS_Z);
		}
		else
		{
			allocateBuffs(POS_Z, effectiveX * effectiveY,
				rankTable[tableInd(coords[Z] - 1, coords[Y], coords[X])], i);
			i++;
		}
		if (lastZ)
		{
			nullBuffs(NEG_Z);
		}
		else
		{
			allocateBuffs(NEG_Z, effectiveX * effectiveY,
				rankTable[tableInd(coords[Z] + 1, coords[Y], coords[X])], i);
			i++;
		}

		if (firstY)
		{
			nullBuffs(POS_Y);
		}
		else
		{
			allocateBuffs(POS_Y, effectiveX * effectiveZ,
				rankTable[tableInd(coords[Z], coords[Y] - 1, coords[X])], i);
			i++;
		}
		if (lastY)
		{
			nullBuffs(NEG_Y);
		}
		else
		{
			allocateBuffs(NEG_Y, effectiveX * effectiveZ,
				rankTable[tableInd(coords[Z], coords[Y] + 1, coords[X])], i);
			i++;
		}

		if (firstX)
		{
			nullBuffs(POS_X);
		}
		else
		{
			allocateBuffs(POS_X, effectiveY * effectiveZ,
				rankTable[tableInd(coords[Z], coords[Y], coords[X] - 1)], i);
			i++;
		}
		if (lastX)
		{
			nullBuffs(NEG_X);
		}
		else
		{
			allocateBuffs(NEG_X, effectiveY * effectiveZ,
				rankTable[tableInd(coords[Z], coords[Y], coords[X] + 1)], i);
			i++;
		}
	}

	void fillFunctions()
	{
		functions[POS_X] = [this](int z, int y, int x, double value)
		{ sendBuf[POS_X][(z - 1) * effectiveX + y - 1] = value; };

		functions[NEG_X] = [this](int z, int y, int x, double value)
		{ sendBuf[NEG_X][(z - 1) * effectiveX + y - 1] = value; };

		functions[POS_Y] = [this](int z, int y, int x, double value)
		{ sendBuf[POS_Y][(z - 1) * effectiveX + x - 1] = value; };

		functions[NEG_Y] = [this](int z, int y, int x, double value)
		{ sendBuf[NEG_Y][(z - 1) * effectiveX + x - 1] = value; };

		functions[POS_Z] = [this](int z, int y, int x, double value)
		{ sendBuf[POS_Z][(y - 1) * effectiveX + x - 1] = value; };

		functions[NEG_Z] = [this](int z, int y, int x, double value)
		{ sendBuf[NEG_Z][(y - 1) * effectiveX + x - 1] = value; };
	}

public:
	enum Direction { POS_X, NEG_X, POS_Y, NEG_Y, POS_Z, NEG_Z };

	int effectiveSidesNum;

	int coords[DIMS],
		*rankTable,
		*sideSizes,
		*neighbors,
		*effectiveSides;
	
	static constexpr int buffSize = DIMS * 2;
	double *sendBuf[buffSize],
		*recvBuf[buffSize];

	MPI_Comm commCube,
		commLines,
		commColumns,
		commPlanes;

	MPI_Request *requests;
	MPI_Status *statuses;

	int rank, size,
		lineRank, lineSize,
		columnRank, columnSize,
		planeRank, planeSize,
		
		xCoef, yCoef, zCoef,
		xfrom = 1, xtill,
		yfrom = 1, ytill,
		zfrom = 1, ztill,

		xmul, ymul, zmul,
		
		effectiveX, effectiveY, effectiveZ; //last_ = _Coef - 2

	double hx = 2.0 / (kn + 1),
		hy = 2.0 / (jn + 1),
		hz = 2.0 / (in + 1);
	
	bool firstX, firstY, firstZ,
		lastX, lastY, lastZ;

	std::function<void (int, int, int, double)> functions[buffSize];
	int toCall[buffSize];
	int toCallCount = 0;

	ProcessInfo(int _rank, int _size) : rank(_rank), size(_size)
	{
		initComms();

		putFlags();

		std::cout << std::boolalpha << "#" << rank << " fx: " << firstX
			<< ", lx: " << lastX << ", fy: " << firstY << ", ly: " << lastY
			<< ", fz: " << firstZ << ", lz: " << lastZ << std::endl;
		putBorders();

		fillRankTable();

		initBuffs();

		fillFunctions();
	}

	~ProcessInfo()
	{
		delete[] rankTable;
		for (int i = 0; i < buffSize; i++)
		{
			delete[] recvBuf[i];
			delete[] sendBuf[i];
		}

		delete[] sideSizes;
		delete[] effectiveSides;
		delete[] requests;
		delete[] statuses;
		delete[] neighbors;
	}
	
	void refreshBorders(int z, int y, int x, double *F, int k)
	{
		toCallCount = 0;
		if (z == zfrom && !firstZ)
		{
			F[defIndIn(k, z - 1, y, x)] = recvBuf[POS_Z][(y - 1) * effectiveX + x - 1];
			toCall[toCallCount++] = POS_Z;
		}
		if (z == effectiveZ && !lastZ)
		{
			F[defIndIn(k, z + 1, y, x)] = recvBuf[NEG_Z][(y - 1) * effectiveX + x - 1];
			toCall[toCallCount++] = NEG_Z;
		}
		if (y == yfrom && !firstY)
		{
			F[defIndIn(k, z, y - 1, x)] = recvBuf[POS_Y][(z - 1) * effectiveX + x - 1];
			toCall[toCallCount++] = POS_Y;
		}
		if (y == effectiveY && !lastY)
		{
			F[defIndIn(k, z, y + 1, x)] = recvBuf[NEG_Y][(z - 1) * effectiveX + x - 1];
			toCall[toCallCount++] = NEG_Y;
		}
		if (x == xfrom && !firstX)
		{
			F[defIndIn(k, z, y, x - 1)] = recvBuf[POS_X][(z - 1) * effectiveX + y - 1];
			toCall[toCallCount++] = POS_X;
		}
		if (x == effectiveX && !lastX)
		{
			F[defIndIn(k, z, y, x + 1)] = recvBuf[NEG_X][(z - 1) * effectiveX + y - 1];
			toCall[toCallCount++] = NEG_X;
		}
	}

	void fillSendBuff(int z, int y, int x, double value)
	{
		while (toCallCount)
		{
			functions[toCall[--toCallCount]](z, y, x, value);
		}
	}

	void exchangeBorders()
	{
		for (int i = 0; i < effectiveSidesNum; i++)
			MPI_Isend(sendBuf[effectiveSides[i]], sideSizes[i], MPI_DOUBLE,
				neighbors[i], rank + neighbors[i], MPI_COMM_WORLD, &requests[i]);

		for (int i = 0; i < effectiveSidesNum; i++)
			MPI_Irecv(recvBuf[effectiveSides[i]], sideSizes[i], MPI_DOUBLE,
				neighbors[i], rank + neighbors[i], MPI_COMM_WORLD, &requests[effectiveSidesNum + i]);
	}

	void waitBorders()
	{
		MPI_Waitall(effectiveSidesNum, requests + effectiveSidesNum, statuses);
	}

	template <class T>
	MPI_Datatype getDatatype()
	{
		return MPI_CHAR;
	}

	template <>
	MPI_Datatype getDatatype<long>()
	{
		return MPI_LONG;
	}

	template <>
	MPI_Datatype getDatatype<double>()
	{
		return MPI_DOUBLE;
	}

	template <class T>
	void findMaxValue(T &flag)
	{
		T *flags = nullptr;
		if (columnRank == 0)
			flags = new T[columnSize];

		MPI_Gather(&flag, 1, getDatatype<T>(), flags, 1, getDatatype<T>(), 0, commColumns);

		if (columnRank == 0)
		{
			for (int i = 0; i < planeSize; i++)
			{
				if (flags[i] > flag)
					flag = flags[i];
			}

			delete[] flags;
			flags = nullptr;
			
			if (lineRank == 0)
				flags = new T[lineSize];

			MPI_Gather(&flag, 1, getDatatype<T>(), flags, 1, getDatatype<T>(), 0, commLines);

			if (lineRank == 0)
			{
				for (int i = 0; i < lineSize; i++)
				{
					if (flags[i] > flag)
						flag = flags[i];
				}
				delete[] flags;
				flags = nullptr;

				if (planeRank == 0)
					flags = new T[planeSize];

				MPI_Gather(&flag, 1, getDatatype<T>(), flags, 1, getDatatype<T>(), 0, commPlanes);

				if (planeRank == 0)
				{
					for (int i = 0; i < planeSize; i++)
					{
						if (flags[i] > flag)
							flag = flags[i];
					}
					delete[] flags;
					flags = nullptr;
				}
			}
		}
		MPI_Bcast(&flag, 1, getDatatype<T>(), 0, MPI_COMM_WORLD);
	}
};

double *initF(ProcessInfo &pinfo)
{
	double *partF = new double[2 * pinfo.xCoef * pinfo.yCoef * pinfo.zCoef];

	for (int i = 0; i < pinfo.zCoef; i++)
	{
		for (int j = 0; j < pinfo.yCoef; j++)
		{
			for (int k = 0; k < pinfo.xCoef; k++)
			{
				bool borderValue = false;
				if ((pinfo.firstZ && i == 0)
					|| (pinfo.lastZ && i == pinfo.zCoef - 1)
					|| (pinfo.firstY && j == 0)
					|| (pinfo.lastY && j == pinfo.yCoef - 1)
					|| (pinfo.firstX && k == 0)
					|| (pinfo.lastX && k == pinfo.xCoef - 1))
				{
					borderValue = true;
				}

				if (borderValue)
				{
					partF[defInd(0, i, j, k)] =
						Fres((pinfo.lineRank * pinfo.xmul + k) * pinfo.hx,
						(pinfo.columnRank * pinfo.ymul + j) * pinfo.hy,
						(pinfo.planeRank * pinfo.zmul + i) * pinfo.hz);
					partF[defInd(1, i, j, k)] = partF[defInd(0, i, j, k)];
				}
				else
				{
					partF[defInd(0, i, j, k)] = 0;
					partF[defInd(1, i, j, k)] = 0;
				}
			}
		}
	}
	return partF;
}

int main(int argc, char **argv)
{
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	ProcessInfo pinfo{ rank, size };

	double *partF = initF(pinfo);

	double owx = pow(pinfo.hx, 2),
		owy = pow(pinfo.hy, 2),
		owz = pow(pinfo.hz, 2),
		c = 2.0 / owx + 2.0 / owy + 2.0 / owz + a,
		e = 0.00001,
		differ,
		max = 0.0;

	char currIt = 0,
		prevIt = 1,
		flag = 1;

	int iteration = 0;

	do
	{
		std::cout << "#" << rank << ": iteration: " << iteration++ << std::endl;

		max = 0.0;
		flag = 0;
		currIt ^= 1;
		prevIt ^= 1;

		for (int i = pinfo.zfrom; i < pinfo.ztill; i += (pinfo.ztill - pinfo.zfrom - 1))
		{
			for (int j = pinfo.yfrom; j < pinfo.ytill; j++)
			{
				for (int k = pinfo.xfrom; k < pinfo.xtill; k++)
				{
					LOOP_BODY
				}
			}
		}

		for (int j = pinfo.yfrom; j < pinfo.ytill; j += (pinfo.ytill - pinfo.yfrom - 1))
		{
			for (int i = pinfo.zfrom; i < pinfo.ztill; i++)
			{
				for (int k = pinfo.xfrom; k < pinfo.xtill; k++)
				{
					LOOP_BODY
				}
			}
		}

		for (int k = pinfo.xfrom; k < pinfo.xtill; k += (pinfo.xtill - pinfo.xfrom - 1))
		{
			for (int i = pinfo.zfrom; i < pinfo.ztill; i++)
			{
				for (int j = pinfo.yfrom; j < pinfo.ytill; j++)
				{
					LOOP_BODY
				}
			}
		}

		pinfo.exchangeBorders();

		for (int i = pinfo.zfrom + 1; i < pinfo.ztill - 1; i++)
		{
			for (int j = pinfo.yfrom + 1; j < pinfo.ytill - 1; j++)
			{
				for (int k = pinfo.xfrom + 1; k < pinfo.xtill - 1; k++)
				{
					partF[defInd(currIt, i, j, k)] =

						((partF[defInd(prevIt, i + 1, j, k)]
							+ partF[defInd(prevIt, i - 1, j, k)]) / owz

							+ (partF[defInd(prevIt, i, j + 1, k)]
								+ partF[defInd(prevIt, i, j - 1, k)]) / owy

							+ (partF[defInd(prevIt, i, j, k + 1)]
								+ partF[defInd(prevIt, i, j, k - 1)]) / owx

							- Ro((pinfo.lineRank * pinfo.xmul + k) * pinfo.hx,
							(pinfo.columnRank * pinfo.ymul + j) * pinfo.hy,
								(pinfo.planeRank * pinfo.zmul + i) * pinfo.hz)) / c;

					if (fabs(partF[defInd(currIt, i, j, k)]
						- partF[defInd(prevIt, i, j, k)]) > e)
					{
						flag = 1;
					}
					else if ((differ = fabs(partF[defInd(currIt, i, j, k)]
						- Fres((pinfo.lineRank * pinfo.xmul + k) * pinfo.hx,
						(pinfo.columnRank * pinfo.ymul + j) * pinfo.hy,
							(pinfo.planeRank * pinfo.zmul + i) * pinfo.hz))) > max)
					{
						max = differ;
					}
				}
			}
		}

		pinfo.waitBorders();

		pinfo.findMaxValue(flag);

	} while (flag);

	pinfo.findMaxValue(differ);

	if (rank == 0)
	{
		std::cout << "Max differ: " << max << std::endl;
	}

	delete[] partF;
	MPI_Finalize();
	return 0;
}