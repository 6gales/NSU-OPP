#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <pthread.h>

typedef struct ProcessInfo
{
	int rank,
		size;
} ProcessInfo;

int *taskList;
int taskListSize;
int currentSize;
double result;
int finished;
int tasksLeft;

pthread_mutex_t mutex;
pthread_cond_t smallIt,
	taskGiving;

// void *getExtraTask(void *_pinfo)
// {
// 	ProcessInfo *pinfo = (ProcessInfo*)_pinfo;

// 	while (!finished)
// 	{
// 		pthread_cond_wait(&finished, &mutex);

// 		MPI_Request *sreq = new MPI_Request[size], rreq;
// 		for (int i = 0; i < pinfo->size; i++)
// 		{
// 			MPI_Isend(&(pinfo->rank), 1, MPI_INT, i, 0, MPI_COMM_WORLD, &sreq[i]);
// 		}
// 		MPI_Status stat;
// 		MPI_Irecv(&extraTask, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &rreq);
		

// 	}
// }

char giversTurn = 0,
	performersTurn = 0;


void *giveExtra(void *_pinfo)
{
	ProcessInfo *pinfo = (ProcessInfo*)_pinfo;
	std::cout << "#" << pinfo->rank << " in GE\n";
	while (!finished)
	{
		MPI_Request req;
		int freeRank;
		MPI_Irecv(&freeRank, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &req);
		int flag;

		do
		{
			std::cout << "#" << pinfo->rank << ", GE: waiting for mutex\n";
			while (!giversTurn)
			{
				std::cout << "#" << pinfo->rank << ", GE: waiting...\n";
				pthread_cond_wait(&smallIt, &mutex);
			}
			std::cout << "#" << pinfo->rank << ", GE: sheet?>>\n";
			pthread_mutex_lock(&mutex);
			std::cout << "#" << pinfo->rank << ", GE: got the power\n";
			MPI_Status stat;
			MPI_Request_get_status(req, &flag, &stat);
			if (flag)
			{
				std::cout << "#" << pinfo->rank << ", GE: got flag, freeRank: " << freeRank << std::endl;
				if (!(finished = pinfo->rank == freeRank))
				{
					if (tasksLeft)
					{
						MPI_Send(&taskList[--currentSize], 1, MPI_INT, freeRank, 1, MPI_COMM_WORLD);
					}
					else
					{
						MPI_Request rreq;
						MPI_Isend(&freeRank, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size, 0, MPI_COMM_WORLD, &rreq);
					}
				}	
			}
			giversTurn ^= 1;

			pthread_mutex_unlock(&mutex);
			pthread_cond_signal(&taskGiving);
			performersTurn ^= 1;
		} while (flag == 0);
	}
}

void *completeTask(void *_pinfo)
{
	ProcessInfo *pinfo = (ProcessInfo*)_pinfo;
	std::cout << "#" << pinfo->rank << " in CT\n";
	pthread_mutex_lock(&mutex);
		std::cout << "#" << pinfo->rank << ", CT: locked mutex\n";
	for (int taskNum = 0; taskNum < currentSize; taskNum++)
	{
			std::cout << "#" << pinfo->rank << ", CT: iteration: " << taskNum << std::endl;
		for (int i = 0; i < taskList[taskNum]; i++)
		{
			result += sqrt(i);
		}
		std::cout << "#" << pinfo->rank << ", CT: giving mutex: " << std::endl;

		performersTurn ^= 1;
		pthread_mutex_unlock(&mutex);
		pthread_cond_signal(&smallIt);
		giversTurn ^= 1;

		while (!performersTurn)
			pthread_cond_wait(&taskGiving, &mutex);
		pthread_mutex_lock(&mutex);
	}
	tasksLeft = 0;
	std::cout << "#" << pinfo->rank << ", CT: finished to iterate " << std::endl;
	while (!finished)
	{
		MPI_Request sreq, rreq;
		MPI_Isend(&(pinfo->rank), 1, MPI_INT, (pinfo->rank + 1) % pinfo->size, 0, MPI_COMM_WORLD, &sreq);

		int extraTask;
		MPI_Irecv(&extraTask, 1, MPI_INT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &rreq);
		int flag;

		do
		{
			MPI_Status stat;
			MPI_Request_get_status(rreq, &flag, &stat);

			if (flag)
			{
				for (int i = 0; i < extraTask; i++)
				{
					result += sqrt(i);
				}
			}
//			pthread_mutex_unlock(&mutex);
			pthread_cond_signal(&smallIt);
			pthread_cond_wait(&taskGiving, &mutex);
			pthread_mutex_lock(&mutex);
		} while (flag == 0 && !finished);
	}
	pthread_mutex_unlock(&mutex);
}

int main(int argc, char **argv)
{
	int rank, size, level;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	taskListSize = 10;
	taskList = new int[taskListSize];

	ProcessInfo pinfo;
	pinfo.rank = rank;
	pinfo.size = size;

	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&smallIt, NULL);
	pthread_cond_init(&taskGiving, NULL);

	for (int it = 0; it < 20; it++)
	{
		for (int task = 0; task < taskListSize; task++)
		{
			taskList[task] = abs(rank - (it % size)) + 10;
			printf("#%d iter %d task %d filled with %d\n", rank, it, task, taskList[task]);
		}
		currentSize = taskListSize;
		result = 0.0;
		finished = 0;
		tasksLeft = 1;

		pthread_attr_t tgAttrs, tpAttrs;
		pthread_t taskGiver, taskPerformer;
		
		pthread_attr_init(&tgAttrs);
		pthread_attr_setdetachstate(&tgAttrs, PTHREAD_CREATE_JOINABLE);
		pthread_attr_init(&tpAttrs);
		pthread_attr_setdetachstate(&tpAttrs, PTHREAD_CREATE_JOINABLE);

		pthread_create(&taskGiver, &tgAttrs, giveExtra, &pinfo);
		pthread_create(&taskPerformer, &tpAttrs, completeTask, &pinfo);

		pthread_attr_destroy(&tpAttrs);
		pthread_attr_destroy(&tpAttrs);

		pthread_join(taskGiver, NULL);
		pthread_join(taskPerformer, NULL);

		printf("#%d is completed iteration %d with result %f\n", rank, it, result);
	}

	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&smallIt);
	pthread_cond_destroy(&taskGiving);

	delete[] taskList;
	MPI_Finalize();
	return 0;
}