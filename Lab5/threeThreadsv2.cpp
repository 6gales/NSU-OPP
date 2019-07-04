#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <string>
#include <mpi.h>
#include <pthread.h>

typedef struct ProcessInfo
{
	int rank,
		size;
} ProcessInfo;

FILE *out;

int *processStatus,
	*taskList,
	taskListSize,
	currentSize,
	busyProcesses,
	tasksLeft,
	hasExtra,
	extraTask,
	iter,
	correctOrder,
	doneWork;

double result;

pthread_mutex_t mutex;
pthread_cond_t taskDone,
	readyToGet,
	backToWork;

void handleFreeRank(int freeRank)
{
	if (!processStatus[freeRank])
	{
		busyProcesses--;
		processStatus[freeRank] = 1;
	}
	std::cout << "now it is __________" << busyProcesses << std::endl;
}

int asyncNotClosed(MPI_Request *req, int flag)
{
	if (flag != 2)
	{
		MPI_Status stat;
		MPI_Request_get_status(*req, &flag, &stat);
		if (!flag)
		{
			MPI_Cancel(req);
		}
		else return 1;
	}
	return 0;
}

void notBusyProcessMessage(int &nbflag, MPI_Request *req, int &nbrank, ProcessInfo *pinfo)
{
	MPI_Status stat;
	if (nbflag)
		MPI_Irecv(&nbrank, 1, MPI_INT, MPI_ANY_SOURCE, 3003 + iter, MPI_COMM_WORLD, req);

	MPI_Request_get_status(*req, &nbflag, &stat);

	if (nbflag)
	{
		nbflag = 2;
		if (pinfo->rank != nbrank)
		{
			MPI_Request nreq;
			handleFreeRank(nbrank);
			MPI_Isend(&nbrank, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size,
							3003 + iter, MPI_COMM_WORLD, &nreq);
		}
	}
}

void *getExtra(void *_pinfo)
{
	ProcessInfo *pinfo = (ProcessInfo*)_pinfo;
	std::cout << "#" << pinfo->rank << " getter in function\n";
	pthread_mutex_lock(&mutex);
	correctOrder--;
	pthread_cond_signal(&backToWork);
	std::cout << "#" << pinfo->rank << " getter locked\n";
	pthread_cond_wait(&readyToGet, &mutex);
	std::cout << "#" << pinfo->rank << " getter waited\n";

	MPI_Request sreq, freq, treq, rreq, nreq;
	MPI_Status fstat, tstat;

	int freeFlag = 1,
		taskFlag = 1,
		nbflag = 1,
		freeRank,
		notBusy,
		sendRank = (pinfo->rank + 1) % pinfo->size,
		rankToSend,
		localExtra,
		isRreq = 0;

	do
	{
		if (taskFlag)
		{
			MPI_Isend(&(pinfo->rank), 1, MPI_INT, sendRank,
				3344 + iter, MPI_COMM_WORLD, &sreq);
			std::cout << "#" << pinfo->rank << " getter is sending freeRank" << std::endl;
			MPI_Irecv(&localExtra, 1, MPI_INT, MPI_ANY_SOURCE, 3274 + iter, MPI_COMM_WORLD, &treq);
		}

		if (freeFlag)
			MPI_Irecv(&freeRank, 1, MPI_INT, MPI_ANY_SOURCE, 3344 + iter, MPI_COMM_WORLD, &freq);

//		std::cout << "#" << pinfo->rank << " getter launched async operations and waiting\n";
		pthread_cond_signal(&backToWork);
		pthread_cond_wait(&taskDone, &mutex);
			
		MPI_Request_get_status(freq, &freeFlag, &fstat);
		MPI_Request_get_status(treq, &taskFlag, &tstat);
		
		if (freeFlag)
		{
			freeFlag = 2;
			std::cout << "#" << pinfo->rank << " getter is about to handle " << freeRank << std::endl;
			handleFreeRank(freeRank);
			rankToSend = freeRank;
			MPI_Isend(&rankToSend, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size,
								3344 + iter, MPI_COMM_WORLD, &rreq);
			isRreq = 1;
		}

		if (taskFlag)
		{
			taskFlag = 2;
			hasExtra = 1;
			sendRank = tstat.MPI_SOURCE;
			for (int i = pinfo->rank; i != sendRank; i = (i + 1) % pinfo->size)
				handleFreeRank(i);
			extraTask = localExtra;
			std::cout << "#" << pinfo->rank << " getter is got job: " << localExtra << std::endl;
		}

		notBusyProcessMessage(nbflag, &nreq, notBusy, pinfo);		

	} while (busyProcesses);
	
	if (asyncNotClosed(&sreq, taskFlag));
	if (isRreq)
		asyncNotClosed(&rreq, 1);
	if (asyncNotClosed(&freq, freeFlag));
	if (asyncNotClosed(&treq, taskFlag))
	{
		taskFlag = 2;
		hasExtra = 1;
		extraTask = localExtra;
		std::cout << "#" << pinfo->rank << " getter is got job: " << localExtra << std::endl;
	}
	pthread_cond_signal(&backToWork);
	//	MPI_Isend(&rankToSend, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size,
	//							3344 + iter, MPI_COMM_WORLD, &rreq);
	
	correctOrder++;
	pthread_mutex_unlock(&mutex);
}

void *giveExtra(void *_pinfo)
{
	MPI_Request req, nreq, rreq, nrreq;
	MPI_Status stat, nstat;
	int freeRank,
		notBusy,
		flag = 1,
		nbflag = 1;
	ProcessInfo *pinfo = (ProcessInfo*)_pinfo;
	std::cout << "#" << pinfo->rank << " giver in function\n";
	pthread_mutex_lock(&mutex);
	std::cout << "#" << pinfo->rank << " giver locked\n";
	correctOrder--;

	while (tasksLeft)
	{
		if (flag)
			MPI_Irecv(&freeRank, 1, MPI_INT, MPI_ANY_SOURCE, 3344 + iter, MPI_COMM_WORLD, &req);
//		std::cout << "#" << pinfo->rank << " giver is trying to recv freeRank\n";

		
//		std::cout << "#" << pinfo->rank << " giver is trying to recv freeRank\n";

		std::cout << "#" << pinfo->rank << " giver is waiting\n";
		pthread_cond_signal(&backToWork);
		pthread_cond_wait(&taskDone, &mutex);		
		std::cout << "#" << pinfo->rank << " giver is waited\n";
			
		MPI_Request_get_status(req, &flag, &stat);
		

		if (flag)
		{
			flag = 2;
			std::cout << "#" << pinfo->rank << " giver is about to handle #" << freeRank << std::endl;
			handleFreeRank(freeRank);
			if (tasksLeft)
			{
				MPI_Send(&taskList[--currentSize], 1, MPI_INT, freeRank, 3274 + iter, MPI_COMM_WORLD);
				fprintf(out, "#%d giver sent extra %d to #%d and decreased size to %d\n",
					pinfo->rank, taskList[currentSize], freeRank, currentSize);
				MPI_Isend(&freeRank, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size,
							3003 + iter, MPI_COMM_WORLD, &rreq);
			}
			else MPI_Isend(&freeRank, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size,
							3344 + iter, MPI_COMM_WORLD, &rreq);
		}
		
		notBusyProcessMessage(nbflag, &nreq, notBusy, pinfo);
			 
	}

	if (asyncNotClosed(&req, flag))
	{
		MPI_Isend(&freeRank, 1, MPI_INT, (pinfo->rank + 1) % pinfo->size, 3344 + iter, MPI_COMM_WORLD, &rreq);
	}

	pthread_cond_signal(&backToWork);
	pthread_mutex_unlock(&mutex);
}

void *completeTask(void *_pinfo)
{
	ProcessInfo *pinfo = (ProcessInfo*)_pinfo;
	std::cout << "#" << pinfo->rank << " performer in function\n";
	pthread_mutex_lock(&mutex);
	while (correctOrder)
	{
		std::cout << "#" << pinfo->rank << " performer in correct start loop\n";
		pthread_cond_wait(&backToWork, &mutex);
	}
	std::cout << "#" << pinfo->rank << " performer locked\n";
	for (int taskNum = 0; taskNum < currentSize; taskNum++)
	{
		fprintf(out, "#%d performer's task %d\n", pinfo->rank, taskNum);
		std::cout << "#" << pinfo->rank << " performer's task " << taskNum << std::endl;
		doneWork += taskList[taskNum];
		for (int i = 0; i < taskList[taskNum]; i++)
		{
			result += sqrt(i);
		}

		if (taskNum + 1 == currentSize) break;
		
		std::cout << "#" << pinfo->rank << " performer signalling\n";
		pthread_cond_signal(&taskDone);
		std::cout << "#" << pinfo->rank << " performer is waiting\n";
		pthread_cond_wait(&backToWork, &mutex);
		std::cout << "#" << pinfo->rank << " performer is back to work\n";
	}

	tasksLeft = 0;
	handleFreeRank(pinfo->rank);
	std::cout << "#" << pinfo->rank << " performer signaled for no tasks\n";

	pthread_cond_signal(&taskDone);
	std::cout << "#" << pinfo->rank << " performer is waiting for back control\n";
	pthread_cond_wait(&backToWork, &mutex);
	
	std::cout << "#" << pinfo->rank << " performer signaled for ready to get\n";
	pthread_cond_signal(&readyToGet);
	pthread_cond_wait(&backToWork, &mutex);
	std::cout << "#" << pinfo->rank << " performer is back to work\n";
	do
	{
		pthread_cond_signal(&taskDone);	
		pthread_cond_wait(&backToWork, &mutex);
		if (hasExtra)
		{
			fprintf(out, "#%d performer is helping with extra: %d\n", pinfo->rank, extraTask);
			std::cout << "#" << pinfo->rank << " performer is helping\n";
			doneWork += extraTask;
			for (int i = 0; i < extraTask; i++)
			{
				result += sqrt(i);
			}
			hasExtra = 0;
		}

	} while (busyProcesses);

	pthread_mutex_unlock(&mutex);
	if (!correctOrder)
		pthread_cond_signal(&taskDone);
}

int main(int argc, char **argv)
{
	int rank, size, level;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	std::string filename = "out" + std::to_string(rank) + ".txt";
	out = fopen(filename.c_str(), "w");

	taskListSize = 10;// + rank * 10;
	taskList = new int[taskListSize];
	processStatus = new int[size];

	int sumTask = 0;
	ProcessInfo pinfo;
	pinfo.rank = rank;
	pinfo.size = size;

	pthread_mutex_init(&mutex, NULL);
	pthread_cond_init(&taskDone, NULL);
	pthread_cond_init(&backToWork, NULL);
	pthread_cond_init(&readyToGet, NULL);

	for (int it = 0; it < 20; it++)
	{
		for (int task = 0; task < taskListSize; task++)
		{
			taskList[task] = 100 * abs(rank - (it % size)) + 100;
			printf("#%d iter %d task %d filled with %d\n", rank, it, task, taskList[task]);
		}
		for (int i = 0; i < size; i++)
		{
			processStatus[i] = 0;
		}
		busyProcesses = size;
		iter = it;
		currentSize = taskListSize;
		result = 0.0;
		tasksLeft = 1;
		hasExtra = 0;
		correctOrder = 2;

		pthread_attr_t getterAttrs, giverAttrs, performerAttrs;
		pthread_t taskGetter, taskGiver, taskPerformer;

		pthread_attr_init(&getterAttrs);
		pthread_attr_setdetachstate(&getterAttrs, PTHREAD_CREATE_JOINABLE);		
		
		pthread_attr_init(&giverAttrs);
		pthread_attr_setdetachstate(&giverAttrs, PTHREAD_CREATE_JOINABLE);
		
		pthread_attr_init(&performerAttrs);
		pthread_attr_setdetachstate(&performerAttrs, PTHREAD_CREATE_JOINABLE);

		pthread_create(&taskGetter, &getterAttrs, getExtra, &pinfo);
		pthread_create(&taskGiver, &giverAttrs, giveExtra, &pinfo);
		pthread_create(&taskPerformer, &performerAttrs, completeTask, &pinfo);

		pthread_join(taskGiver, NULL);
		pthread_join(taskGetter, NULL);
		pthread_join(taskPerformer, NULL);

		pthread_attr_destroy(&getterAttrs);
		pthread_attr_destroy(&giverAttrs);
		pthread_attr_destroy(&performerAttrs);

		fprintf(out, "#%d is completed wave %d with result %f and done %d iterations\n", rank, it, result, doneWork);
		printf("#%d is completed wave %d with result %f and done %d iterations\n", rank, it, result, doneWork);
		sumTask += doneWork;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	fprintf(out, "#%d all iterations %d, %lf\n", rank, sumTask, (double)sumTask / 20);
	fprintf(stdin, "#%d all iterations %d, %lf\n", rank, sumTask, (double)sumTask / 20);

	pthread_mutex_destroy(&mutex);
	pthread_cond_destroy(&readyToGet);
	pthread_cond_destroy(&taskDone);
	pthread_cond_destroy(&backToWork);

	delete[] taskList;
	MPI_Finalize();
	return 0;
}