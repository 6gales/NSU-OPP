#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>

double global_res = 0;
int *task_list;
int list_size;
int current_i;
pthread_mutex_t mutex;
pthread_cond_t condition;
int size, rank;
int token;
int *array;
int count_process;
int extra_task;
int has_extra;
int iter;
int balance_flag = 0;
int weight;

void hand_free_rank(int free_rank){
    if(array[free_rank]) {
        array[free_rank] = 0;
        count_process--;
    }
}

void *communication_thread(void * ptr){

    MPI_Request request, request_token, request_exit;
    int free_rank = 0;
    int fr_to_send = 0;
    int recieved_flag = 0;
    int request_init = 1;
    int task_flag = 1;
    int got_task_flag = 0;
    int exit_flag = 0;
    int have_lock = 0;
    int dummy = 0;

    pthread_mutex_lock(&mutex);
    have_lock = 1;
    int local_token = token;
    array = malloc(sizeof(int)*size);
    for (int i = 0; i < size; ++i) {
        array[i] = 1;
    }

    MPI_Irecv(&fr_to_send, 1, MPI_INT, MPI_ANY_SOURCE, 777 + iter, MPI_COMM_WORLD, &request_exit);

    while(count_process) {

        if (!local_token && have_lock) {
            pthread_mutex_unlock(&mutex);
            have_lock = 0;
        }

        if (request_init) {
            MPI_Irecv(&free_rank, 1, MPI_INT, MPI_ANY_SOURCE, 333 + iter, MPI_COMM_WORLD,
                      &request);
                      request_init = 0;
        }

        if(task_flag && local_token){
            MPI_Isend(&rank, 1, MPI_INT, (rank+1)%size, 333 + iter, MPI_COMM_WORLD, &request_token);
            MPI_Irecv(&extra_task, 1, MPI_INT, MPI_ANY_SOURCE, 666 + iter, MPI_COMM_WORLD, &request_token);
            task_flag = 0;
        }

///////
        MPI_Request_get_status(request, &recieved_flag, MPI_STATUS_IGNORE);
        MPI_Request_get_status(request_exit, &exit_flag, MPI_STATUS_IGNORE);

        if(local_token){
            MPI_Request_get_status(request_token, &got_task_flag, MPI_STATUS_IGNORE);

        }

        if (exit_flag)
            break;

        /////
        if (!have_lock) {
            pthread_mutex_lock(&mutex);
            have_lock = 1;
        }

        local_token = token;

        if (recieved_flag)
        {
            hand_free_rank(free_rank);

            if(token){
                fr_to_send = free_rank;
                MPI_Isend(&fr_to_send, 1, MPI_INT, (rank+1)%size, 333 + iter, MPI_COMM_WORLD, &request);
            }
            else if (balance_flag) {
                MPI_Send(task_list + current_i++, 1, MPI_INT, free_rank, 666 + iter, MPI_COMM_WORLD);
                printf("I send the task with weight %d, count of working process %d\n", task_list[current_i-1], count_process);
            }//изменила уменьшение размера на увеличение
            else MPI_Send(&dummy, 1, MPI_INT, free_rank, 666 + iter, MPI_COMM_WORLD);
            request_init = 1;

        }

        if (token && got_task_flag)
        {
            has_extra = 1;
            task_flag = 1;
            printf("I recieved task with weight %d, count of working process %d\n", extra_task, count_process);
            pthread_mutex_unlock(&mutex);
            have_lock = 0;
        }

        pthread_cond_signal(&condition);
    }
    printf("exited in comm\n");

    if (!have_lock)
        pthread_mutex_lock(&mutex);

    count_process = 0;
    for (int i = 0; i < size; ++i) {
        if (rank != i)
            MPI_Isend(&rank, 1, MPI_INT, i, 777 + iter, MPI_COMM_WORLD, &request);
    }
    printf("before signal in comm\n");
    free(array);
    pthread_cond_signal(&condition);
    pthread_mutex_unlock(&mutex);
    printf("unlocked in comm\n");

}

void *worker_thread(void *ptr){
    pthread_mutex_lock(&mutex);
    current_i = 0;
    while(current_i < list_size){//может быть false если коммуникатор отошлёт последнюю задачу

        weight+=task_list[current_i];
        //как только выполнили задачу увеличиваем счётчик
        for (int j = 0; j < task_list[current_i]; ++j){
            global_res += sqrt(j);
        }
        current_i++;

        if(current_i == list_size) break;//если выполнили последнее задание, то не посылаем сигнал, а изменяем токен
        pthread_cond_wait(&condition, &mutex);
    }
    token = 1;
    hand_free_rank(rank);
    while(count_process){
        if(has_extra){
            for (int i = 0; i < extra_task; ++i)
                global_res += sqrt(i);
            has_extra = 0;
            weight += extra_task;
            printf("I did the task with weight %d, count of working process %d\n", extra_task, count_process);
        }
        pthread_cond_wait(&condition, &mutex);
    }
    printf("worker out of while\n");
    pthread_mutex_unlock(&mutex);

}

int main(int argc, char **argv)
{
    balance_flag = 1;
    if (argc > 1) {
        balance_flag = argv[1][0] - '0';//отключение балансировки
    }
    srand(time(NULL));
    int level;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &level);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int result = 0;
    int result_in_array = 0;
    int result_in_func = 0;
    list_size = 10;
    task_list = malloc(sizeof(int)*list_size);
    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&condition, NULL);
    pthread_attr_t comm_attr, work_attr;


    pthread_attr_init(&comm_attr);
    pthread_attr_init(&work_attr);

    pthread_attr_setdetachstate(&comm_attr, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setdetachstate(&work_attr, PTHREAD_CREATE_JOINABLE);

    for (iter = 0; iter < list_size; ++iter) {
        token = 0;
        count_process = size;
        has_extra = 0;

        result_in_array = 0;
        result_in_func = 0;
        weight = 0;
        result = 0;

        for (int i = 0; i < list_size; ++i) {
            task_list[i] = rand()%10;
            result+=task_list[i];
        }

        pthread_t worker, communicator;

        pthread_create(&worker, &work_attr, worker_thread, NULL);
        pthread_create(&communicator, &comm_attr, communication_thread, NULL);

        pthread_join(worker, NULL);
        printf("rank #%d worker joined\n", rank);
        pthread_join(communicator, NULL);
        printf("rank #%d, iter; %d, w %d\n", rank, iter, weight);
        printf("All joined\n");

        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&result, &result_in_array, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&weight, &result_in_func, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        if (rank == 0)
            printf("%d == %d", result_in_array, result_in_func);
    }

    pthread_attr_destroy(&comm_attr);
    pthread_attr_destroy(&work_attr);
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&condition);
    MPI_Finalize();
    return 0;
}