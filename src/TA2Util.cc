#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <cstdlib>
#include <string.h>
#include <stdio.h>


#include <pthread.h>

#include "QueueJumper.hh"

#include "TA2Util.hh"
#include "Magnets.hh"

pthread_mutex_t region_mutex = PTHREAD_MUTEX_INITIALIZER;


typedef struct
{
    tau::Queue* queue;
    tau::TAargs args;
} ProdArgs;

tau::Job::Job(ThreeVector x0i, ThreeVector v0i)
{
    x0 = x0i;
    v0 = v0i;
}

std::vector<std::pair<ThreeVector, ThreeVector> >* tau::Job::Run()
{
    std::vector<std::pair<ThreeVector, ThreeVector> >* res;
    
    return res;
}

tau::Queue::Queue()
{}

void tau::Queue::AddJob(Job j)
{
    jobs.push_back(j);
}

tau::Job tau::Queue::GetJob()
{
    tau::Job rjob = jobs[0];
    jobs.erase(jobs.begin());
    return rjob;
}

void * tau::producer(void* arg)
{
    ProdArgs* PA = (ProdArgs*) arg;
    std::cout << PA->args.phasespace.c_str() << std::endl;
// This function reads the phase space file and adds jobs to the queue
    long double tx(0), ty(0), tz(0);
    std::ifstream phase(PA->args.phasespace.c_str());
    if(phase.is_open())
    {
        phase.seekg(0, std::ios::end);
        size_t len = phase.tellg();
        phase.seekg(0, std::ios::beg);// Open and get the length of the file
        
        char* buffer = new char[6*sizeof(long double)];// buffer for a set of threads worth of particles!
        int offset(0);
        std::cout << len << std::endl;
        while(phase.tellg() != len)
        {
            phase.read(reinterpret_cast<char*>(buffer),6*sizeof(long double));
            offset = 0;
            memcpy(&tx, buffer+offset, sizeof(long double));
            offset+= sizeof(long double);
            memcpy(&ty, buffer+offset, sizeof(long double));
            offset+= sizeof(long double);
            memcpy(&tz, buffer+offset, sizeof(long double));
            offset+= sizeof(long double);
            ThreeVector x0(tx,ty,tz);
            memcpy(&tx, buffer+offset, sizeof(long double));
            offset+= sizeof(long double);
            memcpy(&ty, buffer+offset, sizeof(long double));
            offset+= sizeof(long double);
            memcpy(&tz, buffer+offset, sizeof(long double));
            offset+= sizeof(long double);
            ThreeVector v0(tx,ty,tz);
            tau::Job newJob(x0, v0);
            pthread_mutex_lock(&region_mutex);
            PA->queue->AddJob(newJob);
            pthread_mutex_unlock(&region_mutex);
            std::cout << PA->queue->GetSize() << std::endl;
        }
    }
    else
    {
        std::cerr << "Failed to open phase space!" << std::endl;
    }
}

void *tau::consumer(void* arg)
{
    tau::Queue* queue = (Queue*)arg;
    std::cout << "Consumer thread" << std::endl;
    pthread_mutex_lock(&region_mutex);
    if(queue->GetSize() != 0)
    {
        std::cout << queue->GetSize() << std::endl;
        queue->GetJob();
        std::cout << queue->GetSize() << std::endl;
    }
    pthread_mutex_unlock(&region_mutex);
}




void tau::RunTracking(tau::TAargs args, std::vector<std::pair<ThreeVector, ThreeVector> >* res)
{
    tau::Queue queue;// Initialise an empty queue to use...
    ProdArgs PA;
    PA.args = args;
    PA.queue = &queue;
    void* arg = &PA;
    
    pthread_t producer_thread;
    pthread_t* consumer_threads;
    int n(8);
    consumer_threads=(pthread_t *)malloc(n*sizeof(*consumer_threads));
    
    void *producer(void* arg);
    void *consumer(void* arg);
    
    pthread_create(&producer_thread, NULL, producer, arg);
    arg = &queue;
    for(int i=0; i<n; ++i)
    {
        pthread_create(&consumer_threads[i], NULL, consumer, arg);
    }
    for(int i=0; i<n; ++i)
    {
        pthread_join(consumer_threads[i], NULL);
    }
    
    
}


std::pair<ThreeVector, ThreeVector> TrackingStep(tau::TAargs args, std::pair<ThreeVector, ThreeVector> x0v0)
{
    ThreeVector dummy(0,0,0);
    std::pair<ThreeVector, ThreeVector> dummy_pair = std::make_pair(dummy, dummy);
    return dummy_pair;
}