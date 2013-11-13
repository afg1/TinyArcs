#include <iostream>
#include <fstream>
#include <cmath>
#include <utility>
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <vector>


#include <pthread.h>

#include "QueueJumper.hh"

#include "TA2Util.hh"
#include "Magnets.hh"

pthread_mutex_t region_mutex = PTHREAD_MUTEX_INITIALIZER;


typedef struct
{
    ApolloQueue* queue;
    tau::TAargs args;
} ProdArgs;

void *tau::consumer(void* arg)
{
/*
This function will do the actual traking
*/
    tau::TAargs* input = (tau::TAargs*)arg
    std::cout << "Consumer thread" << std::endl;

}




void tau::RunTracking(tau::TAargs args, std::vector<std::pair<ThreeVector, ThreeVector> >* res)
{
    ApolloQueue queue(8);

}


std::pair<ThreeVector, ThreeVector> TrackingStep(tau::TAargs args, std::pair<ThreeVector, ThreeVector> x0v0)
{
    ThreeVector dummy(0,0,0);
    std::pair<ThreeVector, ThreeVector> dummy_pair = std::make_pair(dummy, dummy);
    return dummy_pair;
}