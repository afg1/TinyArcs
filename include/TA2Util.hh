#ifndef TA2UTIL_H
#define TA2UTIL_H 1

#include <utility>
#include <vector>
#include <string>


#include "ThreeVector.hh"
#include "Magnets.hh"


namespace tau
{
    struct TAargs
    {
        std::vector<magnet*> sequence;
        ThreeVector limits;
        std::string phasespace;
        long double step;
        int nsteps;
    };
    
    class Job
    {
        public:
            Job(ThreeVector x0i, ThreeVector v0i);
            std::vector<std::pair<ThreeVector, ThreeVector> >* Run();
        private:
            ThreeVector x0;
            ThreeVector v0;
    };
    
    class Queue
    {
        public:
            Queue();
            void AddJob(Job);
            Job GetJob();
            size_t GetSize(){return jobs.size();}
        private:
            std::vector<Job> jobs;
    };
    
    void* producer(void* arg);
    void* consumer(void* arg);
        
        

    std::pair<ThreeVector, ThreeVector> TrackingStep(TAargs args);
    
    void RunTracking(TAargs args, std::vector<std::pair<ThreeVector, ThreeVector> >* res);
    
    void WriteTrackingData(std::vector<std::pair<ThreeVector, ThreeVector> >& data, std::string outloc);
}
#endif
