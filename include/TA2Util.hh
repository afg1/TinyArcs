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
        TAPhasespace* phasespace;
        long double step;
        int nsteps;
    };
    
    void* producer(void* arg);
    void* consumer(void* arg);
        
        

    std::pair<ThreeVector, ThreeVector> TrackingStep(TAargs args);
    
    void RunTracking(TAargs args, std::vector<std::pair<ThreeVector, ThreeVector> >* res);
    
    void WriteTrackingData(std::vector<std::pair<ThreeVector, ThreeVector> >& data, std::string outloc);
}
#endif
