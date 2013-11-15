#ifndef TA2UTIL_H
#define TA2UTIL_H 1

#include <utility>
#include <vector>
#include <string>


#include "ThreeVector.hh"
#include "TA2ConfigParse.hh"
#include "Magnets.hh"


namespace tau
{
    typedef struct
    {
        TA2ConfigParser* conf;
        std::vector<std::pair<ThreeVector, ThreeVector> >* res;
        std::pair<ThreeVector, ThreeVector> initPair;
        
    } ThreadArgs;
    
    void* producer(void* arg);
    void* consumer(void* arg);
    
    void TrackingStep(ThreadArgs* arg);
    
    void RunTracking(TA2ConfigParser* conf, std::vector<std::vector<std::pair<ThreeVector, ThreeVector> >* >& reslist);
    
    void WriteTrackingData(std::vector<std::pair<ThreeVector, ThreeVector> >& data, std::string outloc);
    
    
    long double gammaFromE(long double T);
    long double betaFromE(long double T);
    long double gammaFromV(long double V);
    long double betaFromV(long double V);
    ThreeVector GenerateBmap(ThreeVector x, std::vector<magnet*> magnets);
    long double GetEndA(ThreeVector x, std::vector<magnet*> magnets);
    ThreeVector GetNormal(ThreeVector x, std::vector<magnet*> magnets);
    ThreeVector GetPlanePoint(ThreeVector x, std::vector<magnet*> magnets);
    bool GetInMagnet(ThreeVector x, std::vector<magnet*> magnets);


    
}
#endif
