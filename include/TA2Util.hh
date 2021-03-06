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
        bool lastStep;
        bool firstStep;
        int partNo;
        
    } ThreadArgs;
    
    void* producer(void* arg);
    void* consumer(void* arg);
    
    void TrackingStep(ThreadArgs* arg);
    
    void RunTracking(TA2ConfigParser* conf, std::vector<std::vector<std::pair<ThreeVector, ThreeVector> >* >& reslist);
    
    void WriteTrackingData(std::vector<magnet*> magnets);
    
    
    long double gammaFromE(long double T);
    long double betaFromE(long double T);
    long double gammaFromV(long double V);
    long double betaFromV(long double V);
    ThreeVector GenerateBmap(ThreeVector x, std::vector<magnet*> magnets, int);
    long double GetEndA(ThreeVector x, std::vector<magnet*> magnets);
    long double GetStartA(ThreeVector x, std::vector<magnet*> magnets);
    ThreeVector GetNormal(ThreeVector x, std::vector<magnet*> magnets);
    ThreeVector GetPlanePoint(ThreeVector x, std::vector<magnet*> magnets);
    bool GetInMagnet(ThreeVector x, std::vector<magnet*> magnets);
    bool EligibleForNR(ThreeVector x, std::vector<magnet*> magnets);
    
    void GenerateFieldMap(std::vector<magnet*> magnets, ThreeVector limits, std::string outloc, long double granularity, int ncores);
    
    void GenerateElemMap(std::vector<magnet*> magnets, ThreeVector limits, std::string outloc, int nbins, int ncores);
    
};
#endif
