#ifndef TA2CONFIGPARSE_H
#define TA2CONFIGPARSE_H 1

#include <vector>


#include "ThreeVector.hh"
#include "Magnets.hh"

class TA2ConfigParser
{
    public:
        TA2ConfigParser(const char* fname);
        long double GetStep(){return step;}
        unsigned int GetNsteps(){return nsteps;}
        std::string GetPhaseSpace(){return phasespace;}
        ThreeVector& GetLimits(){return *limits;}
        std::vector<magnet*> GetMagnets(){return magnets;}

    
    private:
        // Things needed to construct a magnet
        long double innerR;
        long double outerR;
        long double startA;
        long double endA;
        long double gap;
        long double alpha;
        long double betaDS;
        ThreeVector* centre;
        ThreeVector* B0;
    
        // Global options for the algorithm
        long double step;
        int nsteps;
        std::string phasespace;
        bool suppress_output;
        ThreeVector* limits;
        ThreeVector* xi;
        ThreeVector* vi;
        std::vector<magnet*> magnets;
    
        // Constants and things
        bool configGood;
        long double pi;



};


#endif