#ifndef TA2CONFIGPARSE_H
#define TA2CONFIGPARSE_H 1

#include <vector>


#include "ThreeVector.hh"
#include "TAPhasespace.hh"
#include "Magnets.hh"

class TA2ConfigParser
{
    public:
        TA2ConfigParser(const char* fname);
        long double GetStep(){return step;}
        unsigned int GetNsteps(){return nsteps;}
        int GetNcores(){return ncores;}
        ThreeVector& GetLimits(){return *limits;}
        std::vector<magnet*> GetMagnets(){return magnets;}
        TAPhasespace* GetPhaseSpace(){return phasespace;}
        bool DoBmap(){return genBmap;}
        long double GetGranularity(){return granularity;}
        std::string GetBmapOutloc(){return outloc;}
        std::pair<ThreeVector, ThreeVector> GetParticle();
        int GetNpart();
    
    private:
        // Things needed to construct a magnet
        long double innerR;
        long double outerR;
        long double startA;
        long double endA;
        long double gap;
        long double alpha;
        long double betaDS;
    
        long double extX;
        long double extY;
        long double extZ;
        std::string mapfile;
        ThreeVector* centre;
        ThreeVector* B0;
        std::string name;
    
    
    
        // Global options for the algorithm
        long double step;
        int nsteps;
        int ncores;
        TAPhasespace* phasespace;
        bool suppress_output;
        ThreeVector* limits;
        ThreeVector* xi;
        ThreeVector* vi;
        std::vector<ThreeVector> configPhaseSpacexi;
        std::vector<ThreeVector> configPhaseSpacevi;
        std::vector<magnet*> magnets;
    
        // Constants and things
        bool configGood;
        bool usePhasespace;
        bool genBmap;
        long double granularity;
        long double pi;
        std::string outloc;



};


#endif
