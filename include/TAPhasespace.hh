#ifndef TAPHASESPACE_HH
#define TAPHASESPACE_HH

#include <iostream>
#include <string>
#include <utility>
#include <cstdlib>
#include <cstring>

#include "ThreeVector.hh"

class TAPhasespace
{
    public:
        TAPhasespace(){}
        TAPhasespace(std::string);
        ~TAPhasespace();
        std::pair<ThreeVector, ThreeVector> GetParticle();
        int GetNpart();
    private:
        size_t position;
        size_t len;
        std::string fname;
    
    
};

#endif
