#ifndef TAPHASESPACE_HH
#define TAPHASESPACE_HH

#include <iostream>
#include <string>
#include <utility>

class TAPhasespace
{
    public:
        TAPhasespace(char* );
        ~TAPhasespace();
        std::pair<ThreeVector, ThreeVector> GetParticle();
    private:
        size_t positon;
        size_t len;
        std::string fname;
    
    
};

#endif
