#include "TAPhasespace.hh"
#include "ThreeVector.hh"

TAPhasespace::TAPhasespace(char* filename)
{
    std::ifstream phase(filename, ios::binary);
    if(phase.is_open())
    {
        phase.seekg(0, std::ios::end);
        len = phase.tellg();
        position = 0;
        phase.close();
        fname = filename;
    
    }
    else
    {
        std::cerr << "Failed to open phase space!" << std::endl;
    }
}

TAPhasespace::~TAPhasespace()
{
}

std::pair<ThreeVector, ThreeVector> GetParticle()
{
    std::ifstream phase(fname.c_str());
    if(phase.is_open())
    {
        char* buffer = new char[6*sizeof(long double)];// buffer for a set of threads worth of particles!
        int offset(0);
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
            
            position = phase.tellg();
            return std::make_pair(x0, v0);
            
        }
        else
        {
            ThreeVector v0f;
            ThreeVector x0f;
            std::cerr << "Failed to open phase space!" << std::endl;
            return std::make_pair(x0f, v0f);
        }
    }
    std::cerr << "Something went wrong... Returning garbage" << std::endl;
    ThreeVector v0f;
    ThreeVector x0f;
    return std::make_pair(x0f, v0f);
}