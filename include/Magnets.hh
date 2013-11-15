#ifndef MAGNETS_H
#define MAGNETS_H 1

#include <string>

#include "ThreeVector.hh"

class magnet
{
    public:
        virtual ThreeVector B(ThreeVector point)=0;
        virtual bool InMagnet(ThreeVector point)=0;
        virtual ThreeVector* GetNormal()=0;
        virtual ThreeVector* GetPlanePoint()=0;
        virtual long double GetEndA()=0;

};


class HardEdgedArcDipole : public magnet
{
    public:
        HardEdgedArcDipole(std::string , long double , long double, long double, long double, long double, ThreeVector*, ThreeVector*);
        ~HardEdgedArcDipole();
        virtual ThreeVector B(ThreeVector point);
        ThreeVector* GetNormal(){return exit_normal;}
        ThreeVector* GetPlanePoint(){return exit_plane_point;}
        long double GetEndA(){return endA;}
        bool InMagnet(ThreeVector point);
    
    private:
        long double innerR;
        long double outerR;
        long double midR;
        long double startA;
        long double endA;
        long double gap;
        std::string name;
        ThreeVector* B0;
        ThreeVector* centre;
        ThreeVector* exit_normal;
        ThreeVector* exit_plane_point;
    
        long double pi;
    
};

class HardEdged225Spectrometer : public magnet
{
    public:
        HardEdged225Spectrometer(std::string , long double , long double, long double, long double, long double, ThreeVector*, ThreeVector*, long double, long double);
        ~HardEdged225Spectrometer();
        long double Hr(long double rp, long double z);
        virtual ThreeVector B(ThreeVector point);
        ThreeVector* GetNormal(){return exit_normal;}
        ThreeVector* GetPlanePoint(){return exit_plane_point;}
        long double GetEndA(){return endA;}
        bool InMagnet(ThreeVector point);
    
    private:
        long double innerR;
        long double outerR;
        long double midR;
        long double startA;
        long double endA;
        long double gap;
        long double alpha;
        long double beta;
        std::string name;
        ThreeVector* B0;
        ThreeVector* centre;
    
        ThreeVector* exit_normal;
        ThreeVector* exit_plane_point;
    
        long double pi;
    
};


#endif