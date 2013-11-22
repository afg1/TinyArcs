#ifndef MAGNETS_H
#define MAGNETS_H 1

#include <string>

#include "ThreeVector.hh"

class magnet
{
    public:
        virtual ThreeVector B(ThreeVector point)=0;
        virtual bool InMagnet(ThreeVector point)=0;
        virtual ThreeVector GetNormal()=0;
        virtual ThreeVector GetPlanePoint()=0;
        virtual long double GetEndA()=0;
        virtual long double GetStartA()=0;
        // Note these two aren't virtual, all inheriting classes get them!
        long double betaFromV(long double V);
        long double gammaFromV(long double V);

};


class HardEdgedArcDipole : public magnet
{
    public:
        HardEdgedArcDipole(std::string , long double , long double, long double, long double, long double, ThreeVector, ThreeVector);
        ~HardEdgedArcDipole();
        virtual ThreeVector B(ThreeVector point);
        ThreeVector GetNormal(){return exit_normal;}
        ThreeVector GetPlanePoint(){return exit_plane_point;}
        long double GetEndA(){return endA;}
        long double GetStartA(){return startA;}
        bool InMagnet(ThreeVector point);
    
    private:
        long double innerR;
        long double outerR;
        long double midR;
        long double R0;
        long double startA;
        long double endA;
        long double gap;
        std::string name;
        ThreeVector B0;
        ThreeVector centre;
        ThreeVector exit_normal;
        ThreeVector exit_plane_point;
    
        long double pi;
    
};

class HardEdged225Spectrometer : public magnet
{
    public:
        HardEdged225Spectrometer(std::string , long double , long double, long double, long double, long double, ThreeVector, ThreeVector, long double, long double);
        ~HardEdged225Spectrometer();
        long double Hr(long double rp, long double z);
        virtual ThreeVector B(ThreeVector point);
        ThreeVector GetNormal(){return exit_normal;}
        ThreeVector GetPlanePoint(){return exit_plane_point;}
        long double GetEndA(){return endA;}
        long double GetStartA(){return startA;}
        bool InMagnet(ThreeVector point);
    
    private:
        long double innerR;
        long double outerR;
        long double midR;
        long double R0;
        long double startA;
        long double endA;
        long double gap;
        long double alpha;
        long double beta;
        std::string name;
        ThreeVector B0;
        ThreeVector centre;
    
        ThreeVector exit_normal;
        ThreeVector exit_plane_point;
    
        long double pi;
    
};

class HardEdgedGeneralDipole : public magnet
{
    public:
        HardEdgedGeneralDipole(std::string);
        ~HardEdgedGeneralDipole();
        ThreeVector B(ThreeVector point);
    
        bool InMagnet(ThreeVector point);
        ThreeVector GetNormal();
        ThreeVector GetPlanePoint();
        long double GetEndA();
        long double GetStartA();
    
    private:
        long double innerR;
        long double outerR;
        long double startA;
        long double endA;
        long double entryPoleAngle;
        long double exitPoleAngle;
        long double gap;
        std::string name;
        ThreeVector B0;
        ThreeVector centre;
        ThreeVector exit_normal;
        ThreeVector exit_plane_point;
    
        long double pi;
};
    
class GeneralMagnetFromMap : public magnet
{
    public:
        GeneralMagnetFromMap(std::string, long double, long double, long double, ThreeVector, std::string);
        ~GeneralMagnetFromMap(){}
    
        ThreeVector B(ThreeVector point);
        bool InMagnet(ThreeVector point);
        ThreeVector GetNormal();
        ThreeVector GetPlanePoint();
        long double GetEndA();
        long double GetStartA();
    
        long double LinearInterpolate(long double y0, long double y1, long double x0, long double x1, long double x);
        ThreeVector TriLinearInterpolate(ThreeVector point);
    
    private:
        long double extX;
        long double extY;
        long double extZ;
        std::string name;
        ThreeVector centre;
    
    
    private:
        static bool IsZero(ThreeVector i)
        {
            ThreeVector zero;
            return i == zero;
        }
        std::vector<std::pair<ThreeVector, ThreeVector> > mapData;
        long double minx;
        long double maxx;
        long double miny;
        long double maxy;
        long double minz;
        long double maxz;
        long double dx;
        long double dy;
        long double dz;
        int nx;
        int ny;
        int nz;
};






#endif