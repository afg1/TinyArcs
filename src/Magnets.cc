#include <cmath>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "Magnets.hh"

HardEdgedArcDipole::HardEdgedArcDipole(std::string namei, long double inner, long double outer, long double start, long double end, long double gapi, ThreeVector* Bvect, ThreeVector* cent)
{
    pi = boost::math::constants::pi<long double>();
    B0 = Bvect;
    centre = cent;
    
    innerR = inner;
    outerR = outer;
    startA = start;
    midR = innerR +(outerR - innerR)/2;
    endA = end;
    name = namei;
    gap = gapi;
    
    
     exit_normal = new ThreeVector(-sin(endA), cos(endA), 0);
     exit_plane_point = new ThreeVector(midR*cos(endA), midR*sin(endA), gap/2);
}

HardEdgedArcDipole::~HardEdgedArcDipole()
{
    delete B0;
    delete centre;
    delete exit_normal;
    delete exit_plane_point;
}

bool HardEdgedArcDipole::InMagnet(ThreeVector point)
{
    long double theta = std::atan2(point.GetElem(0), point.GetElem(1)) + pi;// By trial and error, this ends up right
    long double rp = std::sqrt(std::pow(point.GetElem(0), 2) + std::pow(point.GetElem(1), 2));
    long double z = point.GetElem(2);
    
    if((z <= gap/2 && z >= -gap/2) && ((theta >= startA && theta <= endA) || (theta <= startA && theta >= endA)) && (rp >= innerR && rp <= outerR))
    {
        return true;
    }
    else
    {
        return false;
    }
}

ThreeVector HardEdgedArcDipole::B(ThreeVector point)
{
    long double theta = std::atan2(point.GetElem(0), point.GetElem(1)) + pi;// By trial and error, this ends up right
    long double rp = std::sqrt(std::pow(point.GetElem(0), 2) + std::pow(point.GetElem(1), 2));
    long double z = point.GetElem(2);
    
    if((z <= gap/2 && z >= -gap/2) && ((theta >= startA && theta <= endA) || (theta <= startA && theta >= endA)) && (rp >= innerR && rp <= outerR))
    {
        return *B0;
    }
    else
    {
        ThreeVector rval(0.0, 0.0, 0.0);//Initialise vector to 0 as well
        return rval;
    }
}






HardEdged225Spectrometer::HardEdged225Spectrometer(std::string namei, long double inner, long double outer, long double start, long double end, long double gapi, ThreeVector* Bvect, ThreeVector* cent, long double alph, long double bet)
{
    pi = boost::math::constants::pi<long double>();
    B0 = new ThreeVector(*Bvect);
    centre = new ThreeVector(*cent);
    
    innerR = inner;
    outerR = outer;
    startA = start;
    midR = innerR +(outerR - innerR)/2;
    endA = end;
    name = namei;
    gap = gapi;
    alpha = alph;
    beta= bet;
    
    
     exit_normal = new ThreeVector(-sin(endA), cos(endA), 0);
     exit_plane_point = new ThreeVector(midR*cos(endA), midR*sin(endA), gap/2);
}

long double HardEdged225Spectrometer::Hr(long double rp, long double z)
{
    long double rval;
    long double H0 = B0->GetElem(2);
    rval = z * H0 * ( (2 * beta * (rp - midR) )/( pow(midR,2) ) - (alpha/midR));
    return rval;
}

bool HardEdged225Spectrometer::InMagnet(ThreeVector point)
{
    long double theta = std::atan2(point.GetElem(0), point.GetElem(1)) + pi;// By trial and error, this ends up right
    long double rp = std::sqrt(std::pow(point.GetElem(0), 2) + std::pow(point.GetElem(1), 2));
    long double z = point.GetElem(2);
    
    if((z <= gap/2 && z >= -gap/2) && ((theta >= startA && theta <= endA) || (theta <= startA && theta >= endA)) && (rp >= innerR && rp <= outerR))
    {
        return true;
    }
    else
    {
        return false;
    }
}

ThreeVector HardEdged225Spectrometer::B(ThreeVector point)
{
    long double theta = std::atan2(point.GetElem(0), point.GetElem(1)) + pi;// By trial and error, this ends up right
    long double rp = std::sqrt(std::pow(point.GetElem(0), 2) + std::pow(point.GetElem(1), 2));
    long double z = point.GetElem(2);
    long double H0 = B0->GetElem(2);
    
    if((z <= gap/2 && z >= -gap/2) && ((theta >= startA && theta <= endA) || (theta <= startA && theta >= endA)) && (rp >= innerR && rp <= outerR))
    {
        long double Hx = -1*cos(theta)*Hr(rp, z);
        long double Hy = sin(theta)*Hr(rp, z);
        long double Hz1 = ( (rp - midR) / midR ) * alpha * H0;
        long double Hz2 = pow( ( (rp - midR) / midR ),2) * beta * H0;
        long double Hz3 = pow( (z/midR) ,2 )* beta * H0;
        
        long double Hz = H0 * rp - Hz1 + Hz2 - Hz3;
        ThreeVector rval(Hx, Hy, Hz);
        return rval;
    }
    else
    {
        ThreeVector rval(0.0, 0.0, 0.0);//Initialise vector to 0 as well
        return rval;
    }
}
