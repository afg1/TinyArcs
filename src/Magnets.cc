#include <cmath>
#include <iostream>

#include <boost/math/constants/constants.hpp>

#include "Magnets.hh"

long double magnet::gammaFromV(long double V)
{
    long double b = betaFromV(V);
    long double g = 1.0/(sqrtl(1.0 - (b*b)));
    return g;
}

long double magnet::betaFromV(long double V)
{
    long double c(299792458);
    long double b = V/c;
    return b;
}



HardEdgedArcDipole::HardEdgedArcDipole(std::string namei, long double inner, long double outer, long double start, long double end, long double gapi, ThreeVector Bvect, ThreeVector cent)
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
    
     exit_normal.SetElem(0, -sin(endA));
     exit_normal.SetElem(1, cos(endA));
     exit_normal.SetElem(2, 0.0);
     exit_plane_point.SetElem(0,midR*cos(endA));
     exit_plane_point.SetElem(1, midR*sin(endA));
     exit_plane_point.SetElem(2, gap/2);
}

HardEdgedArcDipole::~HardEdgedArcDipole()
{

}



bool HardEdgedArcDipole::InMagnet(ThreeVector point)
{
    long double theta = std::atan2(point.GetElem(1), point.GetElem(0));
    if(theta < 0)
    {
        theta += 2*pi;
    }
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
    long double theta = std::atan2(point.GetElem(1), point.GetElem(0));
    if(theta < 0)
    {
        theta += 2*pi;
    }
    long double rp = std::sqrt(std::pow(point.GetElem(0), 2) + std::pow(point.GetElem(1), 2));
    long double z = point.GetElem(2);
    
    if((z <= gap/2 && z >= -gap/2) && ((theta >= startA && theta <= endA) || (theta <= startA && theta >= endA)) && (rp >= innerR && rp <= outerR))
    {
        return B0;
    }
    else
    {
        ThreeVector rval(0.0, 0.0, 0.0);//Initialise vector to 0 as well
        return rval;
    }
}






HardEdged225Spectrometer::HardEdged225Spectrometer(std::string namei, long double inner, long double outer, long double start, long double end, long double gapi, ThreeVector Bvect, ThreeVector cent, long double alph, long double bet)
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
    alpha = alph;
    beta = bet;
    
    
     exit_normal.SetElem(0, -sin(endA));
     exit_normal.SetElem(1, cos(endA));
     exit_normal.SetElem(2, 0.0);
     exit_plane_point.SetElem(0,midR*cos(endA));
     exit_plane_point.SetElem(1, midR*sin(endA));
     exit_plane_point.SetElem(2, gap/2);
}

long double HardEdged225Spectrometer::Hr(long double rp, long double z)
{
    long double rval;
    long double H0 = B0.GetElem(2);
    rval = z * H0 * ( (2 * beta * (rp - midR) )/(midR*midR) - (alpha/midR) );
    return rval;
}

bool HardEdged225Spectrometer::InMagnet(ThreeVector point)
{
    long double theta = std::atan2(point.GetElem(1), point.GetElem(0));// By trial and error, this ends up right
    if(theta < 0)
    {
        theta += 2*pi;
    }
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
    long double theta = std::atan2(point.GetElem(1), point.GetElem(0));// By trial and error, this ends up right
    if(theta < 0)
    {
        theta += 2*pi;
    }
    long double rp = std::sqrt(std::pow(point.GetElem(0), 2) + std::pow(point.GetElem(1), 2));
    long double z = point.GetElem(2);
    long double H0 = B0.GetElem(2);
    
    if((z <= gap/2 && z >= -gap/2) && ((theta >= startA && theta <= endA) || (theta <= startA && theta >= endA)) && (rp >= innerR && rp <= outerR))
    {
        long double Hx = cos(theta)*Hr(rp, z);
        long double Hy = sin(theta)*Hr(rp, z);
        
        long double Hz1 = ( (rp - midR) / midR ) * alpha * H0;
        
        long double Hz2 = ( ( (rp - midR) / midR )*( (rp - midR) / midR ) ) * beta * H0;
        
        long double Hz3 = ( (z/midR)*(z/midR) )* beta * H0;
        
        long double Hz = H0 - Hz1 + Hz2 - Hz3;
        ThreeVector rval(Hx, Hy, Hz);
        return rval;
    }
    else
    {
        ThreeVector rval(0.0, 0.0, 0.0);//Initialise vector to 0 as well
        return rval;
    }
}
