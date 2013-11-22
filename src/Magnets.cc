#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

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
    point -= centre;
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
    point -= centre;
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
    point -= centre;
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
    point -= centre;
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

// Magnet from map...
GeneralMagnetFromMap::GeneralMagnetFromMap(std::string namei, long double extXi, long double extYi, long double extZi, ThreeVector cent, std::string mapfile)
{
    name = namei;
    extX = extXi;
    extY = extYi;
    extZ = extZi;
    centre = cent;
    
    std::ifstream map;
    map.open(mapfile.c_str(), std::ios::binary);
    double tx(0), ty(0), tz(0);
    if(map.is_open())
    {
        map.seekg(0, std::ios::end);
        size_t len = map.tellg();
        map.seekg(0, std::ios::beg);
        char* buffer = new char[6*sizeof(double)];
        int offset(0);
        while(map.tellg() != len && map.tellg() != -1)
        {
            map.read(reinterpret_cast<char*>(buffer),6*sizeof(double));
            offset = 0;
            memcpy(&tx, buffer+offset, sizeof(double));
            offset+= sizeof(double);
            memcpy(&ty, buffer+offset, sizeof(double));
            offset+= sizeof(double);
            memcpy(&tz, buffer+offset, sizeof(double));
            offset+= sizeof(double);
            ThreeVector x(tx,ty,tz);
            
            memcpy(&tx, buffer+offset, sizeof(double));
            offset+= sizeof(double);
            memcpy(&ty, buffer+offset, sizeof(double));
            offset+= sizeof(double);
            memcpy(&tz, buffer+offset, sizeof(double));
            offset+= sizeof(double);
            ThreeVector Bi(tx,ty,tz);
            mapData.push_back(std::make_pair(x, Bi));
        }
    }
    long double testx(mapData[0].first.GetElem(0));
    long double testy(mapData[0].first.GetElem(1));
    long double testz(mapData[0].first.GetElem(2));
    dx = dy = dz = 0;
    bool xdone(false), ydone(false), zdone(false);
    minx = mapData[0].first.GetElem(0);
    miny = mapData[0].first.GetElem(1);
    minz = mapData[0].first.GetElem(2);
    maxx = mapData[0].first.GetElem(0);
    maxy = mapData[0].first.GetElem(1);
    maxz = mapData[0].first.GetElem(2);
    
    std::vector<std::pair<ThreeVector, ThreeVector> >::iterator curr;
    for(curr = mapData.begin(); curr!=mapData.end(); ++curr)
    {
        // find the voxel sizes...
        if(curr->first.GetElem(0) != testx && !xdone)
        {
            dx = curr->first.GetElem(0) - testx;
            xdone = true;
        }
        else if(curr->first.GetElem(1) != testy && !ydone)
        {
            dy = curr->first.GetElem(1) - testy;
            ydone = true;
        }
        else if(curr->first.GetElem(2) != testz && !zdone)
        {
            dz = curr->first.GetElem(2) - testz;
            zdone = true;
        }
        
        //Find the edges of the defined region - we assume its a cube
        if(curr->first.GetElem(0) < minx)
        {
            minx = curr->first.GetElem(0);
        }
        else if(curr->first.GetElem(0) > maxx)
        {
            maxx = curr->first.GetElem(0);
        }
        
        if(curr->first.GetElem(1) < miny)
        {
            miny = curr->first.GetElem(1);
        }
        else if(curr->first.GetElem(1) > maxy)
        {
            maxy = curr->first.GetElem(1);
        }
        
        if(curr->first.GetElem(2) < minz)
        {
            minz = curr->first.GetElem(2);
        }
        else if(curr->first.GetElem(2) > maxz)
        {
            maxz = curr->first.GetElem(2);
        }
        
        nx = int(ceil((maxx - minx)/dx)) + 1;
        ny = int(ceil((maxy - miny)/dy)) + 1;
        nz = int(ceil((maxz - minz)/dz)) + 1;
        
        
    }
    std::cerr << "Debug: magnet constructed" << std::endl;
}


ThreeVector GeneralMagnetFromMap::B(ThreeVector point)
{
    ThreeVector rval = TriLinearInterpolate(point);
    return rval;
}
bool GeneralMagnetFromMap::InMagnet(ThreeVector point)
{
    if(abs(point.GetElem(0)) < extX &&   abs(point.GetElem(1)) < extY  &&  abs(point.GetElem(2)) < extZ)
    {
        return true;
    }
    return false;
}
ThreeVector GeneralMagnetFromMap::GetNormal()
{
    std::cerr << "WARNING: Attempting to get exit normal from a map-element - returning (0,0,0)" << std::endl;
    ThreeVector rval;
    return rval;
}
ThreeVector GeneralMagnetFromMap::GetPlanePoint()
{
    std::cerr << "WARNING: Attempting to get plane point from a map-element - returning (0,0,0)" << std::endl;
    ThreeVector rval;
    return rval;
}
long double GeneralMagnetFromMap::GetEndA()
{
    std::cerr << "WARNING: Attempting to get end angle from a map-element - returning 0" << std::endl;
    return 0.0;
}
long double GeneralMagnetFromMap::GetStartA()
{
    std::cerr << "WARNING: Attempting to get start angle from a map-element - returning 0" << std::endl;
    return 0.0;
}


long double GeneralMagnetFromMap::LinearInterpolate(long double y0, long double y1, long double x0, long double x1, long double x)
{
    if((x1 - x0) != 0.0)
    {
        long double rval = y0 + (y1 - y0)*(x - x0)/(x1 - x0);
        return rval;
    }
    return 0.0;
}


ThreeVector GeneralMagnetFromMap::TriLinearInterpolate(ThreeVector point)
{
    if(!InMagnet(point))
    {
        ThreeVector rval(0.0,0.0,0.0);
        return rval;
    }

    std::vector<std::pair<ThreeVector, ThreeVector> > selected;
    long double bxtemp1(0), bxtemp2(0), bxtemp3(0), bxtemp4(0), bxtemp5(0), bxtemp6(0);
    long double bytemp1(0), bytemp2(0), bytemp3(0), bytemp4(0), bytemp5(0), bytemp6(0);
    long double bztemp1(0), bztemp2(0), bztemp3(0), bztemp4(0), bztemp5(0), bztemp6(0);// Yeah, a shitload of temporary variables
    long double btx(0), bty(0), btz(0);

    int ix(0), jy(0), kz(0);

    
    ix = int(ceil((point.GetElem(0) - minx)/dx));
    jy = int(ceil((point.GetElem(1) - miny)/dy));// Find the index of the lower corner of the voxel!
    kz = int(ceil((point.GetElem(2) - minz)/dz));

    if(ix < nx-1)// Must have stepped outside the region we have a field map for... return zero
    {
        selected.push_back(mapData[ix*ny*nz + jy*nz + kz]);
        selected.push_back(mapData[ix*ny*nz + (jy+1)*nz + kz]);
        selected.push_back(mapData[ix*ny*nz + jy*nz + (kz+1)]);
        selected.push_back(mapData[ix*ny*nz + (jy+1)*nz + (kz+1)]);
        
        selected.push_back(mapData[(ix+1)*ny*nz + jy*nz + kz]);
        selected.push_back(mapData[(ix+1)*ny*nz + (jy+1)*nz + kz]);
        selected.push_back(mapData[(ix+1)*ny*nz + jy*nz + (kz+1)]);
        selected.push_back(mapData[(ix+1)*ny*nz + (jy+1)*nz + (kz+1)]);// The order of these is very particular!
        
        std::vector<std::pair<ThreeVector, ThreeVector> >::iterator curr;
        std::vector<ThreeVector> Btest;
        ThreeVector zero;
        
        // Interpolate the B-vector in the x-direction...
        bxtemp1 = LinearInterpolate(selected[0].second.GetElem(0), selected[4].second.GetElem(0), selected[0].first.GetElem(0), selected[4].first.GetElem(0), point.GetElem(0));
        bxtemp2 = LinearInterpolate(selected[1].second.GetElem(0), selected[5].second.GetElem(0), selected[1].first.GetElem(0), selected[5].first.GetElem(0), point.GetElem(0));
        bxtemp3 = LinearInterpolate(selected[2].second.GetElem(0), selected[6].second.GetElem(0), selected[2].first.GetElem(0), selected[6].first.GetElem(0), point.GetElem(0));
        bxtemp4 = LinearInterpolate(selected[3].second.GetElem(0), selected[7].second.GetElem(0), selected[3].first.GetElem(0), selected[7].first.GetElem(0), point.GetElem(0));
        // interpolate Bx in x
        
        bytemp1 = LinearInterpolate(selected[0].second.GetElem(1), selected[4].second.GetElem(1), selected[0].first.GetElem(0), selected[4].first.GetElem(0), point.GetElem(0));
        bytemp2 = LinearInterpolate(selected[1].second.GetElem(1), selected[5].second.GetElem(1), selected[1].first.GetElem(0), selected[5].first.GetElem(0), point.GetElem(0));
        bytemp3 = LinearInterpolate(selected[2].second.GetElem(1), selected[6].second.GetElem(1), selected[2].first.GetElem(0), selected[6].first.GetElem(0), point.GetElem(0));
        bytemp4 = LinearInterpolate(selected[3].second.GetElem(1), selected[7].second.GetElem(1), selected[3].first.GetElem(0), selected[7].first.GetElem(0), point.GetElem(0));
        //interpolate By in x
        
        bztemp1 = LinearInterpolate(selected[0].second.GetElem(2), selected[4].second.GetElem(2), selected[0].first.GetElem(0), selected[4].first.GetElem(0), point.GetElem(0));
        bztemp2 = LinearInterpolate(selected[1].second.GetElem(2), selected[5].second.GetElem(2), selected[1].first.GetElem(0), selected[5].first.GetElem(0), point.GetElem(0));
        bztemp3 = LinearInterpolate(selected[2].second.GetElem(2), selected[6].second.GetElem(2), selected[2].first.GetElem(0), selected[6].first.GetElem(0), point.GetElem(0));
        bztemp4 = LinearInterpolate(selected[3].second.GetElem(2), selected[7].second.GetElem(2), selected[3].first.GetElem(0), selected[7].first.GetElem(0), point.GetElem(0));
        // interpolate Bz in x
        
        // Now interpolate in the y-direction...
        bxtemp5 = LinearInterpolate(bxtemp1, bxtemp2, selected[0].first.GetElem(1), selected[1].first.GetElem(1), point.GetElem(1));
        bxtemp6 = LinearInterpolate(bxtemp3, bxtemp4, selected[0].first.GetElem(1), selected[1].first.GetElem(1), point.GetElem(1));

        bytemp5 = LinearInterpolate(bytemp1, bytemp2, selected[0].first.GetElem(1), selected[1].first.GetElem(1), point.GetElem(1));
        bytemp6 = LinearInterpolate(bytemp3, bytemp4, selected[0].first.GetElem(1), selected[1].first.GetElem(1), point.GetElem(1));
        
        bztemp5 = LinearInterpolate(bztemp1, bztemp2, selected[0].first.GetElem(1), selected[1].first.GetElem(1), point.GetElem(1));
        bztemp6 = LinearInterpolate(bztemp3, bztemp4, selected[0].first.GetElem(1), selected[1].first.GetElem(1), point.GetElem(1));
        
        // Finally in the z-direction
        btx = LinearInterpolate(bxtemp5, bxtemp6, selected[0].first.GetElem(2), selected[2].first.GetElem(2), point.GetElem(2));
        bty = LinearInterpolate(bytemp5, bytemp6, selected[0].first.GetElem(2), selected[2].first.GetElem(2), point.GetElem(2));
        btz = LinearInterpolate(bztemp5, bztemp6, selected[0].first.GetElem(2), selected[2].first.GetElem(2), point.GetElem(2));
        
        ThreeVector rval(btx, bty, btz);
        return rval;
    }
    
    ThreeVector rval;
    return rval;
    
}









