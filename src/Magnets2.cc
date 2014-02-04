#include "Magnets.hh"

#include <iostream>
#include <fstream>
#include <sstream>

DetectorMagnet::DetectorMagnet(std::string namei, std::string output, ThreeVector cent, ThreeVector ext)
{
    centre = cent;
    extent = ext;
    name = namei;
    outloc = output;
    
}

DetectorMagnet::~DetectorMagnet()
{
//    delete hits;
}

bool DetectorMagnet::InMagnet(ThreeVector point)
{
    point -= centre;
    if((point.GetElem(0) <= extent.GetElem(0) && point.GetElem(0) >= -extent.GetElem(0)) && (point.GetElem(1) <= extent.GetElem(1) && point.GetElem(1) >= -extent.GetElem(1)) && (point.GetElem(2) <= extent.GetElem(2) && point.GetElem(2) >= -extent.GetElem(2)))
    {
        return true;
    }
    else
    {
        return false;
    }
}

ThreeVector DetectorMagnet::B(ThreeVector point, int n)
{
    if(InMagnet(point))
    {
        RegisterHit(point, n);
    }
    ThreeVector rval(0.0, 0.0, 0.0);
    return rval;
}

ThreeVector DetectorMagnet::GetNormal()
{
    // Define the normal as the normal to the largest surface
    // Assume that its the y direction for now...
    
    ThreeVector rval(0.0, 1.0, 0.0);
    return rval;
}

ThreeVector DetectorMagnet::GetPlanePoint()
{
    // return a point that should be in the centre of the plane of the thing.
    ThreeVector rval(centre.GetElem(0), (centre.GetElem(1) + extent.GetElem(1)), centre.GetElem(2));
    return rval;
}

long double DetectorMagnet::GetEndA()
{
    std::cerr << "WARNING: Attempting to get end angle from a square magnet - returning 0" << std::endl;
    return 0.0;
}
long double DetectorMagnet::GetStartA()
{
    std::cerr << "WARNING: Attempting to get start angle from a square magnet - returning 0" << std::endl;
    return 0.0;
}

void DetectorMagnet::RegisterHit(ThreeVector point, int n)
{
    hits.at(n).push_back(point);
}

void DetectorMagnet::WriteHits()
{
    for(int n=0; n< npart; n++)
    {
        std::vector<ThreeVector>::iterator curr;
        std::stringstream outname;
        std::string outname_str;
        outname << outloc << "_" << name << "_p" << n;
        outname_str = outname.str();
        std::ofstream output(outname_str.c_str(), std::ios::binary);
        for(curr = hits.at(n).begin(); curr!=hits.at(n).end(); ++curr)
        {
            curr->BinaryWriteToStream(output);
        }
    }
    
}

void DetectorMagnet::SetNpart(int nparti)
{
    npart=nparti;
    std::vector<ThreeVector> temp;
    hits.assign(npart, temp);
}

