// Slightly reworked config parser for TA2

// Standard Libraries
#include <string>
#include <vector>
#include <fstream>

// External Libraries
#include <boost/algorithm/string.hpp>// Needed for the string splitting
#include <boost/math/constants/constants.hpp>

// Parts of this code
#include "TA2ConfigParse.hh"
#include "ThreeVector.hh"
#include "Magnets.hh"

class ThreeVector;

TA2ConfigParser::TA2ConfigParser(const char* fname)// Use c-style string soI don't have to include <string> just yet
{
    // Set member data and some utility bits
    pi = boost::math::constants::pi<long double>();
    configGood = true;
    genBmap = false;
    std::string line, plane;
    
    bool detectorsPresent(false);
    std::vector<DetectorMagnet*> dets;
    
    std::ifstream config(fname);// Open the config file
    if(config.is_open())
    {
        while(std::getline(config, line))
        {
            std::vector<std::string> linesplit;
            boost::split(linesplit, line, boost::is_any_of(", "));
            if(line.compare("DIPOLE") == 0)
            {
                while(std::getline(config, line) && line.compare("END") != 0)
                {
                    std::vector<std::string> linesplit;
                    boost::split(linesplit, line, boost::is_any_of(", "));
                    if(linesplit[0].compare("INNERR") == 0)
                    {
                        innerR = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("OUTERR") == 0)
                    {
                        outerR = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("STARTA") == 0)
                    {
                        startA = static_cast<long double>(atof(linesplit[1].c_str())*pi/180.0);
                    }
                    else if(linesplit[0].compare("ENDA") == 0)
                    {
                        endA = static_cast<long double>(atof(linesplit[1].c_str())*pi/180.0);
                    }
                    else if(linesplit[0].compare("GAP") == 0)
                    {
                        gap = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("CENTRE") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        centre = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("B0") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        B0 = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("NAME") == 0)
                    {
                        name = linesplit[1];
                    }
                    else
                    {
                       continue;
                    }
                    
                }
                magnet* dipole = new HardEdgedArcDipole(name, innerR, outerR, startA, endA, gap, *B0, *centre);
                magnets.push_back(dipole);
            }// End Dipole
            else if(line.compare("255SPECTROMETER") == 0)
            {
                while(std::getline(config, line) && line.compare("END") != 0)
                {
                    std::vector<std::string> linesplit;
                    boost::split(linesplit, line, boost::is_any_of(", "));
                    
                    if(linesplit[0].compare("INNERR") == 0)
                    {
                        innerR = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("OUTERR") == 0)
                    {
                        outerR = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("STARTA") == 0)
                    {
                        startA = static_cast<long double>(atof(linesplit[1].c_str())*pi/180.0);
                    }
                    else if(linesplit[0].compare("ENDA") == 0)
                    {
                        endA = static_cast<long double>(atof(linesplit[1].c_str())*pi/180.0);
                    }
                    else if(linesplit[0].compare("GAP") == 0)
                    {
                        gap = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("ALPHA") == 0)
                    {
                    
                        alpha = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("BETA") == 0)
                    {
                        betaDS = static_cast<long double>(atof(linesplit[1].c_str()));

                    }
                    else if(linesplit[0].compare("CENTRE") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        centre = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("B0") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        B0 = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("NAME") == 0)
                    {
                        name = linesplit[1];
                    }
                    else
                    {
                       continue;
                    }
                    
                }
                magnet* Shull_Dennison = new HardEdged225Spectrometer(name, innerR, outerR, startA, endA, gap, *B0, *centre, alpha, betaDS);
                magnets.push_back(Shull_Dennison);
            }// End Spectrometer (Shull-Dennison)
            else if(line.compare("GENFIELDMAP") == 0)
            {
                while(std::getline(config, line) && line.compare("END") != 0)
                {
                    std::vector<std::string> linesplit;
                    boost::split(linesplit, line, boost::is_any_of(", "));
                    if(linesplit[0].compare("EXTX") == 0)
                    {
                        extX = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("EXTY") == 0)
                    {
                        extY = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("EXTZ") == 0)
                    {
                        extZ = static_cast<long double>(atof(linesplit[1].c_str()));
                    }
                    else if(linesplit[0].compare("MAPFILE") == 0)
                    {
                        mapfile = linesplit[1];
                    }
                    else if(linesplit[0].compare("CENTRE") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        centre = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("NAME") == 0)
                    {
                        name = linesplit[1];
                    }
                    else
                    {
                       continue;
                    }
                    
                }
                magnet* GenFieldMappd = new GeneralMagnetFromMap(name, extX, extY, extZ, *centre, mapfile);
                magnets.push_back(GenFieldMappd);
            }// End General field mapped element
            else if(line.compare("DETECTOR") == 0)
            {
                while(std::getline(config, line) && line.compare("END") != 0)
                {
                    std::vector<std::string> linesplit;
                    boost::split(linesplit, line, boost::is_any_of(", "));
                    if(linesplit[0].compare("EXT") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        ext_det = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("OUTLOC") == 0)
                    {
                        outloc_det = linesplit[1];
                    }
                    else if(linesplit[0].compare("CENTRE") == 0)
                    {
                        long double x(0), y(0), z(0);
                        x = static_cast<long double>(atof(linesplit[1].c_str()));
                        y = static_cast<long double>(atof(linesplit[2].c_str()));
                        z = static_cast<long double>(atof(linesplit[3].c_str()));
                        centre = new ThreeVector(x,y,z);
                    }
                    else if(linesplit[0].compare("NAME") == 0)
                    {
                        name = linesplit[1];
                    }
                    else
                    {
                       continue;
                    }
                    
                }
                magnet* Detector = new DetectorMagnet(name, outloc_det, *centre, *ext_det);
                magnets.push_back(Detector);
                dets.push_back((DetectorMagnet*)Detector);
                detectorsPresent = true;
            }// End DetectorMagnet
            else if(linesplit[0].compare("PHASESPACE") == 0)
            {
                phasespace = new TAPhasespace(linesplit[1]);// This is all we will do with it now, read in later...
                usePhasespace = true;
            }
            else if(linesplit[0].compare("BMAP") == 0)
            {
                outloc = linesplit[1];
                granularity = static_cast<long double>(atof(linesplit[2].c_str()));
                genBmap = true;
            }
            else if(linesplit[0].compare("STEP") == 0)
            {
                step = static_cast<long double>(atof(linesplit[1].c_str()));
            }
            else if(linesplit[0].compare("NSTEPS") == 0)
            {
                nsteps = atoi(linesplit[1].c_str());
            }
            else if(linesplit[0].compare("NCORES") == 0)
            {
                ncores = atoi(linesplit[1].c_str());
            }
            else if(linesplit[0].compare("SUPPRESS") == 0)
            {
                suppress_output = true;
            }
            else if(linesplit[0].compare("LIMITS") == 0)
            {
                long double x(0), y(0), z(0);
                x = static_cast<long double>(atof(linesplit[1].c_str()));
                y = static_cast<long double>(atof(linesplit[2].c_str()));
                z = static_cast<long double>(atof(linesplit[3].c_str()));
                limits = new ThreeVector(x,y,z);
            }
            else if(linesplit[0].compare("XI") == 0)
            {
                long double x(0), y(0), z(0);
                x = static_cast<long double>(atof(linesplit[1].c_str()));
                y = static_cast<long double>(atof(linesplit[2].c_str()));
                z = static_cast<long double>(atof(linesplit[3].c_str()));
                xi = new ThreeVector(x,y,z);
                configPhaseSpacexi.push_back(*xi);
            }
            else if(linesplit[0].compare("VI") == 0)
            {
                long double x(0), y(0), z(0);
                x = static_cast<long double>(atof(linesplit[1].c_str()));
                y = static_cast<long double>(atof(linesplit[2].c_str()));
                z = static_cast<long double>(atof(linesplit[3].c_str()));
                vi = new ThreeVector(x,y,z);
                configPhaseSpacevi.push_back(*vi);
            }
            else if(linesplit[0].compare("OUTLOC") == 0)
            {
                outloc = linesplit[1];
            }
            else if(linesplit[0].compare("ELEMMAP") == 0)
            {
                elements = true;
                nbins = atoi(linesplit[1].c_str());
                elmapoutloc = linesplit[2];

            }
            else
            {
                continue;
            }
        }// Close on while for config file
    }// Close config read
    
    if(detectorsPresent)
    {
        std::vector<DetectorMagnet*>::iterator curr;
        for(curr=dets.begin(); curr!=dets.end(); ++curr)// By now, we should know how many particles to expect, so we tell the detectors
        {
            (*curr)->SetNpart(GetNpart());
        }
    }
    config.close();
}

std::pair<ThreeVector, ThreeVector> TA2ConfigParser::GetParticle()
{
    if(usePhasespace)
    {
        return phasespace->GetParticle();
    }
    else if(configPhaseSpacexi.size() > 0)
    {
        std::pair<ThreeVector, ThreeVector> rval;
        rval = std::make_pair(configPhaseSpacexi[0], configPhaseSpacevi[0]);
        configPhaseSpacexi.erase(configPhaseSpacexi.begin());
        configPhaseSpacevi.erase(configPhaseSpacevi.begin());
        return rval;
    }
    std::cerr << "Something went wrong... Returning garbage" << std::endl;
    ThreeVector v0f;
    ThreeVector x0f;
    return std::make_pair(x0f, v0f);
}


int TA2ConfigParser::GetNpart()
{
    if(usePhasespace)
    {
        return phasespace->GetNpart();
    }
    else
    {
        return configPhaseSpacexi.size();
    }
}
