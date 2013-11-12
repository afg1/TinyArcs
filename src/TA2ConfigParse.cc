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
    std::string line;
    
    std::ifstream config(fname);// Open the config file
    if(config.is_open())
    {
        while(std::getline(config, line))
        {
            std::vector<std::string> linesplit;
            boost::split(linesplit, line, boost::is_any_of(", "));
            if(line.compare("DIPOLE") == 0)
            {
                std::cout << "Standard Dipole" << std::endl;
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
                    else
                    {
                       continue;
                    }
                    
                }
                magnet* dipole = new HardEdgedArcDipole("dipole", innerR, outerR, startA, endA, gap, B0, centre);
                magnets.push_back(dipole);
            }// End Dipole
            else if(line.compare("255SPECTROMETER") == 0)
            {
                std::cout << "Spectrometer using Shull-Dennison parameterisation" << std::endl;
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
                    else
                    {
                       continue;
                    }
                    
                }
                magnet* Shull_Dennison = new HardEdged225Spectrometer("Shull-Dennson", innerR, outerR, startA, endA, gap, B0, centre, alpha, betaDS);
                magnets.push_back(Shull_Dennison);
            }// End Spectrometer (Shull-Dennison)
            else if(linesplit[0].compare("PHASESPACE") == 0)
            {
                phasespace = linesplit[1];// This is all we will do with it now, read in later...
            }
            else if(linesplit[0].compare("STEP") == 0)
            {
                step = static_cast<long double>(atof(linesplit[1].c_str()));
            }
            else if(linesplit[0].compare("NSTEPS") == 0)
            {
                nsteps = atoi(linesplit[1].c_str());
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
            }
            else if(linesplit[0].compare("VI") == 0)
            {
                long double x(0), y(0), z(0);
                x = static_cast<long double>(atof(linesplit[1].c_str()));
                y = static_cast<long double>(atof(linesplit[2].c_str()));
                z = static_cast<long double>(atof(linesplit[3].c_str()));
                vi = new ThreeVector(x,y,z);
            }
            else
            {
                continue;
            }
        }// Close on while for config file
    }// Close config read
    config.close();
}
    

