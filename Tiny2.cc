// Rewriting the TinyArcs program to use my own version of a three-vector library
// We will still need boost though...

#include <iostream>
#include <vector>
#include <string>
#include <utility>
//#include <omp.h>


#include "ThreeVector.hh"
#include "TA2ConfigParse.hh"
#include "TA2Util.hh"


using namespace std;
using namespace tau;

int main(int argc, char* argv[])
{
    // Read the commandline, decide what to do with it
    if(argc == 1)// user wants help
    {
        cout << "Help for TA2" << endl;
        cout << "Usage:" << endl;
        cout << "./TA2 <config file>" << endl;
        cout << "config file - A file containing the details of your magnetic system" << endl << endl;
        cout << "Allowed config file parameters:"<< endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "DIPOLE - Beginning of the dipole element definition, must include all the following:" << endl;
        cout << "\t-INNERR -> The inner radius of curvature" << endl;
        cout << "\t-OUTERR -> The outer radius of curvature" << endl;
        cout << "\t-STARTA -> The angle the dipole starts at" << endl;
        cout << "\t-ENDA -> The angle the dipole ends at" << endl;
        cout << "\t-GAP -> The magnet gap in z-direction" << endl;
        cout << "\t-CENTRE -> The 3D coordinates on the centre of the magnet" << endl;
        cout << "\t-B0 -> The 3D coodrinates of the magnetic field" << endl;
        cout << "\t-NAME -> A descriptive name to give the magnet" << endl;
        cout << "\t-END -> Keyword signifying the end of the magnet definition" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "225SPECTROMETER - Beginning of the DennisonShull element definition, must include all the following:" << endl;
        cout << "\t-INNERR -> The inner radius of curvature" << endl;
        cout << "\t-OUTERR -> The outer radius of curvature" << endl;
        cout << "\t-STARTA -> The angle the dipole starts at" << endl;
        cout << "\t-ENDA -> The angle the dipole ends at" << endl;
        cout << "\t-GAP -> The magnet gap in z-direction" << endl;
        cout << "\t-ALPHA -> One of the parameters which can be tuned in this design" << endl;
        cout << "\t-BETA -> The second tuning parameter. See the paper by Dennison + Shull" << endl;
        cout << "\t-CENTRE -> The 3D coordinates on the centre of the magnet" << endl;
        cout << "\t-B0 -> The 3D coodrinates of the magnetic field" << endl;
        cout << "\t-NAME -> A descriptive name to give the magnet" << endl;
        cout << "\t-END -> Keyword signifying the end of the magnet definition" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "GENFIELDMAP - Beginning of the dipole element definition, must include all the following:" << endl;
        cout << "\t-EXTX -> The extent of the map in x" << endl;
        cout << "\t-EXTY -> The extent of the map in x" << endl;
        cout << "\t-EXTZ -> The extent of the map in x" << endl;
        cout << "\t-ENDA -> The angle the dipole ends at" << endl;
        cout << "\t-MAPFILE -> The magnet field data, must be (x,y,z,Bx,By,Bz)" << endl;
        cout << "\t-CENTRE -> The 3D coordinates on the centre of the magnet" << endl;
        cout << "\t-NAME -> A descriptive name to give the magnet" << endl;
        cout << "\t-END -> Keyword signifying the end of the magnet definition" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "DETECTOR:" << endl;
        cout << "\t-EXT -> The extent of the detector, specified as a 3-vector" << endl;
        cout << "\t-OUTLOC -> Specify the output location of tracks through the detector" << endl;
        cout << "\t-CENTRE -> The 3D coordinates on the centre of the detector" << endl;
        cout << "\t-NAME -> A descriptive name to give the detector" << endl;
        cout << "\t-END -> Keyword signifying the end of the detector definition" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        cout << "General options:" << endl;
        cout << "\t-PHASESPACE -> Name of a file containing a set of data, must be (x,y,z, vx,vy,vz)" << endl;
        cout << "\t-BMAP -> Produce a file mapping the magnetic field. Keyword, then output location and granularity (i.e. pixel size in m)" << endl;
        cout << "\t-STEP -> Specify the step length (in m)" << endl;
        cout << "\t-NSTEPS -> Specify the maximum number of steps to take" << endl;
        cout << "\t-NCORES -> Specify the number of threads to launch. Careful, it is possible to melt your computer with this!" << endl;
        cout << "\t-SUPPRESS -> Suppress all command line output" << endl;
        cout << "\t-LIMITS -> Specify the limits to which we will track, outside this, the particle is killed. Should be specified like a vector x,y,z" << endl;
        cout << "\t-XI -> If not using a phase space, you can track individual particles by specifying their initial position (x,y,z) and..." << endl;
        cout << "\t-VI -> ...their initial velocity (vx,vy,vz)" << endl;
        cout << "\t-OUTLOC -> Specify the output name/directory" << endl;
        cout << "\t-ELEMMAP -> Generate a quick map to show where stuff is, but not fields" << endl;
        cout << "------------------------------------------------------------------------------------" << endl;
        
        cout << "For further details, try to find the example config file I made..." << endl;


        return -1;
    }
    else if(argc == 2)// Everything now has to go within this if statement
    {
        TA2ConfigParser config(argv[1]);
        if(config.DoBmap())
        {
            tau::GenerateFieldMap(config.GetMagnets(), config.GetLimits(), config.GetBmapOutloc(), config.GetGranularity(), config.GetNcores());
            return 0;
        }
        else if(config.MakeElemMap())
        {
            tau::GenerateElemMap(config.GetMagnets(), config.GetLimits(), config.GetOutlocForElmap(), config.GetNbins(), config.GetNcores());
            return 0;
        }
        std::vector<std::pair<ThreeVector, ThreeVector> >* resn[config.GetNpart()];
        std::vector<std::vector<std::pair<ThreeVector, ThreeVector> >* > reslist;
        for(int i=0; i<config.GetNpart(); ++i)
        {
            resn[i] = new std::vector<std::pair<ThreeVector, ThreeVector> >;
            reslist.push_back(resn[i]);
        }
        tau::RunTracking(&config, reslist);
        cout << "Tracking finished!" << endl;

    }
    else
    {
        cout << "Incorrect number of arguments specified! Aborting..." << endl;
        return -2;
    }
    return 0;
}