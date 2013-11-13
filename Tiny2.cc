// Rewriting the TinyArcs program to use my own version of a three-vector library
// We will still need boost though...

#include <iostream>
#include <vector>
#include <string>
#include <utility>


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
        cout << "config file - A file containing the details of your magnetic system" << endl;
        return -1;
    }
    else if(argc == 2)// Everything now has to go within this if statement
    {
        TA2ConfigParser config(argv[1]);
        tau::TAargs args;
        args.sequence = config.GetMagnets();
        args.limits = config.GetLimits();
        if(args.UsePhasespace)
        {
            args.phasespace = config.GetPhaseSpace();
        }
        args.step = config.GetStep();
        args.nsteps = config.GetNsteps();// Build the args struct from the config object

        vector<pair<ThreeVector, ThreeVector> > res;
        tau::RunTracking(args, &res);
        
    }
    else
    {
        cout << "Incorrect number of arguments specified! Aborting..." << endl;
        return -2;
    }
    
    return 0;
}