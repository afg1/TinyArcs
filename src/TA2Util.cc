#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <utility>
#include <cstdlib>
#include <string.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <sstream>


#include <pthread.h>
#include <boost/math/constants/constants.hpp>


#include "QueueJumper.hh"

#include "TAPhasespace.hh"
#include "TA2ConfigParse.hh"
#include "TA2Util.hh"
#include "Magnets.hh"

pthread_mutex_t getPartMutex = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t threadWaitMutex = PTHREAD_MUTEX_INITIALIZER;


void *tau::consumer(void* arg)
{
/*
This function will do the actual traking
*/
    tau::ThreadArgs* input = (tau::ThreadArgs*)arg;
    TA2ConfigParser* conf = input->conf;
    
    pthread_mutex_lock(&getPartMutex);
    input->initPair = conf->GetParticle();
    input->res->push_back(input->initPair);
    int nsteps =conf->GetNsteps();
    pthread_mutex_unlock(&getPartMutex);
    
    for(int n=0; n<nsteps; ++n)
    {
        TrackingStep(input);
        if(input->initPair.first.Magnitude() > conf->GetLimits().Magnitude())
        {
            break;
        }
    }
    pthread_exit(NULL);
}




void tau::RunTracking(TA2ConfigParser* conf, std::vector<std::vector<std::pair<ThreeVector, ThreeVector> >* >& reslist)
{
    
    const int npart = conf->GetNpart();
    pthread_t jobs[npart];// array of pthreads, which are what I'm calling "jobs"
    tau::ThreadArgs* t_args[npart];
    // Need to use struct to pass args into a pthread
    
    for(int i=0; i< npart; ++i)// this loop builds the job queue
    {
        t_args[i] = new ThreadArgs;
        t_args[i]->conf = conf;
        t_args[i]->res = reslist[i];
    }
    if(npart > conf->GetNcores())
    {
        for(int i=0; i<npart; i+=conf->GetNcores())
        {
                for(int j=0; j< conf->GetNcores(); ++j)
                {
                    pthread_create(&jobs[i+j], NULL, tau::consumer, (void*)t_args[i+j]);
                }
                for(int j=0; j<conf->GetNcores(); ++j)
                {
                    pthread_join(jobs[i], NULL);
                }
        }
    }
    else
    {
        for(int i=0; i< npart; ++i)// this loop builds the job queue
        {
            pthread_create(&jobs[i], NULL, tau::consumer, (void*)t_args[i]);
        }
        for(int i=0; i< npart; ++i)// this loop joins the treads
        {
            pthread_join(jobs[i], NULL);
        }
    }

    std::stringstream fname;
    for(int i=0; i< npart; ++i)
    {
        fname << "particle_" << i << ".txt";
        std::ofstream out(fname.str().c_str());
        for(int j=0; j< t_args[i]->res->size() ; ++j)
        {
            out << t_args[i]->res->operator[](j).first << "\t" << t_args[i]->res->operator[](j).second << std::endl;
        }
        out.close();
        fname.str("");
    }
}


long double tau::gammaFromE(long double T)
{
    long double MP_MEV(938.27204621);
    long double g = (long double)1 + (T/MP_MEV);// The proton mass in MeV is defined in the header
    return g;
}

long double tau::betaFromE(long double T)
{
    long double g = gammaFromE(T);
    long double b = sqrtl(1 - pow(g,-2));
    return b;
}

long double tau::gammaFromV(long double V)
{
    long double b = tau::betaFromV(V);
    long double g = 1.0/(sqrtl(1.0 - (b*b)));
    return g;
}

long double tau::betaFromV(long double V)
{
    long double c(299792458);
    long double b = V/c;
    return b;
}

ThreeVector tau::GenerateBmap(ThreeVector x, std::vector<magnet*> magnets)
{
    ThreeVector Bret(0.0, 0.0, 0.0);
    std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        Bret += (*curr)->B(x);
    }
    return Bret;
}

void tau::TrackingStep(tau::ThreadArgs* arg)
{
    long double pi = boost::math::constants::pi<long double>();
    std::vector<magnet*> magnets = arg->conf->GetMagnets();
    ThreeVector x0 = arg->initPair.first;
    ThreeVector v0 = arg->initPair.second;
    ThreeVector zero(0.0, 0.0, 0.0);
    long double Q_E (1.60217656535E-19);
    long double MP(1.67262158E-27);
    long double g = gammaFromV(v0.Magnitude());
    
    // Now we can do the tinyarcs method!
    long double dt = arg->conf->GetStep();
    
    ThreeVector x = x0 + (v0*dt/2);
    
    
    ThreeVector B = tau::GenerateBmap(x, magnets);
    if(B == zero)
    {
        std::pair<ThreeVector, ThreeVector> endpair(x, v0);
        arg->initPair = endpair;
        arg->res->push_back(endpair);
        return;
    }

    long double Bmag = B.Magnitude();
    ThreeVector Bu = B.Normalise();

    
    ThreeVector vpar = Bu*Bu.Dot(v0);
    ThreeVector vperp = v0 - vpar;
    ThreeVector vperpu = vperp.Normalise();
    long double vperpMag = vperp.Magnitude();
    
    ThreeVector F = vperp.Cross(B)*Q_E;// F = q(v x B)
    ThreeVector Fu = F.Normalise();
    long double Fmag = F.Magnitude();
    
    long double r = (g * MP * vperpMag)/(Q_E*Bmag);
    long double omega = vperpMag/r;
    long double theta = omega*dt;
    
    ThreeVector xr;
    ThreeVector vr;
    xr = x0 + Fu*(r-r*cos(theta)) + vperpu*(r*sin(theta)) + vpar*dt;// This should be the same as from the paper...
    vr = Fu*(r*omega*sin(theta)) + vperpu*(r*omega*cos(theta)) + vpar;

    long double endA;
    ThreeVector normal;
    ThreeVector planePoint;
    if(tau::GetInMagnet(x0, magnets))
    {
        endA = tau::GetEndA(x0, magnets);
        normal = tau::GetNormal(x0, magnets);
        planePoint = tau::GetPlanePoint(x0, magnets);
    }
    long double phi0 = theta = std::atan2(xr.GetElem(1), xr.GetElem(0));// By trial and error, this ends up right
    if(phi0 < 0)
    {
        phi0 += 2*pi;
    }
    long double thetat(theta);
    long double t0 = dt;
    long double t1(0);
    long double deltat(0);
    int niter(0);
    // Newton Raphston doesn't work yet...
//    if(phi0 + theta > endA)
//    {
//        ThreeVector x1 = x0 + Fu*(r-r*cos(thetat)) + vperpu*(r*sin(thetat)) + vpar*t0;
//        do
//        {
//            ThreeVector brackets = (x1 - planePoint);
//            long double topline = brackets.Dot(normal);
//            long double bottom = (Fu*vperpMag * sin(phi0+thetat) + vperpu*(vperpMag * cos(phi0+thetat))).Dot(normal);
//            t1 = t0 - topline/bottom;
//            deltat = topline/bottom;
//            t0 = t1;
//            thetat = omega*t0;
//            niter++;
//            std::cerr << deltat/dt << "\t" << dt << std::endl;
//        }while(fabs(deltat/dt) > 10);
//        
//        std::pair<ThreeVector, ThreeVector> endpair(x1, vr);
//        arg->initPair = endpair;
//        arg->res->push_back(endpair);
//        return;
//    }
    std::pair<ThreeVector, ThreeVector> endpair(xr, vr);
    arg->initPair = endpair;
    arg->res->push_back(endpair);
}




// Here be functions for getting sense from the magnet deities

long double tau::GetEndA(ThreeVector x, std::vector<magnet*> magnets)
{
        std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        if((*curr)->InMagnet(x))
        {
            return (*curr)->GetEndA();
        }
    }
    return 0;
}

ThreeVector tau::GetNormal(ThreeVector x, std::vector<magnet*> magnets)
{
    std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        if((*curr)->InMagnet(x))
        {
            return (*curr)->GetNormal();
        }
    }
}

ThreeVector tau::GetPlanePoint(ThreeVector x, std::vector<magnet*> magnets)
{
    std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        if((*curr)->InMagnet(x))
        {
            return (*curr)->GetPlanePoint();
        }
    }
}


bool tau::GetInMagnet(ThreeVector x, std::vector<magnet*> magnets)
{
    std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        if((*curr)->InMagnet(x))
        {
            return true;
        }
    }
    return false;
}



