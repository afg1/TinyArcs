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
    
    pthread_mutex_lock(&getPartMutex);
    TA2ConfigParser* conf = input->conf;
    input->initPair = conf->GetParticle();
    input->res->push_back(input->initPair);
    int nsteps = conf->GetNsteps();
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
        t_args[i]->lastStep = false;
        t_args[i]->firstStep = false;
        
    }
    if(npart > conf->GetNcores())
    {
        for(int i=0; i<npart; i+=conf->GetNcores())
        {
            if(npart - i > conf->GetNcores())
            {
                for(int j=0; j< conf->GetNcores(); ++j)
                {
                    pthread_create(&jobs[i+j], NULL, tau::consumer, (void*)t_args[i+j]);
                }
                for(int j=0; j<conf->GetNcores(); ++j)
                {
                    pthread_join(jobs[i+j], NULL);
                }
            }
            else
            {
                for(int j=0; j< (npart - i); ++j)
                {
                    pthread_create(&jobs[i+j], NULL, tau::consumer, (void*)t_args[i+j]);
                }
                for(int j=0; j<(npart - conf->GetNcores()); ++j)
                {
                    pthread_join(jobs[i+j], NULL);
                }
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
        fname << conf->GetOutloc() << "particle_" << i << ".txt";
        std::ofstream out(fname.str().c_str(), std::ios::binary);
        for(int j=0; j< t_args[i]->res->size() ; ++j)
        {
            t_args[i]->res->operator[](j).first.BinaryWriteToStream(out);
            t_args[i]->res->operator[](j).second.BinaryWriteToStream(out);
//            out << std::setprecision(30) << t_args[i]->res->operator[](j).first << "\t" << t_args[i]->res->operator[](j).second << std::endl;
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

    long double endA(0), startA(0);
    bool eligibleForNR(false);
    ThreeVector normal;
    ThreeVector planePoint;
    if(tau::GetInMagnet(x0, magnets) && tau::EligibleForNR(x0, magnets))
    {
        endA = tau::GetEndA(x0, magnets);
        startA = tau::GetStartA(x0, magnets);
        normal = tau::GetNormal(x0, magnets);
        planePoint = tau::GetPlanePoint(x0, magnets);
        eligibleForNR = true;
    }
    long double phi0 = std::atan2(x0.GetElem(1), x0.GetElem(0));// By trial and error, this ends up right
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
    if(phi0 + theta > endA  && !arg->lastStep && eligibleForNR)
    {
        ThreeVector x1 = x0 + Fu*(r-r*cos(thetat)) + vperpu*(r*sin(thetat)) + vpar*t0;
        do
        {
            thetat = omega*t0;
            x1 = x0 + Fu*(r-r*cos(thetat)) + vperpu*(r*sin(thetat)) + vpar*t0;
            
            ThreeVector brackets = (x1 - planePoint);
            long double topline = brackets.Dot(normal);
            ThreeVector bottomBrackets = Fu*vperpMag * sin(thetat) + vperpu*(vperpMag * cos(thetat));
            long double bottom = bottomBrackets.Dot(normal);
            t1 = t0 - topline/bottom;
            deltat = topline/bottom;
            t0 = t1;
            niter++;
//            std::cerr << topline << "\t" << bottom << "\t" << deltat/dt <<  std::endl;
        }while(fabs(deltat/dt) > 1E-5);
        arg->lastStep = true;
        std::pair<ThreeVector, ThreeVector> endpair(x1, vr);
        arg->initPair = endpair;
        arg->res->push_back(endpair);
        return;
    }
    else if(phi0 + theta < startA  && !arg->lastStep && eligibleForNR)
    {
        ThreeVector x1 = x0 + Fu*(r-r*cos(thetat)) + vperpu*(r*sin(thetat)) + vpar*t0;
        do
        {
            thetat = omega*t0;
            x1 = x0 + Fu*(r-r*cos(thetat)) + vperpu*(r*sin(thetat)) + vpar*t0;
            
            ThreeVector brackets = (x1 - planePoint);
            long double topline = brackets.Dot(normal);
            ThreeVector bottomBrackets = Fu*vperpMag * sin(thetat) + vperpu*(vperpMag * cos(thetat));
            long double bottom = bottomBrackets.Dot(normal);
            t1 = t0 - topline/bottom;
            deltat = topline/bottom;
            t0 = t1;
            niter++;
//            std::cerr << topline << "\t" << bottom << "\t" << deltat/dt <<  std::endl;
        }while(fabs(deltat/dt) > 1E-5);
        arg->lastStep = true;
        std::pair<ThreeVector, ThreeVector> endpair(x1, vr);
        arg->initPair = endpair;
        arg->res->push_back(endpair);
        return;
    }

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

long double tau::GetStartA(ThreeVector x, std::vector<magnet*> magnets)
{
        std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        if((*curr)->InMagnet(x))
        {
            return (*curr)->GetStartA();
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

bool tau::EligibleForNR(ThreeVector x, std::vector<magnet*> magnets)
{
    std::vector<magnet*>::iterator curr;
    for(curr = magnets.begin(); curr != magnets.end(); ++curr)
    {
        if((*curr)->InMagnet(x))
        {
            return (*curr)->Eligible();
        }
    }
    return false;
}



void tau::GenerateFieldMap(std::vector<magnet*> magnets, ThreeVector limits, std::string outloc, long double granular)
{
    ThreeVector rp;
    ThreeVector Bp(0.0, 0.0, 1.0);
    long double xp(-limits.GetElem(0)), yp(-limits.GetElem(1)), zp(-limits.GetElem(2));
    double fieldPoint[6] = {0.0,0.0,0.0,0.0, 0.0, 0.0};
    std::ofstream out;
    out.open(outloc.c_str(), std::ios::binary);
    if(out.is_open())
    {
        for(xp=-1*limits.GetElem(0); xp <= limits.GetElem(0); xp+= granular)
        {
            for(yp=-1*limits.GetElem(1); yp <= limits.GetElem(1); yp+= granular)
            {
                for(zp=-1*limits.GetElem(2); zp <= limits.GetElem(2); zp+= granular)
                {
                    rp.SetElem(0, xp);
                    rp.SetElem(1, yp);
                    rp.SetElem(2, zp);
                    Bp = tau::GenerateBmap(rp, magnets);
                    fieldPoint[0] = xp;
                    fieldPoint[1] = yp;
                    fieldPoint[2] = zp;
                    fieldPoint[3] = Bp.GetElem(0);
                    fieldPoint[4] = Bp.GetElem(1);
                    fieldPoint[5] = Bp.GetElem(2);
                    out.write((char*)&fieldPoint, sizeof(fieldPoint));
                }
            }
        }
    }
    else
    {
        std::cerr << "Unable to open file for Bmap output" << std::endl;
    }
}
    

