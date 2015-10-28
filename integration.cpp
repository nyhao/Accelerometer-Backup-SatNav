#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstring>

using namespace std;

double t, dt, tn, accmeter, gps, va, pa, vg, pg, noise;
double F (double x, double y, double zzz) {return 0.;}

// Testing Runge-Kutta 4th Order Integration Method
int main()
{
    double fg, k1g, k2g, k3g, k4g, fa, k1a, k2a, k3a, k4a;
    double filsam, filcut, filsf;       //for FIR filter
    int dropoff = 0, pickup = 0;        //for switching algorithm
    //int q1 = 0, q2 = 0, q3 = 0;
    //int seed = 1;
    srand(time(NULL));
    dropoff = (rand()%50)+1;            //for switching algorithm
    while (pickup <= dropoff) {pickup = (rand()%50)+1;}     //for switching algorithm
    ofstream myfile;
    myfile.open("record.txt");
    t = 0.0;
    tn = 50.0;
    dt = 0.0002;
    accmeter = 0.0;
    gps = 0.0;
    vg = 0.0;
    pg = 0.0;
    va = 0.0;
    pa = 0.0;
    filsam = 50.0;
    filcut = 1.0;
    filsf = 5000.0;
    noise = 0.0;
    myfile << setprecision(6) << t << "     " << gps << "     " << vg <<  "     " << pg << "     " << endl;
    t = t+dt;
    for (; t < tn; t = t+dt)
    {
        gps = F(t,vg,0.);
        //to create histogram
        /*
        if(a <= 0.011 && a >= -0.011)
        {
            q1 += 1;
            q2 += 1;
            q3 += 1;
        }
        else if(a <= 0.022 && a >= -0.022)
        {
            q2 += 1;
            q3 += 1;
        }
        else if(a <= 0.033 && a >= -0.033)
        {
            q3 += 1;
        }
        */
        //gps readings
        fg = gps;
        k1g = fg;
        fg = F(t+(dt/2.0),vg+(dt/2.0)*k1g,0.);
        k2g = fg;
        fg = F(t+(dt/2.0),vg+(dt/2.0)*k2g,0.);
        k3g = fg;
        fg = F(t+dt,vg+dt*k3g,0.);
        k4g = fg;
        vg += (k1g+2.0*k2g+2.0*k3g+k4g)*(dt/6.0);
        k1g = k2g = k3g = k4g = vg;
        pg += (k1g+2.0*k2g+2.0*k3g+k4g)*(dt/6.0);
        myfile << setprecision(6) << t << "     " << gps << "     " << vg <<  "     " << pg << "     " << endl;
    }
    myfile.close();
    //cout << dropoff << endl << pickup;
    //cout << q1 << endl << q2 << endl << q3 << endl;
    return 0;
}
