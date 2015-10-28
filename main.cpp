#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <cstring>

using namespace std;

double t, dt, tn, accmeter, gps, va, pa, vg, pg, noise;
// Example acceleration function
double F (double x, double y, double zzz) {return 3*x+zzz;}
// Uniform distributed random number generator
double unirand () {return (double)(rand())/(double)(RAND_MAX);}
// Gaussian distributed random number generator
double gausrand (double mean, double stddev)
{
    double u1, u2, v1 = 0, v2 = 0, s;
    int IsSpareReady = 0;
    if (IsSpareReady == 1)
    {
        IsSpareReady = 0;
        return mean+stddev*v2;
    }
    else {
        do {
            u1 = 2.0*unirand()-1.0;
            u2 = 2.0*unirand()-1.0;
            s = u1*u1+u2*u2;
        } while (s >= 1.0);
        v1 = sqrt((-2.0*log(s))/s)*u1;
        v2 = sqrt((-2.0*log(s))/s)*u2;
        return mean+stddev*v1;
    }
}

// FIR Window Design Method (***NOT WORKING***)
/*double FIR (double m, double fc, double fs, double p, double q, double buzz)
{
    double op = 0.0;
    double sumh = 0.0;
    for (double n = (-1)*(m/2); n < m/2; n++)
    {
        if (n != 0)
        {
            double h = (fs/n)*(1/M_PI)*sin(2*M_PI*fc*(n/fs));
            op += (h*(F((((t/dt)-n)*dt)+p,va+q,buzz)));
            sumh += h;
        }
        else
        {
            op += (2*fc*(F((((t/dt)-n)*dt)+p,va+q,buzz)));
            sumh += 2*fc;
        }
    }
    return op/sumh;
} */

/*double FIR (double m, double fc, double fs, double p, double q)
{
    double op = 0.0;
    double sumh = 0.0;
    for (double n = (-1)*(m/2); n < m/2; n++)
    {
        if (n != 0)
        {
            double h = (fs/n)*(1/M_PI)*sin(2*M_PI*fc*(n/fs));
            op += (h*(F((((t/dt)-n)*dt)+p,va+q,noise)));
            sumh += h;
        }
        else
        {
            op += (2*fc*(F((((t/dt)-n)*dt)+p,va+q,noise)));
            sumh += 2*fc;
        }
    }
    return op/sumh;
}*/
/*
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
    for (; t < tn; t = t+dt)
    {
        noise = gausrand(0,0.011);
        accmeter = FIR(filsam,filcut,filsf,0.0,0.0);
        gps = F(t,vg,noise);
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
        /*fg = gps;
        k1g = fg;
        fg = F(t+(dt/2.0),vg+(dt/2.0)*k1g,noise);
        k2g = fg;
        fg = F(t+(dt/2.0),vg+(dt/2.0)*k2g,noise);
        k3g = fg;
        fg = F(t+dt,vg+dt*k3g,noise);
        k4g = fg;
        vg += (k1g+2.0*k2g+2.0*k3g+k4g)*(dt/6.0);
        k1g = k2g = k3g = k4g = vg;
        pg += (k1g+2.0*k2g+2.0*k3g+k4g)*(dt/6.0);
        //accelerometer readings
        /*fa = accmeter;
        k1a = fa;
        fa = FIR(filsam,filcut,filsf,(dt/2.0),(dt/2.0)*k1a);
        k2a = fa;
        fa = FIR(filsam,filcut,filsf,(dt/2.0),(dt/2.0)*k2a);
        k3a = fa;
        fa = FIR(filsam,filcut,filsf,dt,dt*k3a);
        k4a = fa;
        va += (k1a+2.0*k2a+2.0*k3a+k4a)*(dt/6.0);
        k1a = k2a = k3a = k4a = va;
        pa += (k1a+2.0*k2a+2.0*k3a+k4a)*(dt/6.0);*/
        //if (t >= dropoff && t < pickup)
        //{
           // myfile << setprecision(6) << t << "     " << accmeter << "     " << va <<  "     " << pa << "     " << gps << endl;
        //}
  /*      else {myfile << setprecision(6) << t << "     " << gps << "     " << vg <<  "     " << pg << endl;}
    }
    myfile.close();
    cout << dropoff << endl << pickup;
    //cout << q1 << endl << q2 << endl << q3 << endl;

    //to check gaussian distribution
    /*
    int N=100000;
    double mean=0.,variance=0.;
    for(int i=0;i<N;i++)
    {
        double z=gausrand(0, 0.011);
        mean = mean + z;
    }
    mean = mean/N;
    cout << " mean = " << mean << endl;

    for(int i=0;i<N;i++)
    {
        double z=gausrand(0, 0.011);
        variance = variance + (mean - z)*(mean - z);
    }
    variance = variance/N;
    cout << " variance = " << variance << endl;
    */
 /*   return 0;
}
*/

// Least Squares Filter
int main()
{
    double a, b;
    double fg, k1g, k2g, k3g, k4g, fa, k1a, k2a, k3a, k4a;
    double noisynoisy, vn, pn, fn, k1n, k2n, k3n, k4n;
    double avefil, avenoisy, totalfil, totalnoisy, crudefil, crudenoisy;
    srand(time(NULL));
    ofstream myfile;
    ofstream fout;
    myfile.open("roar.txt");
    fout.open("acc.csv");
    tn = 50.0;
    dt = 0.0002;
    int sample = (tn/dt)+1;
    double* lol;
    lol = new double[sample];
    memset(lol,0.,sample*sizeof(double));
    double* haha;
    haha = new double[sample];
    memset(haha,0.,sample*sizeof(double));
    double* sum;
    sum = new double[sample];
    memset(sum,0.,sample*sizeof(double));
    double* rms;
    rms = new double[sample];
    memset(rms,0.,sample*sizeof(double));
    int num_per_bracket = 25000;
    double t1 = 0.0;
    double va0 = 0.0;
    double pa0 = 0.0;
    double r;
    double n;
    avefil = avenoisy = totalfil = totalnoisy = crudefil = crudenoisy = 0.0;
    // Monte-Carlo Simulation
    for (int mnm = 0; mnm < 100; mnm++) {
        double* acc;
        acc = new double[sample+1];
        memset(acc,0.,sample*sizeof(double));
        t = 0.0;
        accmeter = 0.0;
        gps = 0.0;
        noisynoisy = 0.0;
        vn = 0.0;
        pn = 0.0;
        vg = 0.0;
        pg = 0.0;
        va = 0.0;
        pa = 0.0;
        t1 = 0.0;
        va0 = 0.0;
        pa0 = 0.0;
        //fout << setprecision(6) << t1 << "      " << acc[0] << endl;
        //myfile << setprecision(6) << t << "     " << accmeter << "     " << va <<  "     " << pa << endl;
        t = t+dt;
        t1 = t1+dt;
        int bign = 0;
        for (int c = 1; c < sample+1; c = c+num_per_bracket)
        {
            for (int sam = 0; sam < num_per_bracket; sam++)
            {
                acc[c+sam] = F(t1,va,gausrand(0,0.011));
                //fout << setprecision(6) << t1 << "      " << acc[c+sam] << endl;
                t1 = t1+dt;
            }
            double sumx = 0., sumy = 0., sumxx = 0., sumxy = 0.;
            // Performing Least Squares Fitting
            for (int i = 0; i < num_per_bracket; i++)
            {
                sumx += (c+i)*dt;
                sumy += acc[c+i];
                sumxx += ((c+i)*dt)*((c+i)*dt);
                sumxy += ((c+i)*dt)*acc[c+i];
            }
            a = (sumy*sumxx-sumx*sumxy)/(num_per_bracket*sumxx-sumx*sumx);
            b = (num_per_bracket*sumxy-sumx*sumy)/(num_per_bracket*sumxx-sumx*sumx);
            // Performing Runge-Kutta 4th Order Integration
            for (; t < (c+num_per_bracket)*dt; t = t+dt)      //condition have to be c-1 if sampling frequency is 5k, c if sampling frequency is 500
            {
                //accelerometer readings
                accmeter = a+b*t;
                fa = accmeter;
                k1a = fa;
                fa = a+b*(t+(dt/2.0));
                k2a = fa;
                fa = a+b*(t+(dt/2.0));
                k3a = fa;
                fa = a+b*(t+dt);
                k4a = fa;
                va += (k1a+2.0*k2a+2.0*k3a+k4a)*(dt/6.0);
                k1a = k2a = k3a = k4a = va;
                pa += (k1a+2.0*k2a+2.0*k3a+k4a)*(dt/6.0);
                //gps readings
                gps = F(t,vg,0.0);
                fg = gps;
                k1g = fg;
                fg = F(t+(dt/2.0),vg+(dt/2.0)*k1g,0.0);
                k2g = fg;
                fg = F(t+(dt/2.0),vg+(dt/2.0)*k2g,0.0);
                k3g = fg;
                fg = F(t+dt,vg+dt*k3g,0.0);
                k4g = fg;
                vg += (k1g+2.0*k2g+2.0*k3g+k4g)*(dt/6.0);
                k1g = k2g = k3g = k4g = vg;
                pg += (k1g+2.0*k2g+2.0*k3g+k4g)*(dt/6.0);
                //noisy reading
                noise = gausrand(0,0.011);
                noisynoisy = F(t,vg,noise);
                fn = noisynoisy;
                k1n = fn;
                fn = F(t+(dt/2.0),vn+(dt/2.0)*k1n,noise);
                k2n = fn;
                fn = F(t+(dt/2.0),vn+(dt/2.0)*k2n,noise);
                k3n = fn;
                fn = F(t+dt,vn+dt*k3n,noise);
                k4n = fn;
                vn += (k1n+2.0*k2n+2.0*k3n+k4n)*(dt/6.0);
                k1n = k2n = k3n = k4n = vn;
                pn += (k1n+2.0*k2n+2.0*k3n+k4n)*(dt/6.0);
                sum[bign] += (pa-pg)*(pa-pg);
                lol[bign] += (pn-pg)*(pn-pg);
                bign++;
                //if (mnm == 99) myfile << setprecision(6) << t << "     " << accmeter << "     " << va <<  "     " << pa << "      " << pg << "      " << pn << endl;
                //totalfil += accmeter;
                //totalnoisy += noisynoisy;
            }
            /*for (double tt = 0.0002; tt < (num_per_bracket)*dt; tt = tt+dt)
            {
                accmeter = a+b*t;
                va = a*tt+0.5*b*tt*tt+va0;
                pa = 0.5*a*tt*tt+(1/6)*b*tt*tt*tt+va0*tt+pa0;
                myfile << setprecision(6) << t << "     " << accmeter << "     " << va <<  "     " << pa << endl;
                t = t+dt;
            }*/
            //va0 = va;
            //pa0 = pa;
        }
        //myfile.close();
        //fout.close();
        delete[] acc;
        /*avefil += (1./sample)*totalfil;
        avenoisy += (1./sample)*totalnoisy;
        crudefil += (0.5)*(1./sample)*totalfil*tn*tn;
        crudenoisy += (0.5)*(1./sample)*totalnoisy*tn*tn;*/
        //if (mnm == 99) cout << bign << endl;
    }
    /*avefil = (1./1000.)*avefil;
    avenoisy = (1./1000.)*avenoisy;
    crudefil = (1./1000.)*crudefil;
    crudenoisy = (1./1000.)*crudenoisy;
    cout << avefil << endl << avenoisy << endl << crudefil << endl << crudenoisy << endl; */
    //cout << sum[sample-1] << endl << sample;
    //myfile.close();

    // Calculate rms errors of position measurement from Monte-Carlo Simulation
    fout << setprecision(6) << 0 << "," << 0 << "," << 0 << "," << 0 << endl;
    for (int rofl = 0; rofl < sample; rofl++)
    {
        rms[rofl] = sqrt((1./100.)*sum[rofl]);
        haha[rofl] = sqrt((1./100.)*lol[rofl]);
        fout << setprecision(6) << (rofl+1)*dt << "," << rms[rofl] << "," << haha[rofl] << "," << sum[rofl] << endl;
    }
    // Output results to a file
    cout << "filtered: " << rms[sample-1] << endl << "noisy: " << haha[sample-1] << endl;
    fout.close();
    return 0;
}

