/* Homework 7
 * Dormand-Prince 4/5 RK method
 * calculating motion of an light object
 * in gravitational field of two heavy objects
 * 
 * Written by Michael Stumpf
 */
#include <cmath>
#include <fstream>
using namespace std;

//-----------------
//declarations of sub-functions
void func(double y0, double y1, double y2, double y3);
void RK45(double& dt, double* const y, double* const y5);
void stepcontrol(double* const y4, double* const y5, double& dt, double& tol);

//-----------------
//main function
int main(){
  double t0 = 0, tEnd = 18;	//starting and end time
  double dt = 1e-3, tol = 1e-5;	//initial timestep+tol
  const int dim = 4;
  double y4[dim], y5[dim];
  y4[0] = 0.994;		//initial values
  y4[1] = 0;
  y4[2] = 0;
  y4[3] = -2.00158510637908;
  double t = t0;		//set time
  ofstream out("data.txt");	//open output file
  out << t << "\t" << y4[0] << "\t" << y4[1] << "\t" << y4[2] << "\t" << y4[3] << "\t"<< dt << endl;
  while(t<tEnd){		//simulate until tEnd
    RK45(dt, y4, y5);		//solve Dormand-Prince method
    t += dt;
    out << t << "\t" << y4[0] << "\t" << y4[1] << "\t" << y4[2] << "\t" << y4[3] << "\t"<< dt << endl;
    stepcontrol(y4, y5, dt, tol);	//control timesteps
  }
  out.close();
  return 0;
}

//-----------------
//implementation of eode
void func(double y0, double y1, double y2, double y3, double* const k){
  double mu = 0.012277471;
  double r = sqrt((y0+mu)*(y0+mu)+y1*y1);
  double s = sqrt((y0-1+mu)*(y0-1+mu)+y1*y1);
  k[0] = y2;							//x'
  k[1] = y3;							//y'
  k[2] = y0+2*y3-(1-mu)*(y0+mu)/(r*r*r)-mu*(y0-1+mu)/(s*s*s);	//x''
  k[3] = y1-2*y2-(1-mu)*y1/(r*r*r)-mu*y1/(s*s*s);  		//y''
}

//------------------
//implementation of Dormand-Prince method
void RK45(double& dt, double* const y, double* const y5){
  const int dim = 4;
  double k1[dim], k2[dim], k3[dim], k4[dim], k5[dim], k6[dim], k7[dim];
  //calculate all k's first
  func(y[0], y[1], y[2], y[3], k1);
  func(y[0]+1.0/5.0*k1[0]*dt, 
       y[1]+1.0/5.0*k1[1]*dt,
       y[2]+1.0/5.0*k1[2]*dt, 
       y[3]+1.0/5.0*k1[3]*dt, k2);
  func(y[0]+(3.0/40.0*k1[0]+9.0/40.0*k2[0])*dt,
       y[1]+(3.0/40.0*k1[1]+9.0/40.0*k2[1])*dt,
       y[2]+(3.0/40.0*k1[2]+9.0/40.0*k2[2])*dt,
       y[3]+(3.0/40.0*k1[3]+9.0/40.0*k2[3])*dt, k3);
  func(y[0]+(44.0/45.0*k1[0]-56.0/15.0*k2[0]+32.0/9.0*k3[0])*dt,
       y[1]+(44.0/45.0*k1[1]-56.0/15.0*k2[1]+32.0/9.0*k3[1])*dt,
       y[2]+(44.0/45.0*k1[2]-56.0/15.0*k2[2]+32.0/9.0*k3[2])*dt,
       y[3]+(44.0/45.0*k1[3]-56.0/15.0*k2[3]+32.0/9.0*k3[3])*dt, k4);
  func(y[0]+(19372.0/6561.0*k1[0]-25360.0/2187.0*k2[0]+64448.0/6561.0*k3[0]-212.0/729.0*k4[0])*dt,
       y[1]+(19372.0/6561.0*k1[1]-25360.0/2187.0*k2[1]+64448.0/6561.0*k3[1]-212.0/729.0*k4[1])*dt,
       y[2]+(19372.0/6561.0*k1[2]-25360.0/2187.0*k2[2]+64448.0/6561.0*k3[2]-212.0/729.0*k4[2])*dt,
       y[3]+(19372.0/6561.0*k1[3]-25360.0/2187.0*k2[3]+64448.0/6561.0*k3[3]-212.0/729.0*k4[3])*dt, k5);
  func(y[0]+(9017.0/3168.0*k1[0]-355.0/33.0*k2[0]+46732.0/5247.0*k3[0]+49.0/176.0*k4[0]-5103.0/18656.0*k5[0])*dt,
       y[1]+(9017.0/3168.0*k1[1]-355.0/33.0*k2[1]+46732.0/5247.0*k3[1]+49.0/176.0*k4[1]-5103.0/18656.0*k5[1])*dt,
       y[2]+(9017.0/3168.0*k1[2]-355.0/33.0*k2[2]+46732.0/5247.0*k3[2]+49.0/176.0*k4[2]-5103.0/18656.0*k5[2])*dt,
       y[3]+(9017.0/3168.0*k1[3]-355.0/33.0*k2[3]+46732.0/5247.0*k3[3]+49.0/176.0*k4[3]-5103.0/18656.0*k5[3])*dt, k6);
  func(y[0]+(35.0/384.0*k1[0]+500.0/1113.0*k3[0]+125.0/192.0*k4[0]-2187.0/6784.0*k5[0]+11.0/84.0*k6[0])*dt,
       y[1]+(35.0/384.0*k1[1]+500.0/1113.0*k3[1]+125.0/192.0*k4[1]-2187.0/6784.0*k5[1]+11.0/84.0*k6[1])*dt,
       y[2]+(35.0/384.0*k1[2]+500.0/1113.0*k3[2]+125.0/192.0*k4[2]-2187.0/6784.0*k5[2]+11.0/84.0*k6[2])*dt,
       y[3]+(35.0/384.0*k1[3]+500.0/1113.0*k3[3]+125.0/192.0*k4[3]-2187.0/6784.0*k5[3]+11.0/84.0*k6[3])*dt, k7);
  //calculate y's in fifth and fourth order
  for(int i=0; i<dim;  i++){
    y5[i] = y[i]+dt*(
     35.0/384.0		*k1[i]
    +500.0/1113.0	*k3[i]
    +125.0/192.0	*k4[i]
    -2187.0/6784.0	*k5[i]
    +11.0/84.0		*k6[i]);
    y[i] += dt*(
     5179.0/57600.0	*k1[i]
    +7571.0/16695.0	*k3[i]
    +393.0/640.0	*k4[i]
    -92097.0/339200.0	*k5[i]
    +187.0/2100.0	*k6[i]
    +1.0/40.0		*k7[i]);
  }
}

//------------------
//implementation of stepsize control
void stepcontrol(double* const y4, double* const y5, double& dt, double& tol){
  const int dim = 4;
  double v, vmax = 0;
  for(int i=0; i<dim; i++){
    v = abs(y4[i]-y5[i]);	//calculate differences
    if (v>vmax){		//calculate max difference
      vmax = v;
    }
  }
  dt *= pow(tol/vmax, 0.2);	//calculate new stepsize
}