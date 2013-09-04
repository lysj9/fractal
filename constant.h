#ifndef _CONSTANT_H_
#define _CONSTANT_H_

#define GRAV 6.67384e-11
#define MSUN 1.98855e30
#define PC   3.08567758e16
#define MYR  31557600e6
#define VSC  65.58143      // velocity scale coefficient = sqrt(G*MSUN/pc) [m/s]
//#define VSC (sqrt(GRAV*MSUN/PC))
#define TSC  14.90959      // time scale coefficient = sqrt(pc^3/(G*MSUN)) [Myr]
//#define TSC  (sqrt(PC*PC*PC/(GRAV*MSUN))/MYR)

#define R_SUN 4.6491e-3
#define PI 3.14159265358979323846
#define TWO_PI (2*PI)

#endif
