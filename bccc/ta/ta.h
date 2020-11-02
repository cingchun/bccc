// Author: Qingchun Wang @ NJU
// E-mail: qingchun720@foxmail.com 


#ifndef TA_H
#define TA_H


#include "../pub.h"
#include "omp.h"


extern "C" {

real* ta(real* t, uint* pul, uint* ul, uint nu, real* r, uint np, real* h, uint no, real* g, uint* lT0, uint n2, uint* p, uint nc);

}


void ta_(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A1B1(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A1B2(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A2B2(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A3B3(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A6B15(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A7B13(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A7B14(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A8B13(real* t, real* r, real* h, real* g, uint* p, uint nc=0);
void ta_A8B14(real* t, real* r, real* h, real* g, uint* p, uint nc=0);






#endif


