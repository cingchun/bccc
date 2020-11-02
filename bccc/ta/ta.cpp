// Author: Qingchun Wang @ NJU
// E-mail: qingchun720@foxmail.com 


#include "ta.h"


real* ta(real* t, uint* pul, uint* ul, uint nu, real* r, uint np, real* h, uint no, real* g, uint* lT0, uint n2, uint* p, uint nc)
{
     //test();
     ::no = no; ::lT0 = lT0; ::n2 = n2;
     ::np = np; ::nc = nc;
     ::pul=pul; ud = l2d(ul, nu); ::nu = nu;
     td = new real[nu]();
     
     ta_(t, r, h, g, p, nc);
     ta_A1B1(t, r, h, g, p, nc);
     ta_A1B2(t, r, h, g, p, nc);
     ta_A2B2(t, r, h, g, p, nc);
     ta_A3B3(t, r, h, g, p, nc);
     ta_A6B15(t, r, h, g, p, nc);
     ta_A7B13(t, r, h, g, p, nc);
     ta_A7B14(t, r, h, g, p, nc);
     ta_A8B13(t, r, h, g, p, nc);
     ta_A8B14(t, r, h, g, p, nc);
     
     return td;
}






int main(int narg, char* argl[])
{
     /* test input */
     printf(" command line: %s\n", join(argl, narg).c_str());
     printf(" \n");
     
     test();
     

}


