/**********************************************
         Author:    Qingchun Wang @ NJU
      File Name:    pub.h
            Des:    public function
         E-mail:    qingchun720@foxmail.com
    Create Time:    10:11 Ê®ÔÂ-25/2019
***********************************************/


#ifndef PUB_H
#define PUB_H


#include <string>
#include <sstream>
#include <map>
#include <set>
#include <array>
#include <vector>


using namespace std;


typedef unsigned uint;
typedef double real;
typedef map<uint,int> dict;
typedef set<uint> us;
typedef vector<uint> v1D;
typedef vector<v1D> v2D;


/* globle variables */
extern uint nbo, nrf, nbs, nbs2, nrs2; // nbs2=nbs*nbs, nrs2 = nrf*nbs2
extern uint no, *lT0, n2; // n2 = n*(n+1)/2
extern uint np, nc, *pul, nu; // pu = np*nbs
extern dict ud;
extern real *td;


template <typename T>
string join(T l[], uint n, string s=" ")
{
    stringstream buf;

    buf<<l[0];
    for (uint i=1; i<n; ++i)
        buf << s << l[i];

    return buf.str();
}
template <typename T>
string join2(T* ll[], uint m, uint n, string s=" ")
{
    stringstream buf;

    buf << join(ll[0], n, s);
    for (uint i=1; i<m; ++i)
        buf << "; " << join(ll[i], n, s);

    return buf.str();
}

inline uint po(uint p, uint o)
{
    return o==0 ? p : (p>=nc ? p+(np-p)*2-1 : 0);
}
inline uint sei(uint i, uint j)
{
    return i*no+j;
}
void lowTri0(uint no);
uint nnhalf(uint no);
inline uint dei(uint i, uint j, uint k, uint l)
{
    uint ij = i<=j ? lT0[j]+i : lT0[i]+j;
    uint kl = k<=l ? lT0[l]+k : lT0[k]+l;

    return ij*n2+kl;
}

void powul(uint np);
inline uint rdm(uint p, uint r, uint i, uint j)
{
    return p*79872 + r*256 + i*16 + j;
}
dict l2d(uint* l, uint n);
uint phiu(v1D u, uint nc=2); // 1D mode
//uint phiu(v2D u);            // 2D mode

uint index(uint* el, uint n, uint e);
uint index(v1D el, uint e);
int sign(uint* lnew, uint n, uint* lold=NULL);
int sign(v1D lnew, v1D lold={});

uint** l2ll(uint* l, uint n, uint nc=2);
uint* ll2l(uint** ll, uint m, uint nc=2);
uint** llbind(uint **ll, uint m, uint nc=2);
uint* v2l(v1D v);
v1D l2v(uint* l, uint n);
uint** vv2ll(v2D vv);
v2D ll2vv(uint** ll, uint m, uint nc=2);
uint** v2ll(v1D v, uint nc=2);
v1D ll2v(uint** ll, uint m, uint nc=2);
uint* vv2l(v2D vv);
v2D l2vv(uint* l, uint n, uint nc=2);
v2D v2vv(v1D v, uint nc=2);
v1D vv2v(v2D vv);
template <typename T>
bool bcmp(T b1, T b2)
{
    if (b1[0]!=b2[0]) return b1[0]<b2[0];
    else return b1[1]<=b2[1];
}
v1D usort(v1D u, uint nc=2); // 1D mode
//v2D usort(v2D);              // 2D mode




void test(void);





#endif
