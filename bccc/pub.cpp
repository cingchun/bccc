/**********************************************
         Author:    Qingchun Wang @ NJU
      File Name:    pub.cpp
            Des:    public function
         E-mail:    qingchun720@foxmail.com
    Create Time:    10:11 十月-25/2019
           comp:    icc pub.cpp -std=c++11 -o pub
***********************************************/

#include <cstdio> // printf
#include <cstring>
#include <algorithm>
#include "pub.h"


/* globle variables */
uint nbo=2, nrf=312, nbs=16, nbs2=256, nrs2=79872; // nbs2=nbs*nbs, nrs2 = nrf*nbs2
uint no, *lT0, n2;
uint np, nc=0, *pul, nu; // power of u
dict ud;
real *td;


uint nnhalf(uint no)
{
    return no*(no+1)/2;
}
void lowTri0(uint no)
{
    lT0 = new uint [no]();
    for(uint i=1; i<no; ++i)
        lT0[i]=lT0[i-1]+i;
}
void powul(uint np)
{
    pul = new uint [4]();
    pul[0]=1; uint pu = np*16;
    for (uint i=1; i<4; ++i)
        pul[i] = pul[i-1]*pu;
}

uint index(uint* el, uint n, uint e)
{
    for(uint i=0; i<n; ++i)
        if (el[i]==e) return i;

    return n;
}
uint index(v1D el, uint e)
{
    v1D::iterator i0 = el.begin();

    return find(i0, el.end(), e)-i0;
}
int sign(uint* lnew, uint n, uint* lold)
{
    if (n<=1) return 1; // Note: int 1 for tm in LCC
    else
    {
        uint* il;
        if (!lold) il = lnew;
        else
        {
            il = new uint[n];
            for(uint i=0; i<n; ++i)
                il[i]= index(lnew, n, lold[i]);
        }

        int dif = 1;
        for(uint i=0; i<n; ++i)
            for(uint j=0; j<i; ++j)
                if (il[j]>il[i]) ++dif;
       if (dif%2) return 1;
       else return -1;
    }
}
int sign(v1D lnew, v1D lold)
{
    uint n = lnew.size();

    if (n<=1) return 1; // Note: int 1 for tm in LCC
    else
    {
        v1D il(n);
        if (lold.empty()) il = lnew;
        else
        {
            v1D::iterator ln0=lnew.begin(), ln1=lnew.end();
            for(uint i=0; i<n; ++i)
                il[i]= find(ln0,ln1,lold[i])-ln0;  // distance(find(ln0,ln1,lold[i]),ln0);
        }

        int dif = 1;
        for(uint i=0; i<n; ++i)
            for(uint j=0; j<i; ++j)
                if (il[j]>il[i]) ++dif;
       if (dif%2) return 1;
       else return -1;
    }
}

dict l2d(uint* l, uint n)
{
    dict d;

    for (uint i=0; i<n; ++i) d[l[i]] = i;

    return d;
}
uint phiu(v1D u, uint nc)  // 1D mode
{
    uint iphi=0, nb=u.size(), b, s;

    for (uint i=0,j=0; i<nb; ++j,i=i+nc)
    {
        b=u[i]; s=u[i+1];
        iphi += (b*16+s)*pul[j];
    }

    return iphi;
}
//uint phiu(v2D u)  // 2D mode
//{
//    uint iphi=0, nb=u.size();
//
//    for (uint i=0; i<nb; ++i)
//        iphi += (u[i][0]*16+u[i][1])*pul[i];
//
//    return iphi;
//}

uint** l2ll(uint* l, uint n, uint nc)
{
    uint m=n/nc, inc, **ll = new uint* [m];

    for (uint i=0; i<m; ++i)
    {
        //ll[i] = &l[i*nc];
        inc = i*nc;
        ll[i] = new uint [nc]{l[inc],l[inc+1]}; // nc=2
        //ll[i] = new uint [nc];
        //for (uint j=0; j<nc; ++j)
        //    ll[i][j] = l[inc+j];
    }

    return ll;
}
uint* ll2l(uint** ll, uint m, uint nc)
{
    uint inc, *l = new uint [m*nc];

    for (uint i=0; i<m; ++i)
    {
        inc = i*nc;
        memcpy(&l[inc], ll[i], sizeof(uint)*nc);
        //for (uint j=0; j<nc; ++j)
        //    l[inc+j] = ll[i][j];
    }

    return l;
}
uint** llbind(uint** ll, uint m, uint nc)
{
    return l2ll(ll2l(ll,m*nc),m,nc);
}
uint* v2l(v1D v)
{
    uint *l = &v[0];
    //uint n = v.size(), *l = new uint[n];
    //
    //for (uint i=0; i<n; ++i)
    //    l[i] = v[i];

    return l;
}
v1D l2v(uint* l, uint n)
{
    return v1D(l,l+n);
}
uint** vv2ll(v2D vv)
{
    uint m=vv.size(), nc, **ll = new uint* [m];
    
    for (uint i=0; i<m; ++i)
    {
        ll[i] = &vv[i][0];
        //nc = vv[i].size();
        //ll[i] = new uint[nc];
        //for (uint j=0; j<nc; ++j)
        //    ll[i][j] = vv[i][j];
    }

    return ll;
}
v2D ll2vv(uint** ll, uint m, uint nc)
{
    v2D vv(m,v1D(nc));

    for (uint i = 0; i < m; i++)
        //vv[i] = &ll[i];
        vv[i] = v1D(ll[i], ll[i]+nc);
        //vv[i].assign(ll[i],ll[i]+nc);

    return vv;
}
uint** v2ll(v1D v, uint nc)
{
    return l2ll(v2l(v), v.size(), nc);
}
v1D ll2v(uint** ll, uint m, uint nc)
{
    return l2v(ll2l(ll, m, nc), m*nc);
}
uint* vv2l(v2D vv)
{
    return ll2l(vv2ll(vv), vv.size(), vv[0].size());
}
v2D l2vv(uint* l, uint n, uint nc)
{
    return ll2vv(l2ll(l,n,nc),n/nc,nc);
}
v2D v2vv(v1D v, uint nc)
{
    uint n = v.size(), m = n/nc;

    return ll2vv(l2ll(v2l(v), n, nc), m, nc);
}
v1D vv2v(v2D vv)
{
    uint m=vv.size(), nc=vv[0].size(), mnc=m*nc;

    return l2v(ll2l(vv2ll(vv), m, nc), mnc);
}
v1D usort(v1D u, uint nc)  // 1D mode
{
    uint nb = u.size()/nc;
    if (nb>1)
    {
        uint **u1 = v2ll(u,nc);
        sort(u1, u1+nb, bcmp<uint*>);
        u = ll2v(u1, nb, nc);
    }

    return u;
}
//v2D usort(v2D u)             // 2D mode
//{
//    if (u.size()>1)
//        sort(u.begin(), u.end(), bcmp<v1D>);
//
//    return u;
//}



void test(void)
{
    printf(" Enter test(): \n");
    printf(" nbo,nrf,nbs, nbs2, nrs2  = %d, %d, %d, %d, %d\n", nbo, nrf, nbs, nbs2, nrs2);
    printf(" \n");

    /* test index, sign */
    printf(" test index, sign\n");
    uint ne=4, el[] = {7,9,6,8}, el0[] = {7,9,8,6}, e=el[2];
    v1D l1=l2v(el, ne), l0=l2v(el0, ne);
    printf(" el = [%s]\n", join(el,ne).c_str());
    printf(" %d in el, index is %d \n", e, index(el, ne, e));
    printf(" sign(el) = %d\n", sign(el,ne));
    printf(" sign(l2v(el)) = %d\n", sign(l1));
    printf(" el0 = [%s]\n", join(el0,ne).c_str());
    printf(" sign(el, el0) = %d\n", sign(el,ne,el0));
    printf(" sign(l2v(el), l2v(el0)) = %d\n", sign(l1,l0));
    printf(" \n");

    /* test lT0, powul */
    printf(" test lT0, powul\n");
    no=8; n2=nnhalf(no); lowTri0(no);
    printf(" no,n2, lT0 = %d,%d, [%s]\n", no,n2, join(lT0,no).c_str());
    np=no/nbo; nc=0; powul(np);
    printf(" np,nc, pul = %d,%d, [%s]\n", np, nc, join(pul,np).c_str());
    printf(" \n");

    /* test l2d */
    printf(" test l2d\n");
    nu = 4; uint ul[] = {123, 456, 789, 1010};
    printf(" nu, ul = %d, [%s]\n", nu, join(ul,nu).c_str());
    ud = l2d(ul,nu);
    for (uint i=0; i<nu; ++i) printf(" d[%d] = %d\n", ul[i], ud[ul[i]]);
    printf(" \n"); 

    /* test phiu, dei,sei,rdm, po */
    printf(" test phiu, po, rdm, sei, dei\n");
    printf(" phiu({0,4,1,5}) = %d\n", phiu({ 0,4,1,5 }));
    //printf(" phiu({}) = %d\n", phiu({}));
    printf(" po(1,1) = %d\n", po(1, 1));
    printf(" rdm(1,3,2,0) = %d\n", rdm(1, 3, 2, 0));
    printf(" sei(3,2) = %d\n", sei(3, 2));
    printf(" dei(2,3,2,0) = %d\n", dei(2,3,2,0)); 
    printf(" \n"); 

    /* test ll2l, l2ll, bcmp, usort */
    printf(" test l2ll, ll2l, bcmp, usort\n");
    uint nc=2, nb=3, mn=nb*nc, b0[]={2,3}, b1[]={1,4}, b2[]={1,2}, *u[]={b0, b1, b2};  // C/C++ 只能在一维上直接赋值
    printf(" u = [b0, b1, b2] = [%s]\n", join2(u, nb, nc).c_str());
    uint *l = ll2l(u, nb, nc);
    printf(" l = ll2l(u) = [%s]\n", join(l, mn).c_str());
    v1D u0 = l2v(l, mn);
    printf(" u0 = l2v(l) = [%s]\n", join(v2l(u0), mn).c_str());
    printf(" l2vv(l) = [%s]\n", join(vv2l(l2vv(l,mn,nc)), mn).c_str());
    uint **ll = l2ll(l, mn, nc);
    printf(" ll = l2ll(l) = [%s]\n", join2(ll, nb, nc).c_str());
    printf(" ll2v(ll) = [%s]\n", join2(v2ll(ll2v(ll, nb, nc), nc), nb, nc).c_str());
    v2D u1 = ll2vv(u, nb, nc);
    printf(" u1 = ll2vv(ll) = [%s]\n", join2(vv2ll(u1), nb, nc).c_str());
    printf(" v2vv(u0) = [%s]\n", join2(vv2ll(v2vv(u0, nc)), nb, nc).c_str());
    printf(" vv2v(u1) = [%s]\n", join(v2l(vv2v(u1)), mn).c_str());
    printf(" bcmp(b1, b2) = %d\n", bcmp(b0, b1));
    printf(" usort(u0) = [%s]\n", join(v2l(usort(u0)),mn).c_str());
    //printf(" usort(u1) = [%s]\n", join2(vv2ll(usort(u1)), nb, nc).c_str());
    //printf(" usort({}) = [%s]\n", join2(vv2ll(usort({})), nb, nc).c_str());
    printf(" \n");
}  





//int main(int narg, char* argl[])
//{
//    /* test input */
//    printf(" command line: %s\n", join(argl,narg).c_str());
//    printf(" \n");
//
//
//    test();
//
//
//
//
//    return 0;
//}

