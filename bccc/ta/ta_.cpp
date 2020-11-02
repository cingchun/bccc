// Author: Qingchun Wang @ NJU 
// E-mail: qingchun720@foxmail.com 


#include "../pub.h"
#include "ta.h"


void ta_(real* t, real* r, real* h, real* g, uint* p, uint nc)
{   
    
    uint u = ud[phiu({})];
    real tu = 0.0, v=0.0;
    us l = {};
    
    for (uint I=nc; I<np; ++I) {
        if (l.find(I) != l.end()) continue;
        us lI = l; lI.insert(I);
        tu += h[sei(po(I,0),po(I,0))]*r[rdm(I,9,1,0)]*t[ud[phiu({I,1})]];
        tu += h[sei(po(I,0),po(I,0))]*r[rdm(I,24,1,0)]*t[ud[phiu({I,1})]];
        tu += h[sei(po(I,1),po(I,1))]*r[rdm(I,39,1,0)]*t[ud[phiu({I,1})]];
        tu += h[sei(po(I,1),po(I,1))]*r[rdm(I,54,1,0)]*t[ud[phiu({I,1})]];
        tu += g[dei(po(I,0),po(I,0),po(I,0),po(I,0))]*r[rdm(I,194,1,0)]*t[ud[phiu({I,1})]];
        tu += g[dei(po(I,0),po(I,1),po(I,0),po(I,1))]*r[rdm(I,199,1,0)]*t[ud[phiu({I,1})]];
        tu += g[dei(po(I,0),po(I,1),po(I,0),po(I,1))]*r[rdm(I,274,1,0)]*t[ud[phiu({I,1})]];
        tu += g[dei(po(I,1),po(I,1),po(I,1),po(I,1))]*r[rdm(I,279,1,0)]*t[ud[phiu({I,1})]];
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,1})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,1})]];
        tu += v; v = 0.0;
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,15,2,0)]*t[ud[phiu({I,2})]];
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,30,2,0)]*t[ud[phiu({I,2})]];
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,33,2,0)]*t[ud[phiu({I,2})]];
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,48,2,0)]*t[ud[phiu({I,2})]];
        tu += g[dei(po(I,0),po(I,0),po(I,0),po(I,1))]*r[rdm(I,210,2,0)]*t[ud[phiu({I,2})]];
        tu += g[dei(po(I,0),po(I,1),po(I,1),po(I,1))]*r[rdm(I,215,2,0)]*t[ud[phiu({I,2})]];
        tu += g[dei(po(I,0),po(I,0),po(I,0),po(I,1))]*r[rdm(I,258,2,0)]*t[ud[phiu({I,2})]];
        tu += g[dei(po(I,0),po(I,1),po(I,1),po(I,1))]*r[rdm(I,263,2,0)]*t[ud[phiu({I,2})]];
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,2})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,2})]];
        tu += v; v = 0.0;
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,15,3,0)]*t[ud[phiu({I,3})]];
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,30,3,0)]*t[ud[phiu({I,3})]];
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,33,3,0)]*t[ud[phiu({I,3})]];
        tu += h[sei(po(I,0),po(I,1))]*r[rdm(I,48,3,0)]*t[ud[phiu({I,3})]];
        tu += g[dei(po(I,0),po(I,0),po(I,0),po(I,1))]*r[rdm(I,210,3,0)]*t[ud[phiu({I,3})]];
        tu += g[dei(po(I,0),po(I,1),po(I,1),po(I,1))]*r[rdm(I,215,3,0)]*t[ud[phiu({I,3})]];
        tu += g[dei(po(I,0),po(I,0),po(I,0),po(I,1))]*r[rdm(I,258,3,0)]*t[ud[phiu({I,3})]];
        tu += g[dei(po(I,0),po(I,1),po(I,1),po(I,1))]*r[rdm(I,263,3,0)]*t[ud[phiu({I,3})]];
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,9,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,39,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,3})]];
        tu += -v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,24,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        for (uint J=0; J<np; ++J)
            if (lI.find(J)==lI.end())
                v += g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,54,0,0)]*t[ud[phiu({I,3})]];
        tu += v; v = 0.0;
        }
        
    for (uint I=nc; I<np; ++I) {
        if (l.find(I) != l.end()) continue;
        us lI = l; lI.insert(I);
        for (uint J=nc; J<np; ++J) {
            if (lI.find(J) != lI.end()) continue;
            us lIJ = lI; lIJ.insert(J);
            tu += 0.5*h[sei(po(I,0),po(J,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,76,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += -0.25*g[dei(po(I,0),po(J,0),po(I,1),po(I,1))]*r[rdm(I,124,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.25*g[dei(po(I,0),po(J,0),po(I,1),po(I,1))]*r[rdm(I,76,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*g[dei(po(I,0),po(J,0),po(I,1),po(I,1))]*r[rdm(I,86,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,124,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,146,12,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*g[dei(po(I,0),po(J,0),po(J,0),po(J,0))]*r[rdm(I,0,12,0)]*r[rdm(J,97,9,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,117,9,0)]*t[ud[phiu({I,12,J,9})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,12,J,9})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,12,J,9})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,12,J,9})]];
            tu += 0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,0),po(J,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,118,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.25*g[dei(po(I,0),po(J,0),po(I,1),po(I,1))]*r[rdm(I,166,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(I,1))]*r[rdm(I,132,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,144,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += 0.25*g[dei(po(I,0),po(J,0),po(I,1),po(I,1))]*r[rdm(I,118,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,166,14,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(J,0),po(J,0))]*r[rdm(I,2,14,0)]*r[rdm(J,65,7,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,85,7,0)]*t[ud[phiu({I,14,J,7})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,14,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,14,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,14,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,14,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,14,J,7})]];
            tu += -0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,76,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += -0.25*g[dei(po(I,0),po(J,1),po(I,1),po(I,1))]*r[rdm(I,124,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.25*g[dei(po(I,0),po(J,1),po(I,1),po(I,1))]*r[rdm(I,76,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(I,1),po(I,1))]*r[rdm(I,86,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,124,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,146,12,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,161,10,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(J,1),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,181,10,0)]*t[ud[phiu({I,12,J,10})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,12,J,10})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,12,J,10})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,12,J,10})]];
            tu += 0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,118,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.25*g[dei(po(I,0),po(J,1),po(I,1),po(I,1))]*r[rdm(I,166,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(I,1))]*r[rdm(I,132,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,144,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += 0.25*g[dei(po(I,0),po(J,1),po(I,1),po(I,1))]*r[rdm(I,118,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,166,14,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,129,8,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(J,1),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,149,8,0)]*t[ud[phiu({I,14,J,8})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,14,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,14,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,14,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,14,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,14,J,8})]];
            tu += -0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,1),po(J,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += -0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,0))]*r[rdm(I,72,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,120,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,70,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,72,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,0))]*r[rdm(I,120,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(I,1),po(J,0))]*r[rdm(I,130,11,0)]*r[rdm(J,1,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*g[dei(po(I,1),po(J,0),po(J,0),po(J,0))]*r[rdm(I,4,11,0)]*r[rdm(J,97,9,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,117,9,0)]*t[ud[phiu({I,11,J,9})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,11,J,9})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,11,J,9})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,11,J,9})]];
            tu += 0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,1),po(J,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,68,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,0))]*r[rdm(I,114,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.5*g[dei(po(I,0),po(I,0),po(I,1),po(J,0))]*r[rdm(I,80,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,162,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,114,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += 0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,0))]*r[rdm(I,162,13,0)]*r[rdm(J,3,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(J,0),po(J,0))]*r[rdm(I,6,13,0)]*r[rdm(J,65,7,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,85,7,0)]*t[ud[phiu({I,13,J,7})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,13,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,13,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,13,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,13,J,7})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,13,J,7})]];
            tu += -0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,1),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += -0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,1))]*r[rdm(I,72,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,120,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,70,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,72,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,1))]*r[rdm(I,120,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(I,1),po(J,1))]*r[rdm(I,130,11,0)]*r[rdm(J,5,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,161,10,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*g[dei(po(I,1),po(J,1),po(J,1),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,181,10,0)]*t[ud[phiu({I,11,J,10})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,11,J,10})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,11,J,10})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,11,J,10})]];
            tu += 0.5*v; v = 0.0;
            tu += 0.5*h[sei(po(I,1),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,68,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,1))]*r[rdm(I,114,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.5*g[dei(po(I,0),po(I,0),po(I,1),po(J,1))]*r[rdm(I,80,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,162,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += 0.25*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,114,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += 0.25*g[dei(po(I,0),po(I,0),po(I,1),po(J,1))]*r[rdm(I,162,13,0)]*r[rdm(J,7,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,129,8,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(J,1),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,149,8,0)]*t[ud[phiu({I,13,J,8})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,13,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,13,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,13,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,13,J,8})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,13,J,8})]];
            tu += -0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*g[dei(po(I,0),po(I,0),po(I,0),po(J,0))]*r[rdm(I,97,9,0)]*r[rdm(J,0,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,117,9,0)]*r[rdm(J,0,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.25*g[dei(po(I,0),po(J,0),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,76,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,86,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,124,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,146,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += 0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,76,12,0)]*t[ud[phiu({I,9,J,12})]];
            tu += 0.25*g[dei(po(I,0),po(J,0),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,124,12,0)]*t[ud[phiu({I,9,J,12})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,9,J,12})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,9,J,12})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,9,J,12})]];
            tu += -0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(I,0),po(J,0))]*r[rdm(I,65,7,0)]*r[rdm(J,2,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,0))]*r[rdm(I,85,7,0)]*r[rdm(J,2,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += -0.25*g[dei(po(I,0),po(J,0),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,118,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += -0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,166,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,118,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.25*g[dei(po(I,0),po(J,0),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,166,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.5*g[dei(po(I,0),po(J,0),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,132,14,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,144,14,0)]*t[ud[phiu({I,7,J,14})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,7,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,7,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,7,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,7,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,7,J,14})]];
            tu += 0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,1),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,161,10,0)]*r[rdm(J,0,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*g[dei(po(I,1),po(I,1),po(I,1),po(J,0))]*r[rdm(I,181,10,0)]*r[rdm(J,0,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.25*g[dei(po(I,1),po(J,0),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,76,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,86,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,124,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,146,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += 0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,76,12,0)]*t[ud[phiu({I,10,J,12})]];
            tu += 0.25*g[dei(po(I,1),po(J,0),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,124,12,0)]*t[ud[phiu({I,10,J,12})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,10,J,12})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,10,J,12})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,10,J,12})]];
            tu += -0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,1),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,0))]*r[rdm(I,129,8,0)]*r[rdm(J,2,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(I,1),po(J,0))]*r[rdm(I,149,8,0)]*r[rdm(J,2,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += -0.25*g[dei(po(I,1),po(J,0),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,118,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += -0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,166,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,118,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.25*g[dei(po(I,1),po(J,0),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,166,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.5*g[dei(po(I,1),po(J,0),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,132,14,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.5*g[dei(po(I,1),po(J,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,144,14,0)]*t[ud[phiu({I,8,J,14})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,8,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,8,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,8,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,8,J,14})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,8,J,14})]];
            tu += 0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*g[dei(po(I,0),po(I,0),po(I,0),po(J,1))]*r[rdm(I,97,9,0)]*r[rdm(J,4,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,117,9,0)]*r[rdm(J,4,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,70,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.25*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,72,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,120,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,130,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += 0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,72,11,0)]*t[ud[phiu({I,9,J,11})]];
            tu += 0.25*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,120,11,0)]*t[ud[phiu({I,9,J,11})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,9,J,11})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,9,J,11})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,9,J,11})]];
            tu += -0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(I,0),po(J,1))]*r[rdm(I,65,7,0)]*r[rdm(J,6,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,1),po(J,1))]*r[rdm(I,85,7,0)]*r[rdm(J,6,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += -0.25*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,114,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += -0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,162,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.5*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,68,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.25*g[dei(po(I,0),po(J,1),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,114,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,80,13,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.25*g[dei(po(I,0),po(J,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,162,13,0)]*t[ud[phiu({I,7,J,13})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,7,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,7,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,7,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,7,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,7,J,13})]];
            tu += 0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,161,10,0)]*r[rdm(J,4,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*g[dei(po(I,1),po(I,1),po(I,1),po(J,1))]*r[rdm(I,181,10,0)]*r[rdm(J,4,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,70,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.25*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,72,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,120,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,130,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += 0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,72,11,0)]*t[ud[phiu({I,10,J,11})]];
            tu += 0.25*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,120,11,0)]*t[ud[phiu({I,10,J,11})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,10,J,11})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,10,J,11})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,10,J,11})]];
            tu += -0.5*v; v = 0.0;
            tu += -0.5*h[sei(po(I,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(I,0),po(J,1))]*r[rdm(I,129,8,0)]*r[rdm(J,6,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(I,1),po(J,1))]*r[rdm(I,149,8,0)]*r[rdm(J,6,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += -0.25*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,114,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += -0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,162,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.5*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,68,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.25*g[dei(po(I,1),po(J,1),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,114,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.5*g[dei(po(I,1),po(J,1),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,80,13,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.25*g[dei(po(I,1),po(J,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,162,13,0)]*t[ud[phiu({I,8,J,13})]];
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,8,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,9,0,0)]*t[ud[phiu({I,8,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,8,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,39,0,0)]*t[ud[phiu({I,8,J,13})]];
            tu += -0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,0,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.5*v; v = 0.0;
            for (uint K=0; K<np; ++K)
                if (lIJ.find(K)==lIJ.end())
                    v += g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,0,0)]*t[ud[phiu({I,8,J,13})]];
            tu += 0.5*v; v = 0.0;
            tu += 0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,11,15,0)]*r[rdm(J,22,6,0)]*t[ud[phiu({I,15,J,6})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,11,15,0)]*r[rdm(J,52,6,0)]*t[ud[phiu({I,15,J,6})]];
            tu += 0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,41,15,0)]*r[rdm(J,22,6,0)]*t[ud[phiu({I,15,J,6})]];
            tu += 0.5*g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,41,15,0)]*r[rdm(J,52,6,0)]*t[ud[phiu({I,15,J,6})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,9,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,39,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,0))]*r[rdm(I,24,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,0))]*r[rdm(I,54,1,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,1),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,9,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,39,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,0),po(J,0),po(J,1))]*r[rdm(I,24,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,1),po(I,1),po(J,0),po(J,1))]*r[rdm(I,54,1,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3,J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,24,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,54,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,9,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,39,1,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,1})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3,J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,48,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,15,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,33,2,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,2})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3,J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,30,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,48,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,15,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += 0.5*g[dei(po(I,0),po(I,1),po(J,0),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,33,3,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,3})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,18,4,0)]*r[rdm(J,27,5,0)]*t[ud[phiu({I,4,J,5})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,18,4,0)]*r[rdm(J,45,5,0)]*t[ud[phiu({I,4,J,5})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,36,4,0)]*r[rdm(J,27,5,0)]*t[ud[phiu({I,4,J,5})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,36,4,0)]*r[rdm(J,45,5,0)]*t[ud[phiu({I,4,J,5})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,27,5,0)]*r[rdm(J,18,4,0)]*t[ud[phiu({I,5,J,4})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,45,5,0)]*r[rdm(J,18,4,0)]*t[ud[phiu({I,5,J,4})]];
            tu += -0.5*g[dei(po(I,0),po(J,0),po(I,1),po(J,1))]*r[rdm(I,27,5,0)]*r[rdm(J,36,4,0)]*t[ud[phiu({I,5,J,4})]];
            tu += -0.5*g[dei(po(I,0),po(J,1),po(I,1),po(J,0))]*r[rdm(I,45,5,0)]*r[rdm(J,36,4,0)]*t[ud[phiu({I,5,J,4})]];
            tu += 0.5*g[dei(po(I,0),po(J,0),po(I,0),po(J,0))]*r[rdm(I,22,6,0)]*r[rdm(J,11,15,0)]*t[ud[phiu({I,6,J,15})]];
            tu += 0.5*g[dei(po(I,1),po(J,0),po(I,1),po(J,0))]*r[rdm(I,52,6,0)]*r[rdm(J,11,15,0)]*t[ud[phiu({I,6,J,15})]];
            tu += 0.5*g[dei(po(I,0),po(J,1),po(I,0),po(J,1))]*r[rdm(I,22,6,0)]*r[rdm(J,41,15,0)]*t[ud[phiu({I,6,J,15})]];
            tu += 0.5*g[dei(po(I,1),po(J,1),po(I,1),po(J,1))]*r[rdm(I,52,6,0)]*r[rdm(J,41,15,0)]*t[ud[phiu({I,6,J,15})]];
            }}
            
    for (uint I=nc; I<np; ++I) {
        if (l.find(I) != l.end()) continue;
        us lI = l; lI.insert(I);
        for (uint J=nc; J<np; ++J) {
            if (lI.find(J) != lI.end()) continue;
            us lIJ = lI; lIJ.insert(J);
            for (uint K=nc; K<np; ++K) {
                if (lIJ.find(K) != lIJ.end()) continue;
                us lIJK = lIJ; lIJK.insert(K);
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,15,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,15,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,30,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,30,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,15,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,15,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,30,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,30,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,9,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,39,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,0,12,0)]*r[rdm(J,9,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,0,12,0)]*r[rdm(J,24,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,39,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,54,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,24,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,54,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,2,14,0)]*r[rdm(J,24,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,54,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,2,14,0)]*r[rdm(J,9,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,39,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,9,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,39,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,0,12,0)]*r[rdm(J,9,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,0,12,0)]*r[rdm(J,24,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,39,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,54,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,24,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,54,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,2,14,0)]*r[rdm(J,24,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,54,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,2,14,0)]*r[rdm(J,9,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,39,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,15,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,33,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,15,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,30,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,33,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,48,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,15,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,33,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,15,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,30,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,33,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,48,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,30,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,48,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,30,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,48,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,15,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,33,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,30,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,48,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,30,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,48,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,15,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,33,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,15,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,33,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,15,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,30,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,33,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,48,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,15,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,33,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,15,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,30,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,33,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,0,12,0)]*r[rdm(J,48,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,12,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,30,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,48,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,30,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,48,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,15,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,33,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,30,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,48,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,30,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,48,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,15,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,2,14,0)]*r[rdm(J,33,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,14,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,9,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,39,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,4,11,0)]*r[rdm(J,9,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,4,11,0)]*r[rdm(J,24,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,39,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,54,1,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,24,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,54,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,6,13,0)]*r[rdm(J,24,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,54,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,6,13,0)]*r[rdm(J,9,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,39,1,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,9,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,39,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,4,11,0)]*r[rdm(J,9,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,4,11,0)]*r[rdm(J,24,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,39,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,54,1,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,24,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,54,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,6,13,0)]*r[rdm(J,24,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,54,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,6,13,0)]*r[rdm(J,9,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,39,1,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,15,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,33,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,15,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,30,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,33,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,48,2,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,15,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,33,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,15,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,30,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,33,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,48,3,0)]*r[rdm(K,1,9,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,30,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,48,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,30,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,48,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,15,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,33,2,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,30,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,48,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,30,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,48,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,15,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,33,3,0)]*r[rdm(K,3,7,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,15,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,33,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,15,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,30,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,33,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,48,2,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,15,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,33,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,15,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,30,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,33,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,4,11,0)]*r[rdm(J,48,3,0)]*r[rdm(K,5,10,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,11,K,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,30,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,48,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,30,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,48,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,15,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,33,2,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,30,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,48,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,30,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,48,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,15,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,6,13,0)]*r[rdm(J,33,3,0)]*r[rdm(K,7,8,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,13,K,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,0))]*r[rdm(I,9,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,0))]*r[rdm(I,39,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,0))]*r[rdm(I,24,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,0))]*r[rdm(I,54,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,30,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,15,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,48,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,0))]*r[rdm(I,33,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,0))]*r[rdm(I,30,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,0),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(I,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,0),po(J,1),po(K,1))]*r[rdm(I,9,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(I,1),po(J,1),po(K,1))]*r[rdm(I,39,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,0),po(K,1))]*r[rdm(I,24,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(I,1),po(K,1))]*r[rdm(I,54,1,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,1})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,33,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,33,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,48,2,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,0),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(I,1),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,0))]*r[rdm(I,48,3,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,33,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,33,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,30,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,48,2,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,2})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,15,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,48,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(I,1),po(J,1),po(K,1))]*r[rdm(I,33,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(I,1),po(K,1))]*r[rdm(I,30,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(I,1),po(J,1))]*r[rdm(I,48,3,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({I,3})]]*t[ud[phiu({J,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,12,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,14,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,9})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,9})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,7})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,7})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,10})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,11,J,10})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,8})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,13,J,8})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,9,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,39,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,9,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,39,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,24,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,54,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,24,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,54,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,24,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,9,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,54,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,39,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,15,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,33,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,15,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,33,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,30,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,48,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,15,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,33,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,15,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,33,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,30,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,48,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,30,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,48,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,30,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,15,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,48,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,33,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,30,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,48,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,30,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,15,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,48,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,33,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,9,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,39,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,9,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,39,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,24,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,54,1,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,24,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,54,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,24,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,9,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,54,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,39,1,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,15,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,33,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,15,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,33,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,30,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,48,2,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,15,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,33,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,15,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,33,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,30,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,48,3,0)]*r[rdm(K,0,12,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,30,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,48,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,30,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,15,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,48,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,33,2,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,30,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,48,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,30,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,15,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,48,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,33,3,0)]*r[rdm(K,2,14,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,9,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,39,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,9,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,39,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,1,9,0)]*r[rdm(J,24,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,54,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,24,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,54,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,24,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,0))]*r[rdm(I,3,7,0)]*r[rdm(J,9,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,54,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,39,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,15,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,33,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,15,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,33,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,30,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,48,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,15,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,33,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,15,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,33,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,30,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,1,9,0)]*r[rdm(J,48,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,9,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,30,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,48,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,30,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,15,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,48,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,33,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,30,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(J,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,48,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,30,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,15,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,48,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(J,1))]*r[rdm(I,3,7,0)]*r[rdm(J,33,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,7,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,9,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,39,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,9,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,39,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,5,10,0)]*r[rdm(J,24,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,54,1,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,24,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,54,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,24,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,0))]*r[rdm(I,7,8,0)]*r[rdm(J,9,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,54,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,39,1,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,1})]]*t[ud[phiu({I,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,15,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,33,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,15,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,33,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,30,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,48,2,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,15,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,33,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,15,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,33,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,30,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,5,10,0)]*r[rdm(J,48,3,0)]*r[rdm(K,4,11,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,10,K,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,30,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,48,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,30,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,15,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,48,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,33,2,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,2})]]*t[ud[phiu({I,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,30,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(J,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,48,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,30,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,15,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,48,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(J,1))]*r[rdm(I,7,8,0)]*r[rdm(J,33,3,0)]*r[rdm(K,6,13,0)]*t[ud[phiu({J,3})]]*t[ud[phiu({I,8,K,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,0),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,12})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,12})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,14})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,0),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,14})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,9,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,0),po(K,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,0),po(J,1),po(K,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,7,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,24,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,54,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,9,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,39,1,0)]*t[ud[phiu({K,1})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,11})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,10,J,11})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,15,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,33,2,0)]*t[ud[phiu({K,2})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,30,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,13})]];
                tu += 0.16666666666666666*g[dei(po(I,1),po(K,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,48,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,15,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,13})]];
                tu += -0.16666666666666666*g[dei(po(I,1),po(J,1),po(K,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,33,3,0)]*t[ud[phiu({K,3})]]*t[ud[phiu({I,8,J,13})]];
                }}}
                
    for (uint I=nc; I<np; ++I) {
        if (l.find(I) != l.end()) continue;
        us lI = l; lI.insert(I);
        for (uint J=nc; J<np; ++J) {
            if (lI.find(J) != lI.end()) continue;
            us lIJ = lI; lIJ.insert(J);
            for (uint K=nc; K<np; ++K) {
                if (lIJ.find(K) != lIJ.end()) continue;
                us lIJK = lIJ; lIJK.insert(K);
                for (uint L=nc; L<np; ++L) {
                    if (lIJK.find(L) != lIJK.end()) continue;
                    us lIJKL = lIJK; lIJKL.insert(L);
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,12,L,9})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,12,L,10})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,7,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,8,K,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,11,L,9})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,7,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,11,L,10})]]*t[ud[phiu({J,8,K,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,J,9})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,J,10})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,J,9})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,J,10})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,K,9})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,0,12,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,12,K,10})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,K,9})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,4,11,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,11,K,10})]]*t[ud[phiu({J,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,14,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,14,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,13,L,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,13,L,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,2,14,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,6,13,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,J,12})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,12,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,14,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,12,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,14,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,0,12,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,J,12})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,J,11})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,11,K,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,13,K,7})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,3,7,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,11,K,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,4,11,0)]*r[rdm(K,7,8,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,J,11})]]*t[ud[phiu({K,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,12,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,0,12,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,12,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,1,9,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,11,K,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,4,11,0)]*r[rdm(K,5,10,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,11,K,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,9,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,10,K,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,14,L,7})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,14,L,8})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,9,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,13,L,7})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,13,L,8})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,12,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,12,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,11,L,9})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,14,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,14,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,3,7,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,13,L,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,7,8,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,13,L,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,12,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,0,12,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,12,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,1,9,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,11,L,9})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,11,L,10})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,4,11,0)]*r[rdm(L,5,10,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,11,L,10})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,K,12})]]*t[ud[phiu({J,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,K,12})]]*t[ud[phiu({J,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,7,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,7,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,8,L,14})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,8,L,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,9,K,11})]]*t[ud[phiu({J,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,7,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,7,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,8,L,13})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,10,K,11})]]*t[ud[phiu({J,8,L,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,9,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,0,12,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,10,K,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,7,L,14})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,2,14,0)]*t[ud[phiu({I,8,L,14})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,9,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,7,L,13})]]*t[ud[phiu({J,10,K,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,4,11,0)]*r[rdm(L,6,13,0)]*t[ud[phiu({I,8,L,13})]]*t[ud[phiu({J,10,K,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,K,7})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,K,8})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,1,9,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,K,7})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,5,10,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,K,8})]]*t[ud[phiu({J,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,9,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,10,L,12})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,J,7})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,2,14,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,14,J,8})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,9,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,3,7,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,J,7})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,6,13,0)]*r[rdm(J,7,8,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,13,J,8})]]*t[ud[phiu({K,10,L,11})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,14,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,2,14,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,14,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,13,K,8})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,3,7,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,13,K,7})]];
                    tu += -0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,6,13,0)]*r[rdm(K,7,8,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,13,K,8})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,0),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,J,14})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,0),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,2,14,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,J,14})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,1,9,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(J,1),po(K,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,J,13})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(J,1),po(K,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,6,13,0)]*r[rdm(K,5,10,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,J,13})]]*t[ud[phiu({K,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,0),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,9,L,12})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,0),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,10,L,12})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,0))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,7,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,0))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,8,K,14})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,0),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(L,1),po(J,1),po(K,1))]*r[rdm(I,1,9,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,9,L,11})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,0),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,3,7,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,7,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(L,1),po(J,1),po(K,1))]*r[rdm(I,5,10,0)]*r[rdm(J,7,8,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,10,L,11})]]*t[ud[phiu({J,8,K,13})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,9,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,0))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,0))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,0,12,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,10,L,12})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,0),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,K,14})]]*t[ud[phiu({J,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,0),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,2,14,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,K,14})]]*t[ud[phiu({J,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,0),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,0),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,1,9,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,9,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,0),po(K,1),po(J,1),po(L,1))]*r[rdm(I,3,7,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,7,K,13})]]*t[ud[phiu({J,10,L,11})]];
                    tu += 0.041666666666666664*g[dei(po(I,1),po(K,1),po(J,1),po(L,1))]*r[rdm(I,7,8,0)]*r[rdm(J,5,10,0)]*r[rdm(K,6,13,0)]*r[rdm(L,4,11,0)]*t[ud[phiu({I,8,K,13})]]*t[ud[phiu({J,10,L,11})]];
                    }}}}
                    
    td[u] = tu;


}   
    

