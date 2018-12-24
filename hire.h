#include "in.h"
#ifndef HIRE_H
#define HIRE_H


class Hire
{
public:
    void interModuleCC(int ll, int nn, iMatrix& aa);
    void intraModuleCHMN1(int l, int nn, iMatrix &aa,double prob,double alpha);
    void interModuleCCT2(int ll, int nn, iMatrix& aa, double avgcon);
    void intraModuleRC(int ll, int nn, iMatrix& aa);
    void intraModuleRCT2(int ll, int nn, iMatrix& aa, int outcon);
    void EinterModuleCC(int ll, int nn, iMatrix& aa);
    void EintraModuleRC(int ll, int nn, iMatrix& aa);
    void intraModuleN0(int ll, int nn, iMatrix& aa,int bn);
    void intraModuleN1(int ll, int nn, iMatrix& aa,int bn);
    void intraModuleFfibo(int ll, int nn, iMatrix &aa);
    void intraModuleFfiboV(int ll, int nn, iMatrix &aa);
    int  NM1Ccheck(int m0,int c1,int c2,iMatrix& aa);
    //----
    void interModuleNM1(int ll, int nn, iMatrix& aa);
    //----
    void intraModuleNM1fibo(int ll, int nn, iMatrix& aa);
    void intraModuleNM1RC(int ll, int nn, iMatrix& aa, int avgnl);
    void intraModuleNM1_MHMN(int ll, int nn, iMatrix& aa,double s,double alpha);
    void intraModuleNM1_SR_MHMN(int ll, int nn, iMatrix& aa,double s,double alpha);
    void intraModuleNM2_MHMN(int ll, int nn, iMatrix& aa,double s,double alpha);
    //----
    void intraModuleHMN1(int l, int nn, iMatrix& aa, double prob, double alpha);
    void intraModuleHMN2(int l, int nn, iMatrix& aa, double prob, double alpha);
    void intraModuleHMN3(int l, int nn, iMatrix& aa, double prob, double alpha);
    //----
    void Hanoi01(int m0, int l, int nn, iMatrix& aa);
    //
    void intraModuleHMNfibolink(int l, int nn, iMatrix &aa,double prob,double alpha);
};

#endif // HIRE_H
