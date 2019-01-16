#include "in.h"
#include "global_variables.h"
#ifndef IZHIRK_H
#define IZHIRK_H


class IzhiRK
{
private:
    int state;//chek parameter for initialize
    int Spike,myNum,spikeNumber,spikeCon;
    double v,u,I,rI,D,cashI,phi;//cashI is sum of input
    //std::vector<int> inC;
    //std::vector<double> inW;
    std::vector<double> spikeTime;
    std::vector<int> outC;
    //std::vector<double> outW;    
public:
    IzhiRK();
    double fv(double v,double u,double I);
    double fu(double v,double u);
    double k1v(double v,double u,double I,double h);
    double k1u(double v,double u,double h);
    double k2v(double v,double u,double I,double h);
    double k2u(double v,double u,double I,double h);
    double k3v(double v,double u,double I,double h);
    double k3u(double v,double u,double I,double h);
    double k4v(double v,double u,double I,double h);
    double k4u(double v,double u,double I,double h);
    double RKv(double v,double u,double I);
    double RKu(double v, double u,double I);
    void Initialize(double ri);
    void Connect(std::vector<int>  & cv);
    void UpdateI();
    void Run(double step_time);
    void get_spikeTime(vector<double> &spikeT);
    int get_outC_Size();
};

#endif // IZHIRK_H
