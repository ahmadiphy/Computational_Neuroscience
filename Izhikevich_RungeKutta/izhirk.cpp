#include "izhirk.h"

using namespace std;

IzhiRK::IzhiRK()
{
    spikeNumber=0;
    state=0;
    spikeCon=0;
    cashI=0;
    D=0;
    rI=0;
    I=0;
}
//----------------------------------------------------------------------------------------
//--------------------------   Izhikevich Model v' Function   ----------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::fv(double v,double u,double I)
{
    double result=0.04*pow(v,2)+5*v+140-u+I;
    return result;
}
//----------------------------------------------------------------------------------------
//--------------------------   Izhikevich Model u' Function   ----------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::fu(double v,double u)
{
    double result=a*(b*v-u);
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   step 1 in runge kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::k1v(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fv(v,u,I);
    return result;
}
//
double IzhiRK::k1u(double v,double u,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fu(v,u);
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   Step 2 in Runge Kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::k2v(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fv(v+0.5*k1v(v,u,I,h),u+0.5*k1u(v,u,h),I);
    return result;
}
//
double IzhiRK::k2u(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fu(v+0.5*k1v(v,u,I,h),u+0.5*k1u(v,u,h));
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   Step 3 in Runge Kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::k3v(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fv(v+0.5*k2v(v,u,I,h),u+0.5*k2u(v,u,I,h),I);
    return result;
}
//
double IzhiRK::k3u(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fu(v+0.5*k2v(v,u,I,h),u+0.5*k2u(v,u,I,h));
    return result;
}
//----------------------------------------------------------------------------------------
//---------------------------   Step 4 in Runge Kutta   ----------------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::k4v(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fv(v+k3v(v,u,I,h),u+k3u(v,u,I,h),I);
    return result;
}
//
double IzhiRK::k4u(double v,double u,double I,double h)//fr is result of f(x,...) function and h is time step length
{
    double result=h*fu(v+k3v(v,u,I,h),u+k3u(v,u,I,h));
    return result;
}
//----------------------------------------------------------------------------------------
//-----------------------------   Runge Kutta Result   -----------------------------------
//----------------------------------------------------------------------------------------
double IzhiRK::RKv(double v,double u,double I)
{
    double result=v+0.166666667*(k1v(v,u,I,h)+2*k2v(v,u,I,h)+2*k3v(v,u,I,h)+k4v(v,u,I,h));
    return result;
}
//
double IzhiRK::RKu(double v, double u, double I)
{
    double result=u+0.166666667*(k1u(v,u,h)+2*k2u(v,u,I,h)+2*k3u(v,u,I,h)+k4u(v,u,I,h));
    return result;
}
//----------------------------------------------------------------------------------------
//-----------------------------   Initialize Neurons   -----------------------------------
//----------------------------------------------------------------------------------------
void IzhiRK::Initialize(double ri)
{
    rI=ri;
    cashI=0;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<double> vd(-60,30);
    uniform_real_distribution<double> ud(-60,30);
    v=vd(gen);
    u=b*ud(gen);
    D=outC.size();
    state=1;
    I=rI;
}
//----------------------------------------------------------------------------------------
//-------------------------------   Connect Neurons   ------------------------------------
//----------------------------------------------------------------------------------------
void IzhiRK::Connect(std::vector<int> &cv)
{
    for(int i=0;i<cv.size();++i)
    {
        outC.push_back(cv[i]);
    }
    outC.shrink_to_fit();
}

//----------------------------------------------------------------------------------------
//-----------------------------   Run Neural Dynamics   ----------------------------------
//----------------------------------------------------------------------------------------
void IzhiRK::Run(double step_time)
{
    if(state==1)
    {
        double tempV,tempU;
        tempV=RKv(v,u,I);//solve v equation RK4
        tempU=RKu(v,u,I);//solve u equation RK4
        v=tempV;
        u=tempU;
        if(v>=30)
        {
            v=c;
            u=u+d;
            //Spike=1;
            //spikeNumber++;
            //spikeCon++;
            if(step_time>=10000)
                spikeTime.push_back(step_time);
            //cout<<Spike<<" ";
        }//else
        //{
        //    Spike=0;
        //}
        //result[0]=v;
        //result[1]=u;
    }else
    {
        cout<<"   ERROR! (initial value problem)  "<<endl;
        exit (EXIT_FAILURE);
    }
}
//----------------------------------------------------------------------------------------
//------------------------------   Get Random Current   ----------------------------------
//----------------------------------------------------------------------------------------
void IzhiRK::UpdateI()
{
    I=rI+(cashI/D);
    cashI=0;
    //cout<<rI<<" "<<I<<endl;

}
//----------------------------------------------------------------------------------------
//-------------------------------   Phase of Neuron   ------------------------------------
//----------------------------------------------------------------------------------------
