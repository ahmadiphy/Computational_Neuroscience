#include "in.h"
#include "global_variables.h"
#include "izhirk.h"
#include "hire.h"
#include "matrixf.h"
#include "regular_net.h"
using namespace std;
double a;
double b;
double c;
double d;
double h;

int main()
{
    //
    //omp_set_num_threads(3);
    //----------------------------------------------------------------------------------------
    //-----------------------------  core dump error checking --------------------------------
    //----------------------------------------------------------------------------------------
    struct rlimit core_limit;
    core_limit.rlim_cur = RLIM_INFINITY;
    core_limit.rlim_max = RLIM_INFINITY;
    if (setrlimit(RLIMIT_CORE, &core_limit) < 0) {
        /* ERROR */
    }
    //----------------------------------------------------------------------------------------
    //---------------------------- geting system time and date -------------------------------
    //----------------------------------------------------------------------------------------
    time_t now = time(0)+16200;
    char* dt = ctime(&now);
    tm *gmtm = gmtime(&now);
    dt = asctime(gmtm);
    //----------------------------------------------------------------------------------------
    //----------------------------------  make directory -------------------------------------
    //----------------------------------------------------------------------------------------
    const int dir_err = mkdir("data", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    if (-1 == dir_err)
    {
        cout<<"Error Creating Directory"<<endl;
        cout<<"-------------- No Problem ---------------"<<endl;
    }
    //----------------------------------------------------------------------------------------
    //------------------------------------ Make Network --------------------------------------
    //----------------------------------------------------------------------------------------
    Regular_net rg1;
    Hire myh;
    int l=9;
    int m0=2,ll=pow(2,l);
    int N=m0*ll;
    Matrixf mf1;
    iMatrix mat1(N,iRow());
    //rg1.regT1(N,1,mat1);
    //rg1.regR1(N,0.5,mat1);
    myh.interModuleNM1(ll,m0,mat1);//sakht module haye shabake
    myh.intraModuleHMN1(l,m0,mat1,0.25,1);//sakht kole shabake
    //myh.intraModuleHMNfibolink(l,m0,mat1,0.75,1);//sakht kole shabake
    mf1.mat_shrink_to_fit(N,mat1);//behineh saziye ram
    //mf1.CnMatrix(N,mat1);
    //cin.get();
    //----------------------------------------------------------------------------------------
    //-------------------------------- info file write ---------------------------------------
    ostringstream info;
    info<<"./data/info.txt";
    ofstream infoPrint(info.str().c_str(),ios_base::binary);
    //infoPrint<<"."<< endl;
    infoPrint<<"------------------------------------------------------------"<<endl;
    infoPrint<<"--------- Simulation of Izhikevich Neurons on HMN ----------"<<endl;
    infoPrint<<"------------------------------------------------------------"<<endl;
    //----------------------------------------------------------------------------------------
    cout<<"The Runing Start date and time is: "<<dt <<endl;
    infoPrint<<"The Runing Start date and time is: "<<dt<<endl;
    //----------------------------------------------------------------------------------------
    //---------------------------- Initialize Global Variables -------------------------------
    //----------------------------------------------------------------------------------------
    a=0.02;
    b=0.2;
    c=-65;
    d=8;
    h=0.1;//Runge Kutta time step length or delta t
    double g=0.9;//synaptic weight
    infoPrint<<"Network Parameters m0 = "<< m0 <<", l = "<< l << endl;
    infoPrint<<"Model   Parameters a  = "<< a << ", b = "<< b <<", c = "<< c <<", d = "<<d<< endl;
    infoPrint<<"Runge-Kutta 4th    h  = "<< h << endl;
    //----------------------------------------------------------------------------------------
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    poisson_distribution<int> pdist(10);
    normal_distribution<double> distribution(10.0,3);
    uniform_real_distribution<double> udistribution(0.0,10.0);
    gamma_distribution<double> gdistribution(1,6);
    //----------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------
    int loopTry=2000;
    vector<IzhiRK> izhi(N);
    ostringstream outData,vuData,stimes,vD;
    vD<<"./data/vDegree.dat";
    ofstream vdPrint(vD.str().c_str(),ios_base::binary);
    double ri=0;
    for(int i=0;i<N;++i)
    {
        izhi[i].Connect(mat1[i]);
        izhi[i].myNum=i;
        ri=pdist(gen);
        izhi[i].Initialize(ri);
        vdPrint<<i<<" "<<izhi[i].outC.size()<<endl;
        //cout<<ri<<" "<<izhi[i].D<<endl;
    }
    //cin.get();
    //----------------------------------------------------------------------------------------
    //-------------------------------------  file writing  -----------------------------------
    //----------------------------------------------------------------------------------------
    //outData<<"./data/data_activity.dat";
    vuData<<"./data/data_vu.dat";
    stimes<<"./data/stimes_g"<<g<<".dat";
    //ofstream outPrint(outData.str().c_str(),ios_base::binary);
    ofstream vuPrint(vuData.str().c_str(),ios_base::binary);
    ofstream stPrint(stimes.str().c_str(),ios_base::binary);
    //----------------------------------------------------------------------------------------
    mf1.DeleteMatElement(N,mat1);
    cout<<"Dynamics... "<<endl;
    double step_time=0;
    cout<<step_time<<endl;
    for(int i=0;i<loopTry;++i)
    {
        step_time=step_time+h;
        cout<<fixed<<step_time<<endl;
        for(int j=0;j<N;++j)
        {
            //ri=pdist(gen);
            izhi[j].UpdateI();
            izhi[j].Run(step_time);
        }
        //
        for(int j=0;j<N;++j)
        {
            //if(izhi[j].Spike==1)
            //{
            for(int k=0;k<izhi[j].outC.size();++k)
            {
                int tar=izhi[j].outC[k];
               //cout<<izhi[tar].cashI<<endl;
                izhi[j].cashI+=g*(izhi[tar].v-izhi[j].v);
               //cout<<izhi[j].v<<" "<<izhi[tar].v<<" "<<g<<endl;
               // cout<<izhi[tar].cashI<<endl;
            }
            //}
            //outPrint<<izhi[j].Spike<<" ";
            //cout<<izhi[j].Spike<<" ";
        }
        //cout<<izhi[10].I<<endl;
        //outPrint<<endl;
        //cout<<endl;
        if(i>=190000)
        {
            for(int j=0;j<N;++j)
                vuPrint<<izhi[j].v<<" ";
            vuPrint<<endl;
        }
             //if(i>10000)
             //   cin.get();
    }
    for(int i=0;i<N;++i)
    {
        izhi[i].spikeTime.shrink_to_fit();
        for(int j=0;j<izhi[i].spikeTime.size();++j)
        {
            stPrint<<fixed<<izhi[i].spikeTime[j]<<" ";
        }
        stPrint<<endl;
    }
    //----------------------------------------------------------------------------------------
    //------------------------------------  Time printing  -----------------------------------
    //----------------------------------------------------------------------------------------
    now = time(0)+16200;
    dt = ctime(&now);
    gmtm = gmtime(&now);
    dt = asctime(gmtm);
    cout<<"The Runing End   date and time is: "<<dt<< endl;
    infoPrint<<"The Runing End   date and time is: "<<dt<< endl;
    infoPrint<<"------------------------------------------------------------"<<endl;
    //----------------------------------------------------------------------------------------

}
