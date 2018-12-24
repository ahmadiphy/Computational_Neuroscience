#include "hire.h"
#include "ern.h"
#include "matrixf.h"
//
using namespace::std;
//
//
void Hire::interModuleCC(int ll, int nn, iMatrix& aa)
{
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*nn;
        modulE=(i+1)*nn;
        for(int j=modulB;j<modulE;j++)
        {
            for(int k=modulB;k<modulE;k++)
            {
                aa[j][k]=1;
            }
            aa[j][j]=0;
        }
    }
}
//
void Hire::intraModuleCHMN1(int l, int nn, iMatrix &aa,double prob,double alpha)
{
    //
    int ll=pow(2,l);
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int hi;
    hi=l;
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        double hLevel=ii+1;
        double lProbe=0;
        lProbe=alpha*(pow(prob,hLevel));
        cout<<lProbe<<endl;
        //
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            //cout<<endl;
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            //cout<<dis1B<<" - "<<dis1E<<endl;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            //cout<<dis2B<<" - "<<dis2E<<endl;
            int tch=0;
            while(tch==0){
                for(int ch1=dis1B;ch1<=dis1E;++ch1)
                {
                    for(int ch2=dis2B;ch2<=dis2E;++ch2)
                    {
                        double myprob=0;
                        myprob=dist(gen);
                        if(myprob<=lProbe)
                        {
                            aa[ch1][ch2]=1;
                            aa[ch2][ch1]=1;
                            tch=1;
                        }
                    }
                }
            }
        }

    }
}
//
void Hire::intraModuleRC(int ll, int nn, iMatrix& aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0,mch=0;//mc nubmer of module out connectivity
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                if((r1*r2)!=mch)
                {
                    aa[r1][r2]=1;
                    aa[r2][r1]=1;
                    mch=r1*r2;
                    mc++;
                }
            }while (mc<2);
        }

    }
}
//
void Hire::EinterModuleCC(int ll, int nn, iMatrix& aa)
{
    int nl=0;
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*nn;
        modulE=(i+1)*nn;
        for(int ii=0;ii<nn;++ii)
        {
            for(int j=modulB;j<nl;++j)
            {
                aa[nl].push_back(j);
            }
            for(int j=nl+1;j<modulE;++j)
            {
                aa[nl].push_back(j);
            }
            nl++;
        }

    }
}
//
void Hire::EintraModuleRC(int ll, int nn, iMatrix& aa)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0,mch=0;//mc nubmer of module out connectivity
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                if((r1*r2)!=mch)
                {
                    aa[r1].push_back(r2);
                    aa[r2].push_back(r1);
                    mch=r1*r2;
                    mc++;
                }
            }while (mc<2);
        }

    }
}
//
void Hire::interModuleCCT2(int ll, int nn, iMatrix &aa,double avgcon)
{
    if(avgcon<nn)
    {
        ERn eer1;
        Matrixf mm1;
        iMatrix a0(nn,iRow(nn));
        int modulB,modulE;
        for(int i=0;i<ll;++i)
        {
            modulB=i*nn;
            modulE=(i+1)*nn;
            int jj=0;
            mm1.zMatrix(nn,a0);
            eer1.ERnetwork(avgcon,nn,a0);
            for(int j=modulB;j<modulE;j++)
            {
                int kk=0;
                for(int k=modulB;k<modulE;k++)
                {
                    aa[j][k]=a0[jj][kk];
                    kk++;
                }
                aa[j][j]=0;
                jj++;
            }
        }
    }else
    {
        cout<<"matrix not changed "<<endl<<"the avrage connectivity is biger than size of network"<<endl;
        //
    }
}
//
void Hire::intraModuleRCT2(int ll, int nn, iMatrix &aa,int outcon)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0;//mc nubmer of module out connectivity
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                if(aa[r1][r2]==0)
                {
                    aa[r1][r2]=1;
                    aa[r2][r1]=1;
                    mc++;
                }
            }while (mc<outcon);
        }

    }
}
//
void Hire::intraModuleN0(int ll, int nn, iMatrix& aa, int bn)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0;//mc nubmer of module out connectivity
            do{
                r1=0;
                r2=0;
                r1=dis1(gen);
                r2=dis2(gen);
                if(aa[r1][r2]==0)
                {
                    aa[r1][r2]=1;
                    aa[r2][r1]=1;
                    mc++;
                }
            }while (mc<bn);
        }

    }
}
//
void Hire::intraModuleN1(int ll, int nn, iMatrix& aa, int bn)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0;//mc nubmer of module out connectivity
            do{
                r1=0;
                r2=0;
                r1=dis1(gen);
                r2=dis2(gen);
                if(aa[r1][r2]==0)
                {
                    aa[r1][r2]=1;
                    aa[r2][r1]=1;
                    mc++;
                }
            }while (mc<bn);
        }

    }
}
//
void Hire::intraModuleFfibo(int ll, int nn, iMatrix &aa)
{
    //netType="fibo_Hirecial_network";
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    //----
    vector<int> fibo(hi+1);
    int outcon=0;
    fibo[0]=1;
    fibo[1]=1;
    for(int f=2;f<hi;++f)
    {
        fibo[f]=fibo[f-1]+fibo[f-2];
    }
    //----
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        outcon=fibo[(hi-ii)-1];
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0;//mc nubmer of module out connectivity
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                if(aa[r1][r2]==0)
                {
                    aa[r1][r2]=1;
                    aa[r2][r1]=1;
                    mc++;
                }
            }while (mc<outcon);
        }

    }
}
//
void Hire::intraModuleFfiboV(int ll, int nn, iMatrix &aa)
{
    //netType="fibo_Hirecial_network";
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    //----
    vector<int> fibo(hi+1);
    int outcon=0;
    fibo[0]=1;
    fibo[1]=1;
    for(int f=2;f<hi;++f)
    {
        fibo[f]=fibo[f-1]+fibo[f-2];
    }
    //----
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        outcon=fibo[ii];
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0;//mc nubmer of module out connectivity
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                if(aa[r1][r2]==0)
                {
                    aa[r1][r2]=1;
                    aa[r2][r1]=1;
                    mc++;
                }
            }while (mc<outcon);
        }

    }
}
//
void Hire::interModuleNM1(int ll, int nn, iMatrix& aa)
{
    int nl=0;
    int modulB,modulE;
    for(int i=0;i<ll;++i)
    {
        modulB=i*nn;
        modulE=(i+1)*nn;
        for(int ii=0;ii<nn;++ii)
        {
            for(int j=modulB;j<nl;++j)
            {
                aa[nl].push_back(j);
            }
            for(int j=nl+1;j<modulE;++j)
            {
                aa[nl].push_back(j);
            }
            nl++;
        }

    }
}
//
int Hire::NM1Ccheck(int m0,int c1,int c2,iMatrix& aa)
{
    int ch=0;
    for(int i=m0-1;i<aa[c1].size();++i)
    {
        if(aa[c1][i]==c2)
            ch=1;
    }
    return ch;
}
//
void Hire::intraModuleNM1fibo(int ll, int nn, iMatrix &aa)// Error!!! have a problem
{
    cout<<"ERROR!"<<endl;
    //netType="fibo_Hirecial_network";
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    //----
    vector<int> fibo(hi+1);
    int outcon=0;
    fibo[0]=1;
    fibo[1]=1;
    //cout<<"NO ERROR 2"<<endl;
    for(int f=2;f<hi;++f)
    {
        fibo[f]=fibo[f-1]+fibo[f-2];
    }
    //----
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        outcon=fibo[(hi-ii)-1];
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<int> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<int> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0;//mc nubmer of module out connectivity
            //cout<<"NO ERROR 3"<<endl;
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                //check r1,r2
                //cout<<"NO ERROR 4"<<endl;
                int ccheck=NM1Ccheck(nn,r1,r2,aa);//check that r1 is not connected to r2 and V
                //cout<<"NO ERROR 5"<<endl;
                cout<<dis1B<<" "<<dis1E<<"__"<<dis2B<<" "<<dis2E<<endl;
                cout<<r1<<"_"<<r2<<"--"<<aa.size()<<endl;
                if(ccheck==0)
                {
                    cout<<"NO ERROR 6"<<endl;

                    aa[r1].push_back(r2);
                    aa[r2].push_back(r1);
                    mc++;
                }
                //cout<<"NO ERROR 7"<<endl;

            }while (mc<outcon);
            cout<<"NO ERROR 8"<<endl;

        }

    }
}
//
void Hire::intraModuleNM1RC(int ll, int nn, iMatrix& aa,int avgnl)
{
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            int r1=0,r2=0;
            int mc=0,mch=0;//mc nubmer of module out connectivity
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                int ccheck=NM1Ccheck(nn,r1,r2,aa);//check that r1 is not connected to r2 and V
                if(ccheck==0)
                {
                    aa[r1].push_back(r2);
                    aa[r2].push_back(r1);
                    mc++;
                }
            }while (mc<avgnl);
        }

    }
}
//
void Hire::intraModuleNM1_MHMN(int ll, int nn, iMatrix& aa,double s,double alpha)
{
    double phi=2;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        double hLevel=ii+1;
        double sl=-1*s*hLevel;
        double lProbe=0;
        lProbe=alpha*(pow(phi,sl));
        //
        int il=ip;
        ip=2*ip;
        //
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            for(int ch1=dis1B;ch1<=dis1E;++ch1)
            {
                for(int ch2=dis2B;ch2<=dis2E;++ch2)
                {
                    double myprob=0;
                    myprob=dist(gen);
                    if(myprob<=lProbe)
                    {
                        aa[ch1].push_back(ch2);
                        aa[ch2].push_back(ch1);
                    }
                }
            }
        }

    }
}
//
void Hire::intraModuleNM1_SR_MHMN(int ll, int nn, iMatrix &aa, double s, double alpha)
{
    //cout<<"fin1"<<endl;
    double phi=2;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    double hir,hirs;
    hirs=ll;
    int mo=4;
    hir=log2(hirs)/2;
    int hi;
    hi=int(hir);
    int ip=1;
    //cout<<"fin2"<<endl;
    for(int ii=0;ii<hi;ii++)
    {
        //cout<<"fin3"<<endl;
        double hLevel=ii+1;
        double sl=-1*s*hLevel;
        double lProbe=0;
        lProbe=alpha*(pow(phi,sl));
        //
        int il=ip;
        ip=mo*ip;
        //
        for(int i=0;i<ll-il;i+=ip)
        {
            //cout<<i;
            //cout<<"----"<<endl;
            for(int j=0;j<mo-1;++j)
            {
                int dis1B=0,dis1E=0;
                dis1B=(i+j*il)*nn;
                dis1E=((i+j*il+il)*nn)-1;
                //cout<<"1beg:"<<dis1B<<endl;
                //cout<<"1end:"<<dis1E<<endl;
                //
                //cout<<"*****"<<endl;
                for(int jj=j;jj<mo-1;++jj)
                {
                    int dis2B=0,dis2E=0;
                    dis2B=(jj*il+i+il)*nn;
                    dis2E=((jj*il+i+(2*il))*nn)-1;
                    //cout<<"2beg:"<<dis2B<<endl;
                    //cout<<"2end:"<<dis2E<<endl;
                    for(int ch1=dis1B;ch1<=dis1E;++ch1)
                    {
                        for(int ch2=dis2B;ch2<=dis2E;++ch2)
                        {
                            double myprob=0;
                            myprob=dist(gen);
                            if(myprob<=lProbe)
                            {
                                aa[ch1].push_back(ch2);
                                aa[ch2].push_back(ch1);
                            }
                        }
                    }
                }
            }
            //
        }
    }

}
//
void Hire::intraModuleNM2_MHMN(int ll, int nn, iMatrix& aa,double s,double alpha)
{
    double phi=2;
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    double hir,hirs;
    hirs=ll;
    hir=log2(hirs);
    int hi;
    hi=int(hir);
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {

        double hLevel=ii+1;
        double sl=-1*s*hLevel;
        double lProbe=0;
        lProbe=alpha*(pow(phi,sl));
        //
        int il=ip;
        ip=2*ip;
        //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        if(ii>=0)
        {
            int lN=0;
            lN=lProbe*pow(nn,(2*hLevel));
            for(int i=0;i<ll-il;i+=ip)
            {
                int dis1B=0,dis1E=0;
                dis1B=i*nn;
                dis1E=((i+il)*nn)-1;
                uniform_int_distribution<> dis1(dis1B, dis1E);
                int dis2B=0,dis2E=0;
                dis2B=(i+il)*nn;
                dis2E=((i+(2*il))*nn)-1;
                uniform_int_distribution<> dis2(dis2B, dis2E);
                int r1=0,r2=0;
                int mc=0;//mc nubmer of module out connectivity
                do{
                    r1=dis1(gen);
                    r2=dis2(gen);
                    int ccheck=NM1Ccheck(nn,r1,r2,aa);//check that r1 is not connected to r2 and V
                    if(ccheck==0)
                    {
                        aa[r1].push_back(r2);
                        aa[r2].push_back(r1);
                        mc++;
                    }
                }while (mc<lN);
            }
            //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
        }else
        {
            for(int i=0;i<ll-il;i+=ip)
            {
                int dis1B=0,dis1E=0;
                dis1B=i*nn;
                dis1E=((i+il)*nn)-1;
                int dis2B=0,dis2E=0;
                dis2B=(i+il)*nn;
                dis2E=((i+(2*il))*nn)-1;
                for(int ch1=dis1B;ch1<=dis1E;++ch1)
                {
                    for(int ch2=dis2B;ch2<=dis2E;++ch2)
                    {
                        double myprob=0;
                        myprob=dist(gen);
                        if(myprob<=lProbe)
                        {
                            aa[ch1].push_back(ch2);
                            aa[ch2].push_back(ch1);
                        }
                    }
                }
            }
        }

    }
}
//
void Hire::intraModuleHMN1(int l, int nn, iMatrix &aa,double prob,double alpha)
{
    //
    int ll=pow(2,l);
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int hi;
    hi=l;
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        double hLevel=ii+1;
        double lProbe=0;
        lProbe=alpha*(pow(prob,hLevel));
        //cout<<lProbe<<endl;
        //
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            //cout<<endl;
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            //cout<<dis1B<<" - "<<dis1E<<endl;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            //cout<<dis2B<<" - "<<dis2E<<endl;
            int tch=0;
            //cout<<endl;
            while(tch==0){
                for(int ch1=dis1B;ch1<=dis1E;++ch1)
                {
                    for(int ch2=dis2B;ch2<=dis2E;++ch2)
                    {
                        double myprob=0;
                        myprob=dist(gen);
                        if(myprob<=lProbe)
                        {
                            //cout<<ch1<<" - "<<ch2<<endl;
                            aa[ch1].push_back(ch2);
                            aa[ch2].push_back(ch1);
                            tch=1;
                        }
                    }
                }
            }
        }

    }
}
//
void Hire::Hanoi01(int m0,int l, int nn, iMatrix& aa)
{
    int ll=pow(2,l);
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int hi;
    hi=l;// tedade marhaleha der maghaleye Li haman "g" ast
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        double hLevel=ii+1;
        double lProbe=0;
        lProbe=1/(pow(double(m0),2)*(pow(4,double(hLevel))-pow(2,double(hLevel+1))+1));
        //cout<<lProbe<<endl;
        //
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            //cout<<endl;
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            //cout<<dis1B<<" - "<<dis1E<<endl;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            //cout<<dis2B<<" - "<<dis2E<<endl;
            int tch=0;
            //cout<<endl;
            while(tch==0){
                for(int ch1=dis1B;ch1<=dis1E;++ch1)
                {
                    for(int ch2=dis2B;ch2<=dis2E;++ch2)
                    {
                        double myprob=0;
                        myprob=dist(gen);
                        if(myprob<=lProbe)
                        {
                            //cout<<ch1<<" - "<<ch2<<endl;
                            aa[ch1].push_back(ch2);
                            aa[ch2].push_back(ch1);
                            tch=1;
                        }
                    }
                }
            }
        }

    }
}
//
void Hire::intraModuleHMNfibolink(int l, int nn, iMatrix &aa,double prob,double alpha)
{
    //
    int ll=pow(2,l);
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int hi;
    hi=l;
    int ip=1;
    vector<int> fibo(hi+1);
    int outcon=0;
    fibo[0]=1;
    fibo[1]=1;
    //cout<<"NO ERROR 2"<<endl;
    for(int f=2;f<hi;++f)
    {
        fibo[f]=fibo[f-1]+fibo[f-2];
    }
    for(int ii=0;ii<hi;ii++)
    {
        double hLevel=ii+1;
        double lProbe=0;
        //lProbe=alpha*(pow(prob,hLevel));
        outcon=fibo[(hi-ii)-1];
        //cout<<lProbe<<endl;
        //
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int r1=0,r2=0;
            //cout<<endl;
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            uniform_int_distribution<> dis1(dis1B, dis1E);
            //cout<<dis1B<<" - "<<dis1E<<endl;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            uniform_int_distribution<> dis2(dis2B, dis2E);
            //cout<<dis2B<<" - "<<dis2E<<endl;
            int tch=0;
            //cout<<endl;
            //
            int mc=0;
            do{
                r1=dis1(gen);
                r2=dis2(gen);
                //check r1,r2
                //cout<<"NO ERROR 4"<<endl;
                int ccheck=NM1Ccheck(nn,r1,r2,aa);//check that r1 is not connected to r2 and V
                //cout<<"NO ERROR 5"<<endl;
                //cout<<dis1B<<" "<<dis1E<<"__"<<dis2B<<" "<<dis2E<<endl;
                //cout<<r1<<"_"<<r2<<"--"<<aa.size()<<endl;
                if(ccheck==0)
                {
                    //cout<<"NO ERROR 6"<<endl;
                    aa[r1].push_back(r2);
                    aa[r2].push_back(r1);
                    mc++;
                }
                //cout<<"NO ERROR 7"<<endl;

            }while (mc<outcon);
            //
        }

    }
}
//
void Hire::intraModuleHMN2(int l, int nn, iMatrix &aa,double prob,double alpha)
{
    //
    int ll=pow(2,l);
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int hi;
    hi=l;
    int ip=1;
    double lProbe=prob;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            //cout<<endl;
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            //cout<<dis1B<<" - "<<dis1E<<endl;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            //cout<<dis2B<<" - "<<dis2E<<endl;
            int tch=0;
            //cout<<endl;
            while(tch==0){
                for(int ch1=dis1B;ch1<=dis1E;++ch1)
                {
                    for(int ch2=dis2B;ch2<=dis2E;++ch2)
                    {
                        double myprob=0;
                        myprob=dist(gen);
                        if(myprob<=lProbe)
                        {
                            //cout<<ch1<<" - "<<ch2<<endl;
                            aa[ch1].push_back(ch2);
                            aa[ch2].push_back(ch1);
                            tch=1;
                        }
                    }
                }
            }
        }

    }
}

//
void Hire::intraModuleHMN3(int l, int nn, iMatrix &aa,double prob,double alpha)
{
    //
    int ll=pow(2,l);
    random_device rd;
    mt19937 gen(rd());  // to seed mersenne twister.
    uniform_real_distribution<> dist(0,1);
    int hi;
    hi=l;
    int ip=1;
    for(int ii=0;ii<hi;ii++)
    {
        int il=ip;
        ip=2*ip;
        for(int i=0;i<ll-il;i+=ip)
        {
            int dis1B=0,dis1E=0;
            dis1B=i*nn;
            dis1E=((i+il)*nn)-1;
            int dis2B=0,dis2E=0;
            dis2B=(i+il)*nn;
            dis2E=((i+(2*il))*nn)-1;
            aa[dis1B].push_back(dis2E);
            aa[dis2E].push_back(dis1B);
            if(ii==0)
            {
                aa[dis1B].push_back(dis2B);
                aa[dis2B].push_back(dis1B);
            }
        }
    }
}
//
