#include "regular_net.h"

void Regular_net::regT1(int N,int M,iMatrix &a)
{
    for(int i=0;i<N;++i)
    {
        for(int j=1;j<=M;++j)
        {
            int k=i+j;
            if(i+j>=N)
            {
                k=k-N;
            }
            a[i].push_back(k);
            a[k].push_back(i);
        }
    }
}
void Regular_net::regR1(int N,double r,iMatrix & a)
{
    double n=N;
    double rLinks=r*n/2;
    int links=rLinks;
    std::random_device rd;
    std::mt19937 gen(rd());  // to seed mersenne twister.
    std::uniform_int_distribution<> distall(0,N-1);
    for(int i=0;i<links;++i)
    {
        int r1=distall(gen);
        int r2=distall(gen);
        int ch=0;
        for(int j=0;j<a[r1].size();++j)
        {
            if(a[r1][j]==r2)
            {
                ch=1;
                i--;
            }
        }
        if(ch==0)
        {
            a[r1].push_back(r2);
            a[r2].push_back(r1);
        }

    }

}
