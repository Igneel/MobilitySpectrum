#include <iostream>

using namespace std;


#include <math.h>
#include <vector>
#include "mobilityspectrum.h"


int main()
{
    vector<long double> B={0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2};
    vector<long double> sxx={44.359779, 44.311464, 44.182902, 43.987464, 43.74224, 43.466406, 43.177335, 42.888688, 42.609968, 42.346955, 42.102523};
    vector<long double> sxy={0, 0.6208139, 1.2227963, 1.7801449, 2.279831, 2.7172189, 3.0942103, 3.417015, 3.6939559, 3.9337932, 4.1446831};

    freopen("outPut.txt","w",stdout);

    cout << "Hello Log!" << endl;
    mobilitySpectrum c(B,sxx,sxy,11);
    int s=c.getResultSize();

    for(int i=0;i<s;++i)
    {
        cout << c.getResultEX(i)<<"\t"<< c.getResultEY(i)<<"\t"<< c.getResultHY(i) <<endl;
       // cout << c.getResultHX(i)<<"\t"<< c.getResultHY(i)<<endl;
    }

    return 0;
}

