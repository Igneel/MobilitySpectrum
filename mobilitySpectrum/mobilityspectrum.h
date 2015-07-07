#ifndef MOBILITYSPECTRUM_H
#define MOBILITYSPECTRUM_H

#include <math.h>
#include <vector>
#include <iostream>

using namespace std;
typedef vector< pair<long double, long double> > TLineSeries;

class mobilitySpectrum
{
private:
    const int PointPerInt=50;

    typedef vector <long double> Data_spektr ;

    typedef vector <long double> ImageDat;
    typedef vector <vector <long double> > mat ;
    typedef vector <long double> Dat1 ;
    typedef vector <vector <long double> > Dat2 ;
    typedef vector <vector <long double> > Dat3 ;



    int NumberOfPoints,Power_spektr,GridPoints;

    Data_spektr MagField_spektr,GxxExp,GxyExp;

    void logger();

    ImageDat IntGxx;
    ImageDat IntGxy;
    ImageDat IntMagField;
    ImageDat Spectr_e;
    ImageDat Spectr_p;
    ImageDat Mobility;

    long int SizeData;
    Dat1 B_spektr,Gxx_sp,Gxx_MC, Gxy_MC, Gxy_sp,Xr,Lv,Xv,Mv,Vpr;
    Dat2 Am,Qm,Cl,Cr,Cl_t,Cr_t,Cm,Cm_t          ;
    bool bulua;
    int  MSRight,MSLeft;

    long double W,F_s,A1,An,B1,Bn, Ves1, Ves2,Mu_max, Min_Spectr,Coef1,Coef2,Mu_min;

    TLineSeries electronMobilitySpectrum, holeMobilitySpectrum;

    int  MaxPoints;

public:

    size_t getResultSize();

    void  InitArray();
    void  InitArray2();

    void  GetCoef(const Data_spektr &A,const Data_spektr &X, const long double b,
        long double &p0, long double &p1, long double &p2);
    void  GetLnLimits(int &Xmin, int &Xmax );

    void SetLength (vector<vector<long double> > &v,const size_t size1,const size_t size2);

    void  gram(const int N,const int m,const int l, Data_spektr & x, Data_spektr & f, mat & a);
    void  gauss(const int N, mat & a, Data_spektr & x);

    void  fi(const int n,const int m,const int l, Data_spektr & c, Data_spektr & x,
        const long double x1, long double &s);

    void  BAS(const int n,const int m,const int L,const long double x1, Data_spektr & x, Data_spektr & t);

    void  AddExpPoints(TLineSeries &ExpXX, TLineSeries &ExpXY);

    void  CS(Data_spektr& X,Data_spektr& F,Data_spektr& C,const long double p1,const long double pn);
    long double Sp(Data_spektr & X,Data_spektr &F, Data_spektr & C,const long double x1);
    void  MakeInterpolate(TLineSeries &Gxx, TLineSeries &Gxy,
        TLineSeries &ExpXX, TLineSeries &ExpXY);

    void  MakeMNK(const bool a, TLineSeries &Gxx, TLineSeries &Gxy, TLineSeries &ExpXX, TLineSeries &ExpXY);
    void  MakeLagranj();
    void  Tred2(const int n, Dat1 & d, Dat1 & e,
                     Dat2 & a, Dat2 & z, bool & fail);
    void  Imtql2(const int n,const long double macheps, Dat1 & d, Dat1 & e,
                     Dat2 & z, bool & fail);
    long double GetElem(const int j1, const int k1, const int i1);
    void  MakeMatrC();
    void  MakeMatrA();
    void  InverseMatrC(Dat2 & Ci,Dat2 & C,long double & Su,const int NP);
    long double S_s(const long double Mi);
    void  MobilitySpectrumFunc(TLineSeries &LineSeries1, TLineSeries &Series5);

    long double getResultEY(const int i);

    long double getResultEX(const int i);

    long double getResultHY(const int i);

    long double getResultHX(const int i);

    mobilitySpectrum(Data_spektr &MagneticFieldP, Data_spektr &Exx,
                     Data_spektr &Exy,const int size);
};

#endif // MOBILITYSPECTRUM_H
