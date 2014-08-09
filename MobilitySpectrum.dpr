library MobilitySpectrum;
uses
  Windows,  SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Grids, Clipbrd, Buttons, StdCtrls, StrUtils, ExtCtrls, Menus,
  CheckLst, Math, ComCtrls, Gauges,
  Series,Chart,TeEngine, TeeProcs, Messages,Dialogs;

  // optim1

  {
  LineSeries1
  Series5

  Это то, куда выводятся спектры подвижности (электроны и дырки соотв.)

  Gxx
  Gxy

  Это графики расчитанных тензоров проводимости

  ExpXX
  ExpXY

  Экспериментальные значения компонент тензора проводимости


  Chart5 и Chart6 - в них лежат графики компонент тензора проводимости.

  }



{$R *.res}
const
  //MaxPoints=100;
  PointPerInt=50;
  //MaxParameters=8;    // Максимальное число параметров


  {  Mobility                         [ m**2/V*s ]     ;
      Magnetic Field                   [ Tesla ]        ;
      Electric Field                   [ V/m ]          ;
      Surface Recombination Velocity   [ m/s ]          ;
      All other values in SI                             ;   }
  type
    Data_spektr=array of extended;
    

    ImageDat=array of extended;  // (PointPerInt-1)*PointPerInt
  PImageDat=^ImageDat;
  mat= array of array of Extended;
  Dat1=array of Extended;
  Dat2= array of array of Extended;
  Dat3= array of array of Extended;

  {
  mat=array[0..MaxPoints,0..MaxPoints] of extended;
  Dat1=array[1..MaxPoints] of extended;
  Dat2=array[1..MaxPoints,1..MaxPoints] of extended;
  Dat3=array[1..MaxPoints,1..2*MaxPoints] of extended;

   }

  var
    /////////////// спектр подвижности///////////
  NumberofPoints,Power_spektr,GridPoints :word;
  MagField_spektr,GxxExp,GxyExp  :Data_spektr;


  outSpectrEX,outSpectrEY,outSpectrHX,outSpectrHY:Data_spektr;

  IntGxx,IntGxy,IntMagField,Spectr_e,Spectr_p,Mobility:ImageDat;
  
  QSpectr_e,QSpectr_p,Axx,Axy,Axx_d,Axy_d :PImageDat;

  SizeData                               :longint;
  B_spektr,Gxx_sp,Gxx_MC,Gxy_MC,
  Gxy_sp,Xr,Lv,Xv,Mv,Vpr                 :Dat1;
  Am,Qm,Cl,Cr,Cl_t,Cr_t,Cm,Cm_t          :Dat2;
  bulua                                  :boolean;
  MSRight,MSLeft                         :integer;
  Min_Spectr,Coef1,Coef2,Mu_min,Mu_max,
  W,F_s,A1,An,B1,Bn, Ves1, Ves2          :extended;

  electronMobilitySpectrum,holeMobilitySpectrum:TLineSeries;

  MaxPoints : Integer;
{ Important note about DLL memory management: ShareMem must be the
  first unit in your library's USES clause AND your project's (select
  Project-View Source) USES clause if your DLL exports any procedures or
  functions that pass strings as parameters or function results. This
  applies to all strings passed to and from your DLL--even those that
  are nested in records and classes. ShareMem is the interface unit to
  the BORLNDMM.DLL shared memory manager, which must be deployed along
  with your DLL. To avoid using BORLNDMM.DLL, pass string information
  using PChar or ShortString parameters. }



  ////////////////////////////////////////////////////////////////////////////
/////////////////////// НАЧАЛО "ХОЛЛ. ПОДВИЖНОСТЬ" /////////////////////////
////////////////////////////////////////////////////////////////////////////


 procedure InitArray;
begin

   SetLength(IntMagfield,SizeData);

   SetLength(IntGxx,SizeData);

   SetLength(IntGxy,SizeData);
end;

procedure InitArray2;
begin

   SetLength(Spectr_e,SizeData);

   SetLength(Spectr_P,SizeData);
   
   SetLength(Mobility,SizeData);
end;

 

// В общем это освобождение памяти
procedure ProcDLL(Reason: Integer);
begin
  if Reason = DLL_PROCESS_DETACH then
  begin
    IntMagField:=nil;
    IntGxx:=nil;
    IntGxy:=nil;
    Spectr_e:=nil;
    Spectr_p:=nil;
    Mobility:=nil;
    Axx:=nil;
    Axy:=nil;
    Axx_d:=nil;
    Axy_d:=nil;
  end;
end;




procedure GetCoef(var A,X:Data_spektr; b:extended;var P0,P1,P2:extended);
var i,j:word;
     r,r1,r2,q,s:extended;
 begin
    p0:=0;p1:=0;p2:=0;
    for i:=0 to NumberOfPoints do
     begin
      r:=1;r1:=0;r2:=0;
      for j:=0 to NumberOfPoints do
       if i<>j then begin
        q:=X[i]-X[j];
        s:=b-X[j];r2:=(r2*s+2*r1)/q;
        r1:=(r1*s+r)/q; r:=r*s/q;
       end;
      p0:=p0+r*A[i];p1:=p1+r1*A[i];p2:=p2+r2*A[i];
     end;
 end;

// Эта функция также нужна разве что мне для справки, т.к. вывод на график будет
// в основной программе.
// Или нет, не совсем понятно - возможно программа использует построенные графики
// для дальнейших расчетов, что тоже не совсем правильно.
// строит экспериментальные точки на графиках компонент тензора
procedure AddExpPoints(ExpXX:TPointSeries;ExpXY:TPointSeries);
var i:word;
begin
 with ExpXX as TPointSeries do
  begin
   Clear;
   Pointer.HorizSize:=2;
   Pointer.VertSize:=2;
   for i:=0 to NumberOfPoints do
    if MagField_spektr[i]=0 then AddXY(0.001,GxxExp[i],'',clTeeColor)
     else AddXY(MagField_spektr[i],GxxExp[i],'',clTeeColor);
   end;
  with ExpXY as TPointSeries do
   begin
    Clear;
    Pointer.HorizSize:=2;
    Pointer.VertSize:=2;
    for i:=0 to NumberOfPoints do
     if MagField_spektr[i]=0 then AddXY(0.001,GxyExp[i],'',clTeeColor)
       else AddXY(MagField_spektr[i],GxyExp[i],'',clTeeColor);
   end;
end;

procedure GetLnLimits(var Xmin,Xmax:integer);
var f:extended;
begin
  //ShowMessage('Inside Limits');
   if MagField_spektr[0]>0 then
    begin
      //ShowMessage('Go to if first');
     f:=ln(MagField_spektr[0])/ln(10);
     Xmin:=trunc(f);
    end
     else XMin:=-3;

    // ShowMessage('End up with choises');
    // ShowMessage('NumberOfPoints is'+IntToStr(NumberofPoints));
    f:=ln(MagField_spektr[NumberOfPoints])/ln(10);
    //ShowMessage('Step very carefully');
    Xmax:=trunc(f);
    //ShowMessage('Run if');
    if frac(f)>0.001 then Xmax:=Xmax+1;
    //ShowMessage('End GetLnLimits');
end;




procedure MakeInterpolate(Gxx:TLineSeries;Gxy:TLineSeries;ExpXX,ExpXY:TPointSeries);
var temp_l,temp_t,AGxx,AGxy,AField:Data_spektr;
    sf,lm,p,p1:extended;
    i,j:word;
    Lmin,Lmax,k:integer;
procedure CS(var X,F,C:Data_spektr;p1,pn:extended);
var i,j,m:word;
    K:Data_spektr;
    A,B,R:extended;
begin
   SetLength(K,MaxPoints);

   K[1]:=0;
   C[1]:=P1;
   A:=X[1]-X[0];
   B:=X[2]-X[1];
   k[2]:=(3*((F[2]-F[1])/(X[2]-X[1])-(F[1]-F[0])/(X[1]-X[0]))-
       (X[1]-X[0])*C[1])/2/(A+B);
   C[2]:=B/2/(A+B);
   for i:=3 to NumberOfPoints do
    begin
     j:=i-1;m:=j-1;
     A:=X[i]-X[j];B:=X[j]-X[m];R:=2*(A+B)-B*C[j];C[i]:=A/R;
     K[i]:=(3*((F[i]-F[j])/A-(F[j]-F[m])/B)-B*K[j])/R;
    end;
   C[NumberOfPoints]:=K[NumberofPoints]-C[NumberofPoints]*pn;
   for i:=NumberOfPoints-1 downto 2 do C[i]:=K[i]-C[i]*C[i+1];
end;
function Sp(var X,F,C:Data_spektr;x1:extended):extended;
var i,j:word;
    A,B,D,Q,R,P:extended;
begin

i:=1; while (X1>X[i]) do inc(i);
 j:=i-1;A:=F[j];B:=X[j];Q:=X[i]-B;
 R:=X1-B;P:=C[i];D:=C[i+1];
 B:=(F[i]-A)/q-(D+2*P)*Q/3;D:=(D-P)/3/Q;
 SP:=A+R*(B+R*(P+D*R));
end;
begin
    SetLength(temp_l,MaxPoints);
    SetLength(temp_t,MaxPoints);
    SetLength(AGxx,MaxPoints);
    SetLength(AGxy,MaxPoints);
    SetLength(AField,MaxPoints);

  //Формируем новые матрицы для расчета производных в точке В=0
  AField[0]:=-MagField_spektr[1];
  AGxx[0]:=GxxExp[1];
  AGxy[0]:=-GxyExp[1];
  for i:=1 to NumberOfPoints do
   begin
    AField[i]:=MagField_spektr[i-1];
    AGxx[i]:=GxxExp[i-1];
    AGxy[i]:=GxyExp[i-1];
   end;

   //Вычисляем производные в точке В=0
   dec(NumberOfPoints);
   GetCoef(AGxx,AField,AField[1],p,p1,A1);
   p1:=0;
   GetCoef(AGxy,AField,AField[1],p,p1,B1);
   B1:=0;

   //Вычисляем производные в точке В=Вмах
   GetCoef(GxxExp,MagField_spektr,MagField_spektr[NumberOfPoints],p,p1,An);
   GetCoef(GxyExp,MagField_spektr,MagField_spektr[NumberOfPoints],p,p1,Bn);
{   An:=0;}
   cs(MagField_spektr,GxxExp,temp_l,A1,An);
   cs(MagField_spektr,GxyExp,temp_t,B1,Bn);

   AddExpPoints(ExpXX,ExpXY);

   Gxx.Clear;
   Gxy.Clear;

   GetLnLimits(Lmin,Lmax);
   SizeData:=(Lmax-Lmin+1)*(PointPerInt-1)*PointPerInt;

   InitArray; 
   k:=0;
   for i:=0 to (lmax-lmin) do
    begin
     lm:=exp((lmin+i)*ln(10));
     sf:=lm;
     for j:=1 to PointPerInt-1 do
      begin
       IntMagField[k]:=sf;
       IntGxx[k]:=sp(MagField_spektr,GxxExp,temp_l,sf);
       IntGxy[k]:=sp(MagField_spektr,GxyExp,temp_t,sf);
       Gxx.AddXY(IntMagField[k],IntGxx[k],'',clTeeColor);
       Gxy.AddXY(IntMagField[k],IntGxy[k],'',clTeeColor);
       sf:=lm*exp(j/PointPerInt*ln(10));
       if sf>MagField_spektr[NumberofPoints] then break;
       inc(k);
      end;
     if sf>MagField_spektr[NumberofPoints] then break;
    end;
   
end;


procedure GetNoise(var a,b:Data_spektr;err:extended);
 var i:word;
     t1,t2,v1,v2:extended;
 begin

   randomize;
   for i:=1 to NumberOfPoints do begin
    t1:=random;
    t2:=random;
    v1:=sqrt(2*sqr(err)*ln(1/t2))*cos(2*pi*t1);
    a[i]:=a[i]*(1+v1);
    v2:=sqrt(2*sqr(err)*ln(1/t2))*cos(2*pi*t1);
    b[i]:=b[i]*(1+v2);
   end;
 end;

procedure MakeMNK( a:boolean;Gxx:TLineSeries;Gxy:TLineSeries;ExpXX,ExpXY:TPointSeries);
var tmp_m:mat;
    coef_t,coef_l:Data_spektr;
    kind,i,j,k,Lmin,Lmax:integer;
    lm,sf:extended;
 function pow(x:extended;y:word):extended;
 var p:extended;
     i:word;
  begin

   p:=1;
   for i:=1 to y do p:=p*x;
   pow:=p;
  end;
procedure BAS(n,M,L:word;X1:extended;var X,T:Data_spektr);
var
  k:word;z,r,denominator:extended;
begin
 //ShowMessage('Inside BAS n='+FloatToStr(n));
 denominator:=(x[n]-x[0]);
 if denominator=0 then
 z:=2.0*(x1-x[0])-1.0         //--------------------АХТУНГ!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else
 z:=2.0*(x1-x[0])/(x[n]-x[0])-1.0;
 t[0]:=1.0;
 t[1]:=z;

 //ShowMessage('BAS part 2');
 for k:=1 to m-1 do
 begin
 r:=z*t[k];
  case L of
    1: r:=r-t[k-1]/4;
    2: r:=2*r-t[k-1];
    3: r:=((k+k+1)*r-k*t[k-1])/(k+1);
  end;
  t[k+1]:=r;
 end;
end;
procedure gram(N,M,L:word;var x,f:Data_spektr; var a:mat);
var i,j,k:integer;
    q,r,s:extended;
    t:Data_spektr;
    p:array of array of extended;
begin
  SetLength(p,MaxPoints,5*MaxPoints);
  SetLength(t,MaxPoints);
  //ShowMessage('Inside gram memory allocated');
 for i:=0 to N do
 begin
   bas(n,m,l,x[i],x,t);
   for j:=0 to m do
     p[j,i]:=t[j]
 end;
 for k:=0 to m do
 begin
  for j:=k to m do
  begin
  s:=0;
  r:=0;
   for i:=0 to n do
   begin
    q:=p[k,i];
    s:=s+q*p[j,i];
     if j=m then
      r:=r+q*f[i];
   end;
   a[k,j]:=s;
   a[j,k]:=s;
  end;
  a[k,m+1]:=r;  // Warning you are not right.
 end;
end;

procedure gauss(N:word;var a:mat; var x:Data_spektr);
 var i,j,k,k1,n1:word;
 r,s:extended;
 begin

   //ShowMessage('in gauss s='+FloatToStr(a[k][k]));

   n1:=N+1;
   for k:=0 to n do
   begin
    k1:=k+1;
    s:=a[k,k];
    //ShowMessage('in gauss s='+FloatToStr(a[k][k]));

    if s = 0 then s:=1;

    for j:=k1 to n1 do
      a[k,j]:=a[k,j]/s;
      for i:=k1 to n do
      begin
          r:=a[i,k];
          for j:=k1 to n1 do
          a[i,j]:=a[i,j]-a[k,j]*r;
      end;
    end;
     for i:=n downto 0 do
     begin
      s:=a[i,n1];
      for j:=i+1 to n do
      s:=s-a[i,j]*x[j];
      x[i]:=s;
    end;
  end;

procedure fi(n,m,l:word;var c,x:Data_spektr;var x1,s:extended);
var i:word;
t:Data_spektr;
begin
 SetLength(t,MaxPoints);
 s:=c[0];
 bas(n,m,l,x1,x,t);
 for i:=1 to m do s:=s+c[i]*t[i]
end;
begin
  //ShowMessage('Inside MakeMNK');
  SetLength(tmp_m,MaxPoints,MaxPoints);
  SetLength(coef_t,MaxPoints);
  SetLength(coef_l,MaxPoints);
  //ShowMessage('Memory allocated');
 {if not(test) then
  begin}
   Power_spektr:=3;  // какая-то степень, или мощность, сильно похоже на кол-во типов носителей
   Kind:=2;  // разновидности, тоже пока не ясно зачем и что
   // эти вызовы заполняют матрицу tmp_m
   // которая кстати уже не пустая О_о
   //ShowMessage('Run gram 1');
   gram(NumberOfPoints,Power_spektr,Kind,MagField_spektr,GxxExp,tmp_m);
   // Гаусс - судя по всему решает матрицу методом Гаусса, сохраняет всё в coef-l
   //ShowMessage('Run gauss 1');
   gauss(Power_spektr,tmp_m,coef_l);
   //ShowMessage('Run gram 2');
   gram(NumberOfPoints,Power_spektr,Kind,MagField_spektr,GxyExp,tmp_m);
   //ShowMessage('Run gauss 2');
   gauss(Power_spektr,tmp_m,coef_t);

   Gxx.Clear; // чистим графики компонент
   Gxy.Clear;
  // ShowMessage('Goodbye gauss and gram!');

   if a then
     begin

   //   ShowMessage('AddExpPoints');
      AddExpPoints(ExpXX,ExpXY); // добавляет точки на график, экспериментальные
   //   ShowMessage('GetLnLimits');
      GetLnLimits(Lmin,Lmax); // получает пределы, логарифмические
      SizeData:=(Lmax-Lmin+1)*(PointPerInt-1)*PointPerInt; // считаем размер данных
      InitArray;  // выделяем его---------------------------------------------------------------
   //   ShowMessage('After init array in MakeMNK');
      k:=0;
      for i:=0 to (lmax-lmin) do
       begin
       lm:=exp((lmin+i)*ln(10));
       sf:=lm;
       for j:=1 to PointPerInt-1 do
        begin
        IntMagField[k]:=sf;
{        IntGxxExp[k]:=fi(Power,coef_l,sf);
        IntGxyExp[k]:=fi(Power,coef_t,sf);}
        // а тут происходит самое страшное - считаются два основных графика
        fi(NumberOfPoints,Power_spektr,Kind,coef_l,MagField_spektr,sf,IntGxx[k]);
        fi(NumberOfPoints,Power_spektr,Kind,coef_t,MagField_spektr,sf,IntGxy[k]);
        // и судя по всему ось х - действительно подвижность
        // а ось у - это величины компонент тензора проводимости
        // но они как-то модифицированы
        // видимо в предыдущей функции, надо уточнить
        Gxx.AddXY(IntMagField[k],IntGxx[k],'',clTeeColor);
        Gxy.AddXY(IntMagField[k],IntGxy[k],'',clTeeColor);
        sf:=lm*exp(j/PointPerInt*ln(10));
        if sf>MagField_spektr[NumberofPoints] then break;
        inc(k);
       end;
      if sf>MagField_spektr[NumberofPoints] then break;
     end;

    end else
    begin
     for i:=0 to NumberOfPoints do
      begin
        fi(NumberOfPoints,Power_spektr,Kind,coef_l,MagField_spektr,sf,GxxExp[i]);
        fi(NumberOfPoints,Power_spektr,Kind,coef_t,MagField_spektr,sf,GxyExp[i]);
      end;
    end;

end;

procedure MakeLagranj;
var X,Y,Y1,Y2:extended;
begin
  X:=MagField_spektr[NumberOfPoints]-(MagField_spektr[NumberOfPoints]-
        MagField_spektr[NumberOfPoints-1])/4;
  GetCoef(GxxExp,MagField_spektr,X,Y,Y1,Y2);
  GxxExp[NumberOfPoints+1]:=GxxExp[NumberOfPoints];
  GxxExp[NumberOfPoints]:=Y;
  GetCoef(GxyExp,MagField_spektr,X,Y,Y1,Y2);
  GxyExp[NumberOfPoints+1]:=GxyExp[NumberOfPoints];
  GxyExp[NumberOfPoints]:=Y;
  MagField_spektr[NumberOfPoints+1]:=MagField_spektr[NumberOfPoints];
  MagField_spektr[NumberOfPoints]:=X;
  inc(NumberOfPoints);
end;



procedure Tred2(n:word; tol:extended; var d,e:dat1;
                 var a,z:dat2; var fail:boolean);
{label Skip;}
var
  i,j,k,l  :integer;
  f,g,h,hh :extended;

begin
fail:=false;
for i:=1 to n do
  for j:=1 to i do z[i,j]:=a[i,j];

for i:=n downto 2 do     { ? }
begin
  l:=i-2; f:=z[i,i-1]; g:=0;
  for k:=1 to l do g:=g+sqr(z[i,k]);
{ if g is small and do not sure of orthogonalization then procedure fails }
  h:=g+f*f;
{  if g<=tol then
    begin fail:=true; e[i]:=f; h:=0.; goto Skip end;}

  l:=l+1;
  if f>=0. then begin e[i]:=-sqrt(h); g:=-sqrt(h) end
    else begin e[i]:=sqrt(h); g:=sqrt(h) end;
  h:=h-f*g; z[i,i-1]:=f-g; f:=0.;

  for j:=1 to l do
  begin
    z[j,i]:=z[i,j]/h; g:=0.;
    for k:=1 to j do g:=g+z[j,k]*z[i,k];
    for k:=j+1 to l do g:=g+z[k,j]*z[i,k];
    e[j]:=g/h; f:=f+g*z[j,i]
  end;

  hh:=f/(h+h);
{ matrix reduce }
  for j:=1 to l do
  begin
    f:=z[i,j]; g:=e[j]-hh*f; e[j]:=e[j]-hh*f;
    for k:=1 to j do
      z[j,k]:=z[j,k]-f*e[k]-g*z[i,k]
  end;

{Skip:}
  d[i]:=h
end; { i }

  d[1]:=0; e[1]:=0;
{ store transformation matrix }
  for i:=1 to n do
  begin
    l:=i-1;
    if d[i]<>0. then
      for j:=1 to l do
      begin
        g:=0.;
        for k:=1 to l do g:=g+z[i,k]*z[k,j];
        for k:=1 to l do z[k,j]:=z[k,j]-g*z[k,i]
      end;
      d[i]:=z[i,i]; z[i,i]:=1;
      for j:=1 to l do begin z[i,j]:=0; z[j,i]:=0 end;
  end;
end; { Tred2 }


{ **************************************************************** }

{ J.H.Wilkinson, C.Reinch. Handbook for Automatic Computation.
  Linear Algebra. Springer Verlag: Heidelberg, NewYork, Berlin.
  algorithm II-4.
  Procedure calculate all eigenvalues and eigenvectors
  of real tridiagonal symmetric matrix by QL-algorithm
  with implicit shift.

  Data:
  n - matrix order;
  macheps - smallest possible value such thet 1+macheps>1;
  d - array [1..n] - containts diagonal elements of matrix;
  e - array [1..n] - containts underdiagonal elements of matrix;
   e[1] - arbitraray;
  z - array [1..n,1..n] - unit matrix or matrix of reduction
   to tridiagonal form by Tred2.

  Results:
  d - eigenvalues in increasing order;
  z - eigenvectors (z[i,k],k=1..n);
  e - information is losed.

  Test:
  d=(1,1e2,1e4,1e6,1e8,1e10,1e12);
  e=(0,10,1e3,1e5,1e7,1e9,1e11);
   eigenvalues:         eigenvector for min eigenvalue

  -9.46347415648e8      9.99949513816e-1
  -9.46346919727e2     -1.00974711993e-5
   9.99899020189e-1    -9.99849548745e-3
   1.04633721478e3      9.99850548602e-4
   1.00989903020e6     -1.00974710982e-14
   1.04633771269e9     -9.99850548500e-6
   1.01000000980e12     9.99850548498e-7

  d=(1e12,1e10,1e8,1e6,1e4,1e2,1);
  e=(0,1e11,1e9,1e7,1e5,1e3,10);
   eigenvalues:         eigenvector for min eigenvalue

  -9.46347415707e8      9.99850548432e-7
  -9.46346919876e2     -9.99850548417e-6
   9.99899020157e-1    -1.12526974437e-14
   1.04633721466e3      9.99850548546e-4
   1.00989903020e6     -9.99849548703e-3
   1.04633771265e9     -1.00974740462e-5
   1.01000000980e12     9.99949513813e-1                   }

procedure Imtql2(n:word; macheps:extended; var d,e:dat1;
                 var z:dat2; var fail:boolean);
label Test,Cont,Fail_exit;
var
  i,ia,j,k,m,its  :integer;
  h,c,p,q,s,t,u   :extended;
begin
fail:=false;
for i:=2 to n do e[i-1]:=e[i];
e[n]:=0.; k:=n-1;
for j:=1 to n do
begin
  its:=0;

{ searching of negligibly small underdiagonal element}
Test:
  for m:=j to k do
    if abs(e[m])<=macheps*(abs(d[m])+abs(d[m+1])) then goto Cont;
  m:=n;

Cont:
  u:=d[j];
  if m<>j then
  begin
    if its=30 then begin fail:=true; goto Fail_exit end;
    its:=its+1;
{ formation of shift }
    q:=(d[j+1]-u)/(2*e[j]);
    t:=sqrt(1.+sqr(q));
    if q<0. then q:=d[m]-u+e[j]/(q-t) else q:=d[m]-u+e[j]/(q+t);
    u:=0.; s:=1.; c:=1.;

    for i:=m-1 downto j do
    begin
      p:=s*e[i]; h:=c*e[i];
      if abs(p)>=abs(q)
      then
      begin
        c:=q/p; t:=sqrt(sqr(c)+1.0);
        e[i+1]:=p*t; s:=1./t; c:=c*s;
      end
      else
      begin
        s:=p/q; t:=sqrt(sqr(s)+1.0);
        e[i+1]:=q*t; c:=1./t; s:=s*c
      end;

      q:=d[i+1]-u; t:=(d[i]-q)*s+2.*c*h;
      u:=s*t; d[i+1]:=q+u; q:=c*t-h;

{ calculation of eigenvector }
      for ia:=1 to n do
      begin
        p:=z[ia,i+1];
        z[ia,i+1]:=s*z[ia,i]+c*p;
        z[ia,i]:=c*z[ia,i]-s*p
      end { ia }
    end { i };
    d[j]:=d[j]-u; e[j]:=q; e[m]:=0.; goto Test
  end { m<>j }
end; { j }

{ sorting of eigenvalues and eigenvectors }
for i:=1 to n do
begin
  k:=i; p:=d[i];
  for j:=i+1 to n do
    if d[j]<p then begin k:=j; p:=d[j] end;
  if k<>i then
  begin
    d[k]:=d[i]; d[i]:=p;
    for j:=1 to n do
    begin p:=z[j,i]; z[j,i]:=z[j,k]; z[j,k]:=p end
  end
end; { i }
Fail_exit:
end; { Imtql2 }

function GetElem(j1,k1,i1:word):extended;
 var
  s:extended;
  ii:word;
 begin
    s:=0;
    for  ii:=i1 to (NumberOfPoints-(j1-2)) do
    begin
     if j1=2 then
     begin

       if ii<>k1 then
          s:=s+sqr(B_spektr[ii])
     end
     else
        if ii<>k1 then
        s:=s+sqr(B_spektr[ii])*GetElem(j1-1,k1,ii+1);
     end;
    GetElem:=s;
 end;


procedure MakeMatrC;
var
 j,k:word;
begin
 {fillchar(Cr, Sizeof(Cr), 0);
 ShowMessage('after Filling');
 Cr[1,1]:=12;
 Cr[11,11]:=14;
 fillchar(Cl, Sizeof(Cl), 0);
 fillchar(Cr_t,Sizeof(Cr_t),0);
 fillchar(Cl_t,Sizeof(Cl_t),0);
 fillchar(Am, Sizeof(Am), 0); }  // заполнение массивов нулями как я погляжу...
 //ShowMessage('Cr is'+FloatToStr(Cr[5,5]));

 //ShowMessage('End up with fillings, Number of points is'+ IntToStr(NumberofPoints));
 for j:=1 to NumberOfPoints do
  for k:=1 to NumberOfPoints do
  begin
   if j=1 then
   begin
    // ShowMessage('Step here j='+IntToStr(j)+' k='+IntToStr(k));
    Cr[j,k]:=1;
   end  

   else
   begin
     //ShowMessage('One step close');
     Cr[j,k]:=GetElem(j,k,1);
   end;  

     //ShowMessage('Are you there?');
     Cl[j,k]:=-Cr[j,k]*B_spektr[k]; 
    end;
end;

procedure MakeMatrA;
 var i,j,k:word;
 begin
  //fillchar(Am,sizeof(Am),0);
  //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  for i:=1 to NumberOfPoints do
   for j:=1 to NumberOfPoints do
    if odd(i+j) then
      for k:=1 to NumberOfPoints do
           Am[i,j]:=Am[i,j]+Gxy_sp[k]*Cl_t[k,(i+j-1) shr 1]
               else
       for k:=1 to NumberOfPoints do
           Am[i,j]:=Am[i,j]+Gxx_sp[k]*Cr_t[k,(i+j) shr 1]
 end;

  {if su=0 then det(C)=0}
procedure InverseMatrC(var Ci:dat2;var C:dat2;var Su:extended;NP:word);
  var i,j,k:word;
       at:dat3;
       sr:extended;
  begin
    SetLength(at,MaxPoints+1,2*MaxPoints+1);
  for i:=1 to Np do
   for j:=1 to NP do
      at[i,j]:=Ci[i,j];
   for i:=1 to NP do begin
    for j:=NP+1 to 2*NP do At[i,j]:=0; At[i,i+NP]:=1;
   end;
   for k:=1 to NP do begin Su:=At[k,k]; j:=k;
    for i:=k+1 to NP  do begin Sr:=At[i,k];
     if abs(Sr)>abs(Su) then begin su:=sr; j:=i; end;
    end;
    if su=0 then exit;
    if j<>k then for i:=k to 2*NP do begin
     sr:=at[k,i];at[k,i]:=at[j,i];at[j,i]:=sr;end;
    for j:=k+1 to 2*NP do at[k,j]:=at[k,j]/su;
    for i:=k+1 to NP do begin sr:=at[i,k];
    for j:=k+1 to 2*NP do at[i,j]:=at[i,j]-at[k,j]*sr
   end;
  end;
  if su<>0 then
   for j:=NP+1 to 2*NP do
    for i:=NP-1 downto 1 do begin sr:=at[i,j];
    for k:=i+1 to NP do sr:=sr-at[k,j]*at[i,k];
    at[i,j]:=sr;
   end;
   if su<>0 then
    for i:=1 to Np do
     for j:=1 to NP do
      C[i,j]:=at[i,j+NP];
  end;

function S_s(Mi:extended):extended;
   var im,js:word;
       a,bb:extended;
   begin
     { make Vm }
     Mv[1]:=1;
     if abs(Mi)<1e-26 then
      for im:=2 to NumberOfPoints do Mv[im]:=0 else
     for im:=2 to NumberOfPoints do
      if odd(im) then
       Mv[im]:=exp((im-1)*ln(abs(Mi))) else
       if Mi<0 then Mv[im]:=exp((im-1)*ln(abs(Mi))) else
          Mv[im]:=-exp((im-1)*ln(abs(Mi)));
(*     { compute ¦Vm¦¤}
     vm:=0;
     for is:=1 to NPoint do
       vm:=vm+sqr(Mv[is]);
     { compute am/Go for Hp=0}
     a:=0;
     for is:=1 to NPoint do a:=a+sqr(Mv[is]);*)
     a:=1;
     { compute am/Go }
     for im:=1 to NumberOfPoints do a:=a*(1+sqr(B_spektr[im]*Mi));
     for im:=1 to NumberOfPoints do begin Vpr[im]:=0;
     for js:=1 to NumberOfPoints do
          Vpr[im]:=Vpr[im]+Qm[js,im]*Mv[js];
      Vpr[im]:=sqr(Vpr[im]);end;
     bb:=0;
     for im:=1 to NumberOfPoints do begin
{      if (abs(Lv[im])<1e-30)and(Vpr[im]<1e-15) then continue}
       if (abs(Lv[im])=0)and(Vpr[im]>1e-15) then begin a:=0;break;end;
      bb:=bb+Vpr[im]/abs(Lv[im]); end;
     S_s:=(a/bb);
   end;

//////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
procedure MobilitySpectrumFunc(LineSeries1,Series5:TLineSeries);
var
  Sf,lm:extended;
  j,k,i:word;
  Lmin,Lmax:integer;
{  LogAxis,NoLogAxis:boolean;}


begin
 // ShowMessage('Inside MobilitySpectrum');
   for i:=0 to NumberOfPoints do
    begin
     B_spektr[i+1]:=MagField_spektr[i];
     Gxx_sp[i+1]:=GxxExp[i];
     Gxy_sp[i+1]:=GxyExp[i];
    end;
   inc(NumberOfPoints);
  // ShowMessage('Run MakeMatrC');
   MakeMatrC;
  // ShowMessage('RunInverseMatcC');
   InverseMatrC(Cr,Cr_t,Sf,NumberOfPoints);
   if Sf=0 then begin MessageDLG('Определитель равен нулю!',mtError,[mbOK],0);
         dec(NumberOfPoints);
         exit; end;
   if MagField_spektr[0]=0 then
   begin
    for j:=1 to NumberOfPoints-1 do
    for k:=1 to NumberOfPoints-1 do
     Cm[j,k]:=Cl[j,k+1];
     InverseMatrC(Cm,Cm_t,Sf,NumberOfPoints-1);
      for j:=1 to NumberOfPoints do
       for k:=1 to NumberOfPoints do
        if j=1 then Cl_t[1,k]:=0 else
         if k=NumberOfPoints then Cl_t[j,NumberOfPoints]:=0 else
          Cl_t[j,k]:=Cm_t[j-1,k];
    end else
    InverseMatrC(Cl,Cl_t,Sf,NumberOfPoints);

    //ShowMessage('Run MakeMartA. Line 891');
    MakeMatrA;
    // интересно что это за числа ниже?
    //ShowMessage('Run Tred2');
    Tred2(NumberOfPoints,6e-4913,Lv,Xv,am,Qm,bulua);
    //ShowMessage('Run Imtql2');
    Imtql2(NumberOfPoints,5.42e-20,Lv,Xv,Qm,bulua);

     LineSeries1.Clear;
     Series5.Clear;
   Lmin:=MSLeft;LMax:=MSRight;
   SizeData:=(Lmax-Lmin+1)*(PointPerInt-1)*PointPerInt; // возможно тут придется домножать на размер элемента
   // впрочем нет - там же сейчас другая функция.
  // ShowMessage('InitArray2');
   InitArray2;
   k:=0;
  // ShowMessage('Last for...');
   for i:=0 to (lmax-lmin) do
    begin
     lm:=exp((lmin+i)*ln(10));
     sf:=lm;
     for j:=1 to PointPerInt-1 do
      begin
       Mobility[k]:=sf;
       Spectr_e[k]:=S_s(-sf);
       Spectr_p[k]:=S_s(sf);
       LineSeries1.AddXY(Mobility[k],Spectr_e[k],'',clTeeColor);
       Series5.AddXY(Mobility[k],Spectr_p[k],'',clTeeColor);
       sf:=lm*exp(j/PointPerInt*ln(10));
       inc(k);
       if sf>10 then break;
      end;
     if sf>10 then break;
    end;
    GridPoints:=k-1;
   //AddPoints2;
   dec(NumberOfPoints);
end;


function GxxExpi:extended;
var i:word;
begin
  Result:=0;
  for i:=0 to gridPoints do
   result:=result+Axx^[i]*(QSpectr_p^[i]+QSpectr_e^[i]);
end;


    {
      Итак, по сути, здесь происходит тыканье значений максимумов и
      отображение их в таблице под графиком спектра подвижности.
      Т.е. по факту этот код будет уже на сях и здесь он не нужен.
      В комментарии его.


      }
procedure chtSpectrClickSeries(chtSpectr: TCustomChart;
  Series: TChartSeries; ValueIndex: Integer;
  Shift: TShiftState; X, Y: Integer; Series5:TCustomSeries;
  LineSeries1:TLineSeries; StringGrid3:TStringGrid;RowInFocus:Integer);
//var Mu,G_e,G_p,con_p,con_e:extended;
 //     Ind, Ind2:longint;


begin
    {
with chtSpectr do
   begin

    ind:=Series5.GetCursorValueIndex;
    Ind2:= LineSeries1.GetCursorValueIndex;
    if ind<>-1 then
      begin                 // это по дыркам.
       G_p:=Series5.YValue[ind]; // наверное проводимость
       Mu:=Series5.XValue[ind];  // подвижность
       con_p:=g_p/(mu*1.602e-19); // концентрация
       StringGrid3.Cells[1,RowInFocus]:=FloatToStr(con_p); // показываем выбранные результаты
       StringGrid3.Cells[2,RowInFocus]:=FloatToStr(Mu);
       if RowInFocus<3 then
         Inc(RowInFocus) // и да, переключаем фокус таблицы
       else
         RowInFocus:=1;
      end;
      if ind2<>-1 then
      begin
       G_e:=LineSeries1.YValue[ind2]; // это по электронам
       Mu:=Series5.XValue[ind2];
       con_e:=g_e/(mu*1.602e-19);
       StringGrid3.Cells[1,RowInFocus]:=FloatToStr(-con_e);
       StringGrid3.Cells[2,RowInFocus]:=FloatToStr(-Mu);
       if RowInFocus<3 then
         Inc(RowInFocus)
       else
         RowInFocus:=1;
      end;
   end;   }

end;


   // управление фокусом таблицы, особо без надобности.
   {
procedure SetGridFocus(SGrid,StringGrid3: TStringGrid; r:integer);
var
 Cell_Focus:TGridRect;

begin
  Cell_Focus.Top:=r;
  Cell_Focus.Bottom:=r;
  Cell_Focus.Left:=0;
  Cell_Focus.Right:=3;
  StringGrid3.Selection:=Cell_Focus;
end;       }
    {
procedure chtSpectrMouseMove(Sender: TObject; Shift: TShiftState; X,
  Y,RowInFocus: Integer;StringGrid3,SGrid: TStringGrid);
begin
  SetGridFocus(StringGrid3,SGrid,RowInFocus);  // Ахтунг! Странное использование функции!
end;
        }


procedure LoadSpektrResults;//Загружает результаты Спектра подвижности
//var i:Integer;
//   con, mob:Extended;
begin
  {
  con:=0;
  mob:=0;
{$I-}
   { Reset(Data_File);
{$I+}
{
  if IOResult<>0 then
   begin
    ShowMessage('Ошибка открытия файла');
    Exit;
   end
  else
 begin
  try
    Readln(Data_File, NumberOfPoints);
    Readln(Data_File,Max_Value[1],Max_Value[2],Max_Value[3],Max_Value[4],
      Max_Value[5],Max_Value[6]);
    Readln(Data_File,Min_Value[1],Min_Value[2],Min_Value[3],Min_Value[4],
      Min_Value[5],Min_Value[6]);
    

    for i:=1 to StringGrid3.RowCount-1 do
    begin
      Readln(Data_file,con,mob);
      StringGrid3.cells[1,i]:=FloatToStr(con);
      StringGrid3.cells[2,i]:=FloatToStr(mob);
    end;

    
    Label22.Caption:=OpenDialog2.FileName;
    MakeMNK(true);
    MobilitySpectrumFunc;
  except
    ShowMessage('Ошибка загрузки данных');
    Rewrite(Data_File);
  end;
  CloseFile(Data_File);
 end;}
end;


// Загрузка результатов спектра подвижности
// Тут от неё толку нет.
{
procedure btnSpectrResultClick(Sender: TObject);
begin
  
  if ODSpektrRes.Execute then
  begin
     AssignFile(Data_File, ODSpektrRes.FileName);
     ReplacementOfSeparator;
     LoadSpektrResults;
  end;
end;  }

// Запись результатов Спектра подвижности
// Количество точек
// Минимумы и максимумы по носителям зарядов
// И дальше неизвестно что:)
// Но нам надо то, что оно выводило на графики.
procedure WriteSpektrResults;
//var i:Integer;
begin
  {
{$I-}
 { Rewrite(Data_File);
{$I+}
 { if IOResult<>0 then
   begin
    ShowMessage('Ошибка сохранения файла');
    Exit;
   end
  else
  begin
    Writeln(Data_File, NumberOfPoints+1);
    Writeln(Data_File,Max_Value[1],'	',Max_Value[2],'	',Max_Value[3],'	',
    Max_Value[4],'	',Max_Value[5],'	',Max_Value[6]);
    Writeln(Data_File,Min_Value[1],'	',Min_Value[2],'	',Min_Value[3],'	',
    Min_Value[4],'	',Min_Value[5],'	',Min_Value[6]);
    
          // GxxExp и GxxExp в спектре подвижности играют роль экспериментальных
          // данных.
    for i:=1 to StringGrid3.RowCount-1 do
    Writeln(Data_file,StringGrid3.cells[1,i]+'	'+StringGrid3.cells[2,i]);
    CloseFile(Data_File);
  end;
  }
end;


////////////////////////////////////////////////////////////////////////////
/////////////////////// КОНЕЦ "ХОЛЛ. ПОДВИЖНОСТЬ"///////////////////////////
////////////////////////////////////////////////////////////////////////////

function getResultEY(i:Integer):Extended; stdcall;
begin
  Result:=electronMobilitySpectrum.YValue[i];
end;

function getResultEX(i:Integer):Extended; stdcall;
begin
  Result:=electronMobilitySpectrum.XValue[i];
end;

function getResultHY(i:Integer):Extended; stdcall;
begin
  Result:=holeMobilitySpectrum.YValue[i];
end;

function getResultHX(i:Integer):Extended; stdcall;
begin
  Result:=holeMobilitySpectrum.XValue[i];
end;

function getResults(var outElectronSpectrX,  outElectronSpectrY,
outHoleSpectrX,outHoleSpectrY: Data_spektr): Pointer; stdcall;
var i:Integer;
begin
  ShowMessage('Inside getResults. Count is '+ IntToStr(electronMobilitySpectrum.XValues.Count));

  SetLength(outSpectrEX,electronMobilitySpectrum.XValues.Count);
  SetLength(outSpectrEY,electronMobilitySpectrum.XValues.Count);
  SetLength(outSpectrHX,electronMobilitySpectrum.XValues.Count);
  SetLength(outSpectrHY,electronMobilitySpectrum.XValues.Count);
  // Краш на 48 на 4м
  for i:=0 to electronMobilitySpectrum.XValues.Count-1 do
  begin
   // ShowMessage(IntToStr(i));
   // ShowMessage('1');
  outSpectrEX[i]:=electronMobilitySpectrum.XValue[i];
  //ShowMessage('2');
  outSpectrEY[i]:=electronMobilitySpectrum.YValue[i];
  //ShowMessage('3');
  outSpectrHX[i]:=holeMobilitySpectrum.XValue[i];
  //ShowMessage('4');
  //ShowMessage(IntToStr(Length(outSpectrHY)));
  outSpectrHY[i]:=holeMobilitySpectrum.YValue[i];
  end;
  Result:=@outSpectrEX;
  {
  for i:=0 to electronMobilitySpectrum.XValues.Count do
  begin
    ShowMessage('i='+IntToStr(i));

    outElectronSpectrX[i]:=electronMobilitySpectrum.XValue[i];
    ShowMessage('For Y');
    outElectronSpectrY[i]:=electronMobilitySpectrum.YValue[i];
  end;
  ShowMessage('Netx Part of getResults');
  for i:=0 to holeMobilitySpectrum.XValues.Count do
  begin
    ShowMessage('write x');
   outHoleSpectrX[i]:=holeMobilitySpectrum.XValue[i];
   ShowMessage('write y');
   outHoleSpectrY[i]:=holeMobilitySpectrum.YValue[i];
  end;
   }
end;  


// входные параметры есть, надо подумать что мы возвращаем, где оно хранится
// и как мы это будем возвращать.
function RunMobilitySpectrum (MagneticFieldP,Exx,Exy: Data_spektr; size:Integer ):Integer; stdcall;
var i:Integer;
Gxx,Gxy:TLineSeries;
ExpXX,ExpXY:TPointSeries;
//t:Int64;
begin

  MaxPoints:=size;
  ShowMessage(IntToStr(MaxPoints));
  if MaxPoints>1000 then
  MaxPoints:=11;

  NumberOfPoints:=MaxPoints-1;

  //ShowMessage(IntToStr(MaxPoints));

  SetLength(MagField_spektr,MaxPoints);
  SetLength(GxxExp,MaxPoints);
  SetLength(GxyExp,MaxPoints);

                       // тут бы вообще другой размер нужен (см. ниже).
  {
  SetLength(IntGxx,(PointPerInt-1)*PointPerInt);
  SetLength(IntGxy,(PointPerInt-1)*PointPerInt);
  SetLength(IntMagField,(PointPerInt-1)*PointPerInt);
  
  SetLength(Spectr_e,(PointPerInt-1)*PointPerInt);
  SetLength(Spectr_p,(PointPerInt-1)*PointPerInt);
  SetLength(Mobility,(PointPerInt-1)*PointPerInt);
  }

  // где они используются вообще?О_о
  {SetLength(QSpectr_e,(PointPerInt-1)*PointPerInt);
  SetLength(QSpectr_p,(PointPerInt-1)*PointPerInt);
  SetLength(Axx,(PointPerInt-1)*PointPerInt);
  SetLength(Axy,(PointPerInt-1)*PointPerInt);
  SetLength(Axx_d,(PointPerInt-1)*PointPerInt);
  SetLength(Axy_d,(PointPerInt-1)*PointPerInt);
  SetLength(dIntGxx,(PointPerInt-1)*PointPerInt);
  SetLength(dIntGxy,(PointPerInt-1)*PointPerInt);
  SetLength(QGxx,(PointPerInt-1)*PointPerInt);
  SetLength(QGxy,(PointPerInt-1)*PointPerInt);
  SetLength(dQGxx,(PointPerInt-1)*PointPerInt);
  SetLength(dQGxy,(PointPerInt-1)*PointPerInt);  }


  SetLength(B_spektr,MaxPoints+1);

  SetLength(Gxx_sp,MaxPoints+1);
  SetLength(Gxx_MC,MaxPoints+1);
  SetLength(Gxy_MC,MaxPoints+1);
  SetLength(Gxy_sp,MaxPoints+1);
  SetLength(Xr,MaxPoints+1);
  SetLength(Lv,MaxPoints+1);
  SetLength(Xv,MaxPoints+1);
  SetLength(Mv,MaxPoints+1);
  SetLength(Vpr,MaxPoints+1);



  SetLength(Am,MaxPoints+1,MaxPoints+1);
  SetLength(Qm,MaxPoints+1,MaxPoints+1);
  SetLength(Cl,MaxPoints+1,MaxPoints+1);
  SetLength(Cr,MaxPoints+1,MaxPoints+1);
  
  SetLength(Cl_t,MaxPoints+1,MaxPoints+1);
  SetLength(Cr_t,MaxPoints+1,MaxPoints+1);
  SetLength(Cm,MaxPoints+1,MaxPoints+1);
  SetLength(Cm_t,MaxPoints+1,MaxPoints+1);


   {
  mat=array[0..MaxPoints,0..MaxPoints] of extended;
  Dat1=array[1..MaxPoints] of extended;
  Dat2=array[1..MaxPoints,1..MaxPoints] of extended;
  Dat3=array[1..MaxPoints,1..2*MaxPoints] of extended;

   }


  //ShowMessage('Memory allocated successfully');
  
  for i:=0 to MaxPoints-1 do
  begin
    MagField_spektr[i]:=MagneticFieldP[i];
    GxxExp[i]:=Exx[i];
    GxyExp[i]:=Exy[i];
  end;

  //ShowMessage('Data has copied');

  Gxx:=TLineSeries.Create(nil);
  Gxy:=TLineSeries.Create(nil);
  ExpXX:=TPointSeries.Create(nil);
  ExpXY:=TPointSeries.Create(nil);
  electronMobilitySpectrum:=TLineSeries.Create(nil);
  holeMobilitySpectrum:=TLineSeries.Create(nil);

  //ShowMessage('Run MakeMNK');

  MakeMNK(true,Gxx,Gxy,ExpXX,ExpXY);
  ShowMessage('Run MobilitySpectrumFunc');

  MobilitySpectrumFunc(electronMobilitySpectrum,holeMobilitySpectrum);


  // Прибираемся
  ShowMessage('Run clean up');
  SetLength(MagField_spektr,0);
  SetLength(GxxExp,0);
  SetLength(GxyExp,0);

  // где они используются вообще?О_о
  {SetLength(QSpectr_e,0);
  SetLength(QSpectr_p,0);
  SetLength(Axx,0);
  SetLength(Axy,0);
  SetLength(Axx_d,0);
  SetLength(Axy_d,0);
  SetLength(dIntGxx,0);
  SetLength(dIntGxy,0);
  SetLength(QGxx,0);
  SetLength(QGxy,0);
  SetLength(dQGxx,0);
  SetLength(dQGxy,0);
   }

  SetLength(B_spektr,0);
  SetLength(Gxx_sp,0);
  SetLength(Gxx_MC,0);
  SetLength(Gxy_MC,0);
  SetLength(Gxy_sp,0);
  SetLength(Xr,0);
  SetLength(Lv,0);
  SetLength(Xv,0);
  SetLength(Mv,0);
  SetLength(Vpr,0);

  SetLength(Am,0,0);
  SetLength(Qm,0,0);
  SetLength(Cl,0,0);
  SetLength(Cr,0,0);
  SetLength(Cl_t,0,0);
  SetLength(Cr_t,0,0);
  SetLength(Cm,0,0);
  SetLength(Cm_t,0,0);

  Result:=electronMobilitySpectrum.XValues.Count;
end;  

exports

RunMobilitySpectrum,getResults,getResultEY,getResultEX,getResultHY,getResultHX;

begin

     DllProc:=@ProcDLL;

     // Картина такая:
     // Надо выделить память.

   //SetLength(IntMagfield,SizeData);
   //SetLength(IntGxx,SizeData);
   //SetLength(IntGxy,SizeData);

   //SetLength(Spectr_e,SizeData);
   //SetLength(Spectr_P,SizeData);
   //SetLength(Mobility,SizeData);


end.

