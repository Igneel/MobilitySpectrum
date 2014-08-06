library MobilitySpectrum;
uses
  Windows,  SysUtils, Variants, Classes, Graphics, Controls, Forms,
  Grids, Clipbrd, Buttons, StdCtrls, StrUtils, ExtCtrls, Menus,
  CheckLst, Math, ComCtrls, Gauges,
  Series,Chart,TeEngine, TeeProcs, Messages,Dialogs;

  // optim1

{$R *.res}
const
  MaxPoints=100;
  PointPerInt=50;
  //MaxParameters=8;    // Максимальное число параметров


  {  Mobility                         [ m**2/V*s ]     ;
      Magnetic Field                   [ Tesla ]        ;
      Electric Field                   [ V/m ]          ;
      Surface Recombination Velocity   [ m/s ]          ;
      All other values in SI                             ;   }
  type
    Data_spektr=array [0..MaxPoints] of extended;
    PData_spektr=^Data_spektr;

    ImageDat=array [0..(PointPerInt-1)*PointPerInt] of extended;
  PImageDat=^ImageDat;  // указатель на массив из двух элементов.
  mat=array[0..MaxPoints,0..MaxPoints] of extended;
  Dat1=array[1..MaxPoints] of extended;
  Dat2=array[1..MaxPoints,1..MaxPoints] of extended;
  Dat3=array[1..MaxPoints,1..2*MaxPoints] of extended;



  var
    /////////////// спектр подвижности///////////
  NumberofPoints,Power_spektr,GridPoints :word;
  MagField_spektr,Gxx,Gxy,GxxExp,GxyExp  :Data_spektr;
  IntGxx,IntGxy,IntMagField,Spectr_e,Spectr_p,Mobility,
  QSpectr_e,QSpectr_p,Axx,Axy,Axx_d,Axy_d,
  dIntGxx,dIntGxy,QGxx,QGxy,dQGxx,dQGxy  :PImageDat;
  //Peak_e,Peak_p,Valley_e,Valley_p        :PImageDat3;
  SizeData                               :longint;
  B_spektr,Gxx_sp,Gxx_MC,Gxy_MC,
  Gxy_sp,Xr,Lv,Xv,Mv,Vpr                 :Dat1;
  Am,Qm,Cl,Cr,Cl_t,Cr_t,Cm,Cm_t          :Dat2;
  bulua                                  :boolean;
  MSRight,MSLeft                         :integer;
  Min_Spectr,Coef1,Coef2,Mu_min,Mu_max,
  W,F_s,A1,An,B1,Bn, Ves1, Ves2          :extended;

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




  // Выделение памяти нужно будет перенести в инициализацию разумеется.
procedure InitArray;
begin
   if IntMagField<>nil then freemem(IntMagField);
   getmem(IntMagfield,SizeData);
   if IntGxx<>nil then freemem(IntGxx);
   getmem(IntGxx,SizeData);
   if IntGxy<>nil then freemem(IntGxy);
   getmem(IntGxy,SizeData);
end;

procedure InitArray2;
begin
  if Spectr_e<>nil then freemem(Spectr_e);
   getmem(Spectr_e,SizeData);
   if Spectr_P<>nil then freemem(Spectr_P);
   getmem(Spectr_P,SizeData);
   if Mobility<>nil then freemem(Mobility);
   getmem(Mobility,SizeData);
end;

// В общем это освобождение памяти, и его надо бы перетащить в другое место.
// И разумеется это старьё надо переделать на нормальные динамические массивы:).
procedure FormDestroy(Sender: TObject);
begin
if IntMagField<>nil then freemem(IntMagField);
   if IntGxx<>nil then freemem(IntGxx);
   if IntGxy<>nil then freemem(IntGxy);
   if Spectr_e<>nil then freemem(Spectr_e);
   if Spectr_P<>nil then freemem(Spectr_P);
   if Mobility<>nil then freemem(Mobility);
   if Axx<>nil then freemem(Axx);
   if Axy<>nil then freemem(Axy);
   if Axx_d<>nil then freemem(Axx_d);
   if Axy_d<>nil then freemem(Axy_d);
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


 // Вывод точек на график. Пока остаётся закоментированным.

procedure AddPoints2;
//var sf,b1,b2:extended;

begin
  {
   with chtSpectr do    // график спектра как я посмотрю
    begin
      with BottomAxis do // его нижние оси
       begin
        Logarithmic:=true;  // настройка масштаба
        Automatic:=false;
        sf:=0.05*(Series[0].MaxXValue-Series[0].MinXValue);
        if (Series[0].MinXValue-sf)>0 then Minimum:=Series[0].MinXValue-sf
            else Minimum:=0.0005;
       Maximum:=Series[0].MaxXValue+sf;
      end;
     with LeftAxis do
      begin
       Logarithmic:=True;
       Automatic:=False;
       if Series[0].MaxYValue>Series[1].MaxYValue then
         b1:=Series[0].MaxYValue else b1:=Series[1].MaxYValue;
       if Series[0].MinYValue<Series[1].MinYValue then
         b2:=Series[0].MinYvalue else b2:=Series[1].MinYvalue;
        sf:=0.05*(b1-b2);
       if (b2-sf)>Maximum then
        begin
         Maximum:=b1+sf;
         if (b2-sf)>0 then
                       Minimum:=b2-sf
                      else
                       Minimum:=b2*0.95;

        end
        else begin
         if (b2-sf)>0 then
                       Minimum:=b2-sf
                      else
                       Minimum:=b2*0.95;
         Maximum:=b1+sf;
        end;
      end;
    end;  }
end;

procedure Addpoints(Chart:TChart);
//var sf:extended;
begin
     with Chart do
      begin
      with BottomAxis do
       begin
       // Logarithmic:=true;
       // Automatic:=false;
        Automatic:=True;
       { sf:=0.05*(Series[0].MaxXValue-Series[0].MinXValue);
        if (Series[0].MinXValue-sf)>0 then Minimum:=Series[0].MinXValue-sf
            else Minimum:=0.0005;

       Maximum:=Series[0].MaxXValue+sf;      }
      end;
     with LeftAxis do
      begin
    //   Logarithmic:=true;
      // Automatic:=false;
        Automatic:=True;
      { sf:=0.05*(Series[0].MaxYValue-Series[0].MinYvalue);
       if (Series[0].MinYValue-sf)>Maximum then
        begin
         Maximum:=Series[0].MaxYValue+sf;
         if (Series[0].MinYValue-sf)>0 then
                                Minimum:=Series[0].MinYValue-sf
                                        else
                                Minimum:=Series[0].MinYValue*0.95;
        end
        else begin
         if (Series[0].MinYValue-sf)>0 then
                                Minimum:=Series[0].MinYValue-sf
                                        else
                                Minimum:=Series[0].MinYValue*0.95;
         Maximum:=Series[0].MaxYValue+sf;
        end;                                  }
      end;
     end;
end;


// Эта функция также нужна разве что мне для справки, т.к. вывод на график будет
// в основной программе.
// Или нет, не совсем понятно - возможно программа использует построенные графики
// для дальнейших расчетов, что тоже не совсем правильно.
procedure AddExpPoints(Chart5:TChart;Chart6:TChart);
var i:word;
begin
 with Chart5.Series[1] as TPointSeries do
  begin
   Clear;
   Pointer.HorizSize:=2;
   Pointer.VertSize:=2;
   for i:=0 to NumberOfPoints do
    if MagField_spektr[i]=0 then AddXY(0.001,GxxExp[i],'',clTeeColor)
     else AddXY(MagField_spektr[i],GxxExp[i],'',clTeeColor);
   end;
  with Chart6.Series[1] as TPointSeries do
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
   if MagField_spektr[0]>0 then
    begin
     f:=ln(MagField_spektr[0])/ln(10);
     Xmin:=trunc(f);
    end
     else XMin:=-3;
    f:=ln(MagField_spektr[NumberOfPoints])/ln(10);
    Xmax:=trunc(f);
    if frac(f)>0.001 then Xmax:=Xmax+1;
end;




procedure MakeInterpolate(Chart5:TChart;Chart6:TChart);
var temp_l,temp_t,AGxx,AGxy,AField:Data_spektr;
    sf,lm,p,p1:extended;
    i,j:word;
    Lmin,Lmax,k:integer;
procedure CS(var X,F,C:Data_spektr;p1,pn:extended);
var i,j,m:word;
    K:Data_spektr;
    A,B,R:extended;
begin
   K[1]:=0;C[1]:=P1;
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
begin i:=1; while (X1>X[i]) do inc(i);
 j:=i-1;A:=F[j];B:=X[j];Q:=X[i]-B;
 R:=X1-B;P:=C[i];D:=C[i+1];
 B:=(F[i]-A)/q-(D+2*P)*Q/3;D:=(D-P)/3/Q;
 SP:=A+R*(B+R*(P+D*R));
end;
begin
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

   AddExpPoints(Chart5,Chart6);
   With Chart5 do
    begin
     Series[0].Clear;
    end;
   with Chart6 do
    begin
     Series[0].Clear;
    end;
   GetLnLimits(Lmin,Lmax);
   SizeData:=(Lmax-Lmin+1)*sizeof(ImageDat);
   InitArray;
   k:=0;
   for i:=0 to (lmax-lmin) do
    begin
     lm:=exp((lmin+i)*ln(10));
     sf:=lm;
     for j:=1 to PointPerInt-1 do
      begin
       IntMagField^[k]:=sf;
       IntGxx^[k]:=sp(MagField_spektr,GxxExp,temp_l,sf);
       IntGxy^[k]:=sp(MagField_spektr,GxyExp,temp_t,sf);
       Chart5.Series[0].AddXY(IntMagField^[k],IntGxx^[k],'',clTeeColor);
       Chart6.Series[0].AddXY(IntMagField^[k],IntGxy^[k],'',clTeeColor);
       sf:=lm*exp(j/PointPerInt*ln(10));
       if sf>MagField_spektr[NumberofPoints] then break;
       inc(k);
      end;
     if sf>MagField_spektr[NumberofPoints] then break;
    end;
   AddPoints(Chart5);
   AddPoints(Chart6);
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

procedure MakeMNK( a:boolean;Chart5:TChart;Chart6:TChart);
var tmp_m:mat;
    coef_t,coef_l:Data_spektr;
    kind,i,j,k,Lmin,Lmax:integer;
    lm,sf:extended;
 function pow(x:extended;y:word):extended;
 var p:extended;
     i:word;
  begin
  {if x>0 then
   p:=exp(y*ln(x))
         else  p:=1;}
   p:=1;
   for i:=1 to y do p:=p*x;
   pow:=p;
  end;
procedure BAS(N,M,L:word;X1:extended;var X,T:Data_spektr);
var
  k:word;z,r:extended;
begin
 z:=2*(x1-x[0])/(x[n]-x[0])-1;
 t[0]:=1;
 t[1]:=z;
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
    p:array [0..MaxPoints,0..5*MaxPoints] of extended;
begin
 for i:=0 to N do begin bas(n,m,l,x[i],x,t);
  for j:=0 to m do p[j,i]:=t[j]
 end;
 for k:=0 to m do begin
  for j:=k to m do begin s:=0; r:=0;
   for i:=0 to n do begin q:=p[k,i];s:=s+q*p[j,i];
    if j=m then r:=r+q*f[i];
   end; a[k,j]:=s;a[j,k]:=s;
  end; a[k,m+1]:=r;  // Warning you are not right.
 end;
end;

procedure gauss(N:word;var a:mat; var x:Data_spektr);
 var i,j,k,k1,n1:word; r,s:extended;
 begin
   n1:=N+1;
   for k:=0 to n do begin k1:=k+1;s:=a[k,k];
    for j:=k1 to n1 do a[k,j]:=a[k,j]/s;
      for i:=k1 to n do begin r:=a[i,k];
          for j:=k1 to n1 do a[i,j]:=a[i,j]-a[k,j]*r
      end;
    end;
     for i:=n downto 0 do begin s:=a[i,n1];
     for j:=i+1 to n do s:=s-a[i,j]*x[j];
     x[i]:=s;
    end;
  end;

procedure fi(n,m,l:word;var c,x:Data_spektr;var x1,s:extended);
var i:word;t:Data_spektr;
begin
 s:=c[0];
 bas(n,m,l,x1,x,t);
 for i:=1 to m do s:=s+c[i]*t[i]
end;
begin
 {if not(test) then
  begin}
   Power_spektr:=3;  // какая-то степень, или мощность, сильно похоже на кол-во типов носителей
   Kind:=2;  // разновидности, тоже пока не ясно зачем и что
   // эти вызовы заполняют матрицу tmp_m
   // которая кстати уже не пустая О_о
   gram(NumberOfPoints,Power_spektr,Kind,MagField_spektr,GxxExp,tmp_m);
   // Гаусс - судя по всему решает матрицу методом Гаусса, сохраняет всё в coef-l
   gauss(Power_spektr,tmp_m,coef_l);
   gram(NumberOfPoints,Power_spektr,Kind,MagField_spektr,GxyExp,tmp_m);
   gauss(Power_spektr,tmp_m,coef_t);
   With Chart5 do
    begin
     Series[0].Clear; // чистим графики компонент
    end;
   with Chart6 do
    begin
     Series[0].Clear;
    end;
   if a then
     begin
      AddExpPoints(Chart5,Chart6); // добавляет точки на график, не логарифмически
      GetLnLimits(Lmin,Lmax); // получает пределы, логарифмические
      SizeData:=(Lmax-Lmin+1)*sizeof(ImageDat); // считаем размер данных
      InitArray;  // выделяем его
      k:=0;
      for i:=0 to (lmax-lmin) do
       begin
       lm:=exp((lmin+i)*ln(10));
       sf:=lm;
       for j:=1 to PointPerInt-1 do
        begin
        IntMagField^[k]:=sf;
{        IntGxxExp^[k]:=fi(Power,coef_l,sf);
        IntGxyExp^[k]:=fi(Power,coef_t,sf);}
        // а тут происходит самое страшное - считаются два основных графика
        fi(NumberOfPoints,Power_spektr,Kind,coef_l,MagField_spektr,sf,IntGxx^[k]);
        fi(NumberOfPoints,Power_spektr,Kind,coef_t,MagField_spektr,sf,IntGxy^[k]);
        // и судя по всему ось х - действительно подвижность
        // а ось у - это величины компонент тензора проводимости
        // но они как-то модифицированы
        // видимо в предыдущей функции, надо уточнить
        Chart5.Series[0].AddXY(IntMagField^[k],IntGxx^[k],'',clTeeColor);
        Chart6.Series[0].AddXY(IntMagField^[k],IntGxy^[k],'',clTeeColor);
        sf:=lm*exp(j/PointPerInt*ln(10));
        if sf>MagField_spektr[NumberofPoints] then break;
        inc(k);
       end;
      if sf>MagField_spektr[NumberofPoints] then break;
     end;
     AddPoints(Chart5);
     AddPoints(Chart6);
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
     if j1=2 then begin
       if ii<>k1 then s:=s+sqr(B_spektr[ii]) end else
        if ii<>k1 then s:=s+sqr(B_spektr[ii])*GetElem(j1-1,k1,ii+1);
    GetElem:=s;
 end;


procedure MakeMatrC;
var
 j,k{,i1,i2,i3,i4,i5,i6,i7,i8}:word;
begin
 fillchar(Cr, Sizeof(Cr), 0);
 fillchar(Cl, Sizeof(Cl), 0);
 fillchar(Cr_t,Sizeof(Cr_t),0);
 fillchar(Cl_t,Sizeof(Cl_t),0);
 fillchar(Am, Sizeof(Am), 0);
 for j:=1 to NumberOfPoints do
  for k:=1 to NumberOfPoints do begin
   if j=1 then Cr[j,k]:=1 else
     Cr[j,k]:=GetElem(j,k,1);
     Cl[j,k]:=-Cr[j,k]*B_spektr[k];
    end;
end;

procedure MakeMatrA;
 var i,j,k:word;
 begin
  fillchar(Am,sizeof(Am),0);
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
procedure MobilitySpectrumFunc(chtSpectr:TChart);
var
  Sf,lm:extended;
  j,k,i:word;
  Lmin,Lmax:integer;
{  LogAxis,NoLogAxis:boolean;}


begin
   for i:=0 to NumberOfPoints do
    begin
     B_spektr[i+1]:=MagField_spektr[i];
     Gxx_sp[i+1]:=GxxExp[i];
     Gxy_sp[i+1]:=GxyExp[i];
    end;
   inc(NumberOfPoints);
   MakeMatrC;
   InverseMatrC(Cr,Cr_t,Sf,NumberOfPoints);
   if Sf=0 then begin MessageDLG('Определитель равен нулю!',mtError,[mbOK],0);
         dec(NumberOfPoints);
         exit; end;
   if MagField_spektr[0]=0 then begin
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
    MakeMatrA;
    // интересно что это за числа ниже?
    Tred2(NumberOfPoints,6e-4913,Lv,Xv,am,Qm,bulua);
    Imtql2(NumberOfPoints,5.42e-20,Lv,Xv,Qm,bulua);
   with chtSpectr do
    begin
     Series[0].Clear;
     Series[1].Clear;
    end;
   Lmin:=MSLeft;LMax:=MSRight;
   SizeData:=(Lmax-Lmin+1)*sizeof(ImageDat);
   InitArray2;
   k:=0;
   for i:=0 to (lmax-lmin) do
    begin
     lm:=exp((lmin+i)*ln(10));
     sf:=lm;
     for j:=1 to PointPerInt-1 do
      begin
       Mobility^[k]:=sf;
       Spectr_e^[k]:=S_s(-sf);
       Spectr_p^[k]:=S_s(sf);
       chtSpectr.Series[0].AddXY(Mobility^[k],Spectr_e^[k],'',clTeeColor);
       chtSpectr.Series[1].AddXY(Mobility^[k],Spectr_p^[k],'',clTeeColor);
       sf:=lm*exp(j/PointPerInt*ln(10));
       inc(k);
       if sf>10 then break;
      end;
     if sf>10 then break;
    end;
    GridPoints:=k-1;
   AddPoints2;
   dec(NumberOfPoints);
end;


function GxxExpi:extended;
var i:word;
begin
  Result:=0;
  for i:=0 to gridPoints do
   result:=result+Axx^[i]*(QSpectr_p^[i]+QSpectr_e^[i]);
end;

      // Загрузка данных из файла
procedure btnLoadTenzorClick(OpenDialog2:TOpenDialog;Series5:TPointSeries;
LineSeries1:TLineSeries);

//MagField_spektr[i],GxxExp[i],GxyExp[i]
// Это поле, продольная и поперечная компоненты тензора проводимости.
// Вот в эти три массива нужно записать входные данные.


begin
  {
   MakeMNK(true);
   MobilitySpectrumFunc;  // спектр подвижности. Начало.
 end;     }
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
var i:Integer;
   con, mob:Extended;
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


procedure btnSpectrResultClick(Sender: TObject);
begin
  {
  if ODSpektrRes.Execute then
  begin
     AssignFile(Data_File, ODSpektrRes.FileName);
     ReplacementOfSeparator;
     LoadSpektrResults;
  end; }
end;

// Запись результатов Спектра подвижности
procedure WriteSpektrResults;
var i:Integer;
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

procedure btnSpectrResult1Click(Sender: TObject);
begin
 { if SDSpektrRes.Execute then
  begin
    AssignFile(Data_File, SDSpektrRes.FileName);
    if not Interval then  // если таблица не пуста и интервалы рассчитаны, то
    WriteSpektrResults;
  end;}
end;

procedure StringGrid3SelectCell(Sender: TObject; ACol,
  ARow: Integer; var CanSelect: Boolean;StringGrid3:TStringGrid);
begin
 //StringGrid3.RowInFocus:=ARow;
end;    

////////////////////////////////////////////////////////////////////////////
/////////////////////// КОНЕЦ "ХОЛЛ. ПОДВИЖНОСТЬ"///////////////////////////
////////////////////////////////////////////////////////////////////////////

// входные параметры есть, надо подумать что мы возвращаем, где оно хранится
// и как мы это будем возвращать.
procedure RunMobilitySpectrum (MagneticFieldP,Exx,Exy: PData_spektr);
begin
  ;
end;  


{
Как я понял, для построения спектра подвижности нужно:

Вызвать MakeMNK(true)
и  MobilitySpectrumFunc;

}

exports

RunMobilitySpectrum;

begin
     // Картина такая:
     // Надо выделить память.
     // Получить входные данные и сохранить их в нужные массивы.
     // Вызвать две заветные функции.
     // И отдать результаты.
end.



{

unit UnitMain;



const
    EKT=150.574; {77} // используется в ФМЭ
{    EKT=136.4; {85}
{     EKT=115.9; {100}
{    EKT=92.75; {125}
   {  Mobility                         [ m**2/V*s ]     ;
      Magnetic Field                   [ Tesla ]        ;
      Electric Field                   [ V/m ]          ;
      Surface Recombination Velocity   [ m/s ]          ;
      All other values in SI                             ;   }
 { MaxPoints=100;
  PointPerInt=50;
  MaxParameters=8;    // Максимальное число параметров
  MaxRepeat=100;      // Максимальное число повторов
type
  EditArray = array[1..36] of TEdit;
  SeriesArray = array[1..9] of TBarSeries;  
  Holl   = array[1..4, 1..21] of Extended;
  HollB  = array[1..2, 1..11] of Extended;
  HollRh = array[1..11] of Extended;
  vector = array [0..100] of Extended;
  DataValue = array[1..MaxParameters, 1..MaxRepeat] of Extended;
    procedure CreateTabs;
    procedure FormCreate(Sender: TObject);



    procedure Button2Click(Sender: TObject);
    procedure Button4Click(Sender: TObject);
    procedure N2Click(Sender: TObject);
    procedure N100Click(Sender: TObject);
    procedure N8Click(Sender: TObject);
    procedure Gistogram(mass:DataValue; Ser:SeriesArray; m,n:Integer);
    procedure Statistic(mass:DataValue; Edit:EditArray; m,n:Integer);
    procedure N9Click(Sender: TObject);

    procedure btnLoadTenzorClick(Sender: TObject);
    procedure MakeInterpolate;
    procedure MakeMNK(a:boolean);
    procedure MakeLagranj;
    procedure InitArray;
    procedure InitArray2;
    procedure AddExpPoints;
    procedure GetLnLimits(var Xmin,Xmax:integer);
    procedure Addpoints(Chart:TChart);
    procedure AddPoints2;
    procedure MobilitySpectrumFunc;
    procedure chtSpectrClickSeries(Sender: TCustomChart; Series: TChartSeries;
      ValueIndex: Integer; Button: TMouseButton; Shift: TShiftState; X,
      Y: Integer);
    procedure FormDestroy(Sender: TObject);
    procedure chtSpectrMouseMove(Sender: TObject; Shift: TShiftState; X,
      Y: Integer);
    procedure SetGridFocus(SGrid: TStringGrid; r:integer);
    procedure btnClearSpectrTableClick(Sender: TObject);
   
    procedure Button10Click(Sender: TObject);
    procedure Button12Click(Sender: TObject);

    procedure FileOpen_hall8;
    procedure ReplacementOfSeparator;


    procedure LoadConfig;

    procedure LoadSpektrResults;
    procedure LoadFaradData;
    procedure LoadFaradRes;
    procedure LoadHall8Res;
    procedure LoadFMEData;
    procedure LoadFMERes;


    procedure WriteSpektrResults;
    procedure WriteFaradRes;
    procedure WriteHall8Res;
    procedure WriteFMERes;
    function  Interval:Boolean;
    procedure Normir;
    procedure Scale_FME(FpFme:Holl);
    procedure ShowGraphics(Edit:EditArray; n:Integer; Graph:Proc);

    procedure btnSpectrResult1Click(Sender: TObject);
    procedure btnSpectrResultClick(Sender: TObject);
    procedure btnFeatOnceClick(Sender: TObject);
    procedure btnFeatMultiClick(Sender: TObject);
    procedure N32Click(Sender: TObject);

    procedure btnSaveFeatResultsClick(Sender: TObject);
    procedure btnLoadFeatResultsClick(Sender: TObject);
    procedure Button28Click(Sender: TObject);
    procedure N35Click(Sender: TObject);

    procedure StringGrid3SelectCell(Sender: TObject; ACol, ARow: Integer;
      var CanSelect: Boolean);
    procedure Button33Click(Sender: TObject);
    procedure LoadModelData;
    procedure NormModel;
    p
    procedure btnFeatClick(Sender: TObject);


  private
   FNextViewer:HWnd;
    procedure WMChangeCBChain(var Msg: TWMChangeCBChain); message WM_CHANGECBCHAIN;
    procedure WMDrawClipboard(var Msg: TWMDrawClipboard); message WM_DRAWCLIPBOARD;
  public
    { Public declarations }
 { end;
  Splain=record
   A,B,C,D:extended;
  end;
  PeakInfo=record
    Index:word;
    Value:extended;
    end;
  ImageDat=array [0..(PointPerInt-1)*PointPerInt] of extended;
  PImageDat=^ImageDat;  // указатель на массив из двух элементов.
  //ImageDat3=array[0..1] of peakinfo;
  //PImageDat3=^ImageDat3;
  Data_spektr=array [0..MaxPoints] of extended;
  mat=array[0..MaxPoints,0..MaxPoints] of extended;
  Dat1=array[1..MaxPoints] of extended;
  Dat2=array[1..MaxPoints,1..MaxPoints] of extended;
  Dat3=array[1..MaxPoints,1..2*MaxPoints] of extended;

var
  Form1: TForm1;
  e:Extended;
 
  NPoint,NPoint_hm                       :Integer;
  fil, fil_hall8, Data_File, Config_File :text;
  muP,Ex,D,UPh,UGrad,UGrad0,UDif,UDif0,
  Alfa,Alfa0,K_trap,w1,w2,g_s,muE1       :Extended;
  MagField,FpeExp,Fpe,FmExp,Fm,RoExp,Ro,
  FpExp,Fp, MagField_FME                 :vector;
  GraphON_FME, GraphON_hall,
  GraphOn_Farad                          :Boolean;
  d1                                     :DataValue;
  DefaultDir, Recent_File_FME,
  Recent_File_Farad                      :string;
  RowInFocus                             :Word;
 /////////////// спектр подвижности///////////
  NumberofPoints,Power_spektr,GridPoints :word;
  MagField_spektr,Gxx,Gxy,GxxExp,GxyExp  :Data_spektr;
  IntGxx,IntGxy,IntMagField,Spectr_e,Spectr_p,Mobility,
  QSpectr_e,QSpectr_p,Axx,Axy,Axx_d,Axy_d,
  dIntGxx,dIntGxy,QGxx,QGxy,dQGxx,dQGxy  :PImageDat;
  //Peak_e,Peak_p,Valley_e,Valley_p        :PImageDat3;
  SizeData                               :longint;
  B_spektr,Gxx_sp,Gxx_MC,Gxy_MC,
  Gxy_sp,Xr,Lv,Xv,Mv,Vpr                 :Dat1;
  Am,Qm,Cl,Cr,Cl_t,Cr_t,Cm,Cm_t          :Dat2;
  bulua                                  :boolean;
  MSRight,MSLeft                         :integer;
  Min_Spectr,Coef1,Coef2,Mu_min,Mu_max,
  W,F_s,A1,An,B1,Bn, Ves1, Ves2          :extended;
  // Геометрия Фарадея
  FpI                                    :Holl;
  LoadOld                                :Boolean;
  Mp, Post                               :Extended;
  ves:array[1..11] of Extended;
  
implementation

uses Unit2, Unit3, Unit4, Unit5, Unit6, Unit7;


{$R *.dfm}
  {
procedure TForm1.FormCreate(Sender: TObject);// создане формы
begin
  CreateTabs;  // заполнение таблиц
  DefaultDir:=GetCurrentDir;
  e:=1.602176487*Power(10,-19);
  DecimalSeparator:=',';
  FnextViewer:=SetClipboardViewer(Handle);
  GraphON_FME:=True;  // графики включены
  GraphON_hall:=True;
  GraphOn_Farad:=True;
  // из Спектра подвижности
  MSLeft:=-2;
  MSRight:=1;
  Coef1:=2.1;
  Coef2:=1.6;
  Mu_min:=0.01;
  Mu_max:=100;
  Min_Spectr:=1e-4;
  RowInFocus:=1;
  SetGridFocus(StringGrid3,RowInFocus);
  AssignFile(Config_File,DefaultDir+'\Data\Config.txt');
  AssignFile(Data_File,DefaultDir+'\Data\Recent_Data.dat');
  LoadConfig;
  ReplacementOfSeparator;
  
  Memo7.Lines.LoadFromFile(DefaultDir+'\Data\FME_Prototype.dat');
  Label46.Caption:='ПРОТОТИП файла .fm';
  Memo16.Lines.LoadFromFile(DefaultDir+'\Data\Farad_Prototype.dat');
  Label39.Caption:='ПРОТОТИП файла .mpc';
end;

procedure TForm1.CreateTabs;
var j:Integer;
begin


  // Создание таблицы
  StringGrid3.Cells[0,1]:='Тяжёлые дырки';
  StringGrid3.Cells[0,2]:='Лёгкие дырки';
  StringGrid3.Cells[0,3]:='Электроны';
  StringGrid3.Cells[1,0]:='Концентрация';
  StringGrid3.Cells[2,0]:='Подвижность';

end;
// возникает при изменении очереди наблюдателей за буфером обмена
procedure TForm1.WMChangeCBChain(var Msg: TWMChangeCBChain);
begin
 if Msg.Remove = FNextViewer then
  FnextViewer:=Msg.Next
 else
  SendMessage(FNextViewer, WM_CHANGECBCHAIN, Msg.Remove, Msg.Next);
end;

// возникает при изменении буфера обмена
procedure TForm1.WMDrawClipboard(var Msg: TWMDrawClipboard);
begin
  SendMessage(FNextViewer, WM_DRAWCLIPBOARD, 0, 0);
end;




////////////////////////////////////////////////////////////////////////////
///////////////// НАЧАЛО "ХОЛЛ. МНОГОЗОННАЯ ПОДГОНКА" //////////////////////
////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
////////////////////// ПРОЦЕДУРЫ ДЛЯ ЗАГРУЗКИ И СОХРАНЕНИЯ////////////////////
//////////////////////////////////////////////////////////////////////////////
       {
procedure TForm1.ReplacementOfSeparator;   // Убирает запятые из файла
var j,position,size_m:Integer;
    MiniString, MaxiString:string;

begin
  MiniString:='';
  MaxiString:='';
{$I-}
 {   Reset(Data_File);
{$I+}
 { if IOResult<>0 then
   begin
    ShowMessage('Ошибка открытия файла');
    Exit;
   end
  else
    begin
     while not Eof(Data_File) do // Убираем запятые, если они есть,
                            // и заменяем их на '.'
      begin
        Readln(Data_File,  ministring);
        size_m:=SizeOf(MiniString);
        for j:=0 to size_m+1 do
        begin
          position:=Pos(',',MiniString);
          if position<>0 then
          MiniString[position]:='.';
        end;

        MaxiString:=MaxiString+MiniString+#13#10;
      end;
      CloseFile(Data_File);
    end;

  Rewrite(Data_File);    // Переписываем файл без запятых
  Write(Data_File, MaxiString);
  CloseFile(Data_File);
end;

procedure TForm1.LoadConfig;
begin
{$I-}
    Reset(Config_File);
{$I+}
 { if IOResult<>0 then    // обходим ошибку открытия
   begin
    ShowMessage('Ошибка открытия файла конфигураций');
    Exit;
   end
  else
 begin
  try    // обходим ошибку загрузки
   Readln(Config_File, NPoint_hm);
   Readln(Config_File, NPoint);
  except
    ShowMessage('Ошибка загрузки конфигураций');
    Rewrite(Config_File);
    Writeln(Config_File,21);
    Writeln(Config_File,11);
  end;
  CloseFile(Config_File);
 end;
end;

        }

