#include "../stdafx.h"
#include "../nrecip/nr.h"

using namespace std;

const double Err=1.0e-8;//Relative error for Laplace transform inversion
const double pi = 4.0*atan(1.0);//defining the number "pi"

struct Aquif_params{
  char welltype;
  double piez_pos_r, piez_pos_z, aquif_Kz, aquif_Kr, aquif_Ss, aquif_Sy, ini_satu_thickness;
};

DP NR::bessj0(const DP x)
{
	DP ax,z,xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934945152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

DP NR::midpntBM2(DP func(const DP, double, Aquif_params), const DP a, const DP b, const int n, double t, Aquif_params params)
{
	int it,j;
	DP x,tnm,del,ddel,sum;
	static DP s;

	if (n == 1) {
		return (s=(b-a)*func(0.5*(a+b),t,params));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=0;j<it;j++) {
			sum += func(x,t,params);
			x += ddel;
			sum += func(x,t,params);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

DP NR::gammln(const DP xx)
{
	int j;
	DP x,y,tmp,ser;
	static const DP cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
		-0.5395239384953e-5};

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<6;j++) ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

DP NR::factln(const int n)
{
	static DP a[101];

	if (n < 0) nrerror("Negative factorial in routine factln");
	if (n <= 1) return 0.0;
	if (n <= 100)
		return (a[n] != 0.0 ? a[n] : (a[n]=gammln(n+1.0)));
	else return gammln(n+1.0);
}

DP NR::factrl(const int n)
{
	static int ntop=4;
	static DP a[33]={1.0,1.0,2.0,6.0,24.0};
	int j;

	if (n < 0) nrerror("Negative factorial in routine factrl");
	if (n > 32) return exp(gammln(n+1.0));
	while (ntop<n) {
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}

DP NR::bico(const int n, const int k)
{
	return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

DP NR::qromoBM(DP func(const DP, double, Aquif_params), const DP a, const DP b, DP choose(DP (*)(const DP, double, Aquif_params), const DP, const DP, const int, double, Aquif_params), double t, Aquif_params params)
{
	const int JMAX=14, JMAXP=JMAX+1, K=5;
	const DP EPS=3.0e-9;
	int i,j;
	DP ss,dss;
	Vec_DP h(JMAXP),s(JMAX),h_t(K),s_t(K);

	h[0]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j-1]=choose(func,a,b,j,t,params);
		if (j >= K) {
			for (i=0;i<K;i++) {
				h_t[i]=h[j-K+i];
				s_t[i]=s[j-K+i];
			}
			polint(h_t,s_t,0.0,ss,dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j]=h[j-1]/9.0;
	}
	nrerror("Too many steps in routine qromo");
	return 0.0;
}

void NR::gauleg(const DP x1, const DP x2, Vec_O_DP &x, Vec_O_DP &w)
{
	const DP EPS=1.0e-14;
	int m,j,i;
	DP z1,z,xm,xl,pp,p3,p2,p1;

	int n=x.size();
	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=0;i<m;i++) {
		z=cos(3.141592654*(i+0.75)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=0;j<n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j+1.0)*z*p2-j*p3)/(j+1);
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z;
		x[n-1-i]=xm+xl*z;
		w[i]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i]=w[i];
	}

//for(int i=0;i<5;i++)
  //cout<<x[i]<<endl;

}

DP NR::qgausBM(DP func(const DP, double, Aquif_params), const DP a, const DP b,int M, double *x, double *w, double t, Aquif_params params)
{
	//static const DP x[]={0.1488743389816312,0.4333953941292472,
	//	0.6794095682990244,0.8650633666889845,0.9739065285171717};
	//static const DP w[]={0.2955242247147529,0.2692667193099963,
	//	0.2190863625159821,0.1494513491505806,0.0666713443086881};
	int j;
	DP xr,xm,dx,s;

	xm=0.5*(b+a);
	xr=0.5*(b-a);
	s=0;
	for (j=0;j<M/2;j++) {
		dx=xr*x[j];
		s += w[j]*(func(xm+dx,t,params)+func(xm-dx,t,params));
	}
	return s *= xr;
}

void NR::polint(Vec_I_DP &xa, Vec_I_DP &ya, const DP x, DP &y, DP &dy)
{
	int i,m,ns=0;
	DP den,dif,dift,ho,hp,w;

	int n=xa.size();
	Vec_DP c(n),d(n);
	dif=fabs(x-xa[0]);
	for (i=0;i<n;i++) {
		if ((dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=0;i<n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ((den=ho-hp) == 0.0) nrerror("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
	}
}

complex<double> func_u(double a, complex<double> p, Aquif_params params)
{
  return (2.0/(p*(p + pow(a,2))));
}

complex<long double> func_u(long double a, complex<long double> p, Aquif_params params)
{
  return ((long double)2.0/(p*(p + pow(a,2))));
}

complex<long double> delta(long double a, complex<long double> p, Aquif_params params)
{
  long double KzD,aDy,a_r,a_y;
  KzD = (long double)(params.aquif_Kz/params.aquif_Kr);
  a_y = (long double)(params.ini_satu_thickness*params.aquif_Kz/params.aquif_Sy);
  a_r = (long double)(params.aquif_Kr/params.aquif_Ss);
  aDy = a_y/a_r;

  complex<long double> eta = sqrt((p+pow(a,(long double)2.0))/KzD);
  return (cosh(eta) + (eta*aDy/p)*sinh(eta));
}

complex<long double> LH_drawdown(long double a, complex<long double> p, Aquif_params params)
{//Linearized kinematic condition at watertable
  complex<long double> eta,vD;
  long double KzD,zD;

  zD = (long double)(params.piez_pos_z/params.ini_satu_thickness);
  KzD = (long double)(params.aquif_Kz/params.aquif_Kr);

  eta = sqrt((p+pow(a,(long double)2.0))/KzD);
  vD = cosh(eta*zD)/delta(a,p,params);
  return func_u(a,p,params)*((long double)1.0 - vD);
}

complex<long double> LH_mean_drawdown(long double a, complex<long double> p, Aquif_params params)
{//Linearized kinematic condition at watertable
  complex<long double> eta,vD;
  long double KzD,zD;

  zD = (long double)(params.piez_pos_z/params.ini_satu_thickness);
  KzD = (long double)(params.aquif_Kz/params.aquif_Kr);

  eta = sqrt((p+pow(a,(long double)2.0))/KzD);
  vD = sinh(eta)/(eta*delta(a,p,params));
  return func_u(a,p,params)*((long double)1.0 - vD);
}

double fixedTalbot(int M, double t, double rho, complex<double> *Lft)
{
	double theta,sig,cot,sum = 0.0;
	complex<double> s;
	complex<double> I(0.0,1.0);

	for(int k=1;k<M;k++){
		theta = k*pi/M;
		cot = cos(theta)/sin(theta);
		s = rho*theta*(cot+I);
		sig = theta +(theta*cot-1.0)*cot;
		sum += real(exp(t*s)*Lft[k]*(1.0+I*sig));}

	return ((0.5*real(Lft[0])*exp(rho*t)+sum)*rho/M);
}

double lap_invert_u(double tD, double a, Aquif_params params)
{
  int M=8;
  complex<double> I(0.0,1.0);
  complex<double> p,v[M];
  double theta,rho,cot;

  rho = 2.0*M/(5.0*tD);
  for(int i=0;i<M;i++){
    theta = i*pi/M;
    cot = 1.0/tan(theta);
    if(i==0){
      p = complex<double>(rho,0.0);}
    else{
      p = rho*theta*(cot+I);}
    if(params.welltype == 'P')
      v[i] = LH_drawdown(a,p,params);
    else
      v[i] = LH_mean_drawdown(a,p,params);
  }
  return fixedTalbot(M,tD,rho,v);
}

DP func(DP x, double t, Aquif_params params)
{//Simply replace LH_v2 with LH_v of LH_meanv

  double rD = params.piez_pos_r/params.ini_satu_thickness;
  return x*lap_invert_u(t,x,params)*abs(NR::bessj0(x*rD));
}

DP integral_Ik(DP func(DP, double, Aquif_params), int k, double tD, double *j0, Aquif_params params,int M, double *x, double *w)
{//Returns the integral I_k
  double a,b,s;
  double rD = params.piez_pos_r/params.ini_satu_thickness;
  a = j0[k]/rD;
  b = j0[k+1]/rD;

  //if(tD<3000.0)
    return NR::qgausBM(func,a,b,M,x,w,tD,params); 
  //else
    //return NR::qromoBM(func,a,b,NR::midpntBM2,tD,params); 
}

DP Del_kIm(DP func(DP, double, Aquif_params),int k, double tD, double *j0, Aquif_params params,int M, double *x, double *w)
{//This function returns the difference term 
  double sum = 0.0;
  for(int m=0;m<(k+1);m++){
    sum += pow(-1.0,m)*NR::bico(k,m)*integral_Ik(func,k-m,tD,j0,params,M,x,w);}
  return sum;
}

DP totalIntgrl(DP func(DP, double, Aquif_params), int N, int M, double t, double *j0, Aquif_params params,int M2, double *x, double *w)
{//returns the total integral
  double sum = 0.0;
  for(int k=N;k<(M+1);k++){
    sum += pow(-1.0,k)*Del_kIm(func,k,t,j0,params,M2,x,w)/pow(2.0,k+1);}
  return sum;
}

void aquif_parameters(Aquif_params &params)
{
  double kappa;
  ifstream inp("neum_param.dat");
  if(!inp){
    cerr<<"File could not be opened"<<endl;
    exit(1);}
  const int MAX=80;
  char ch,buff[MAX];
  inp.getline(buff,MAX);
  inp>>params.welltype;
  inp>>params.piez_pos_r;
  inp>>params.piez_pos_z;
  inp>>params.aquif_Kr;
  inp>>kappa;
  inp>>params.aquif_Ss;
  inp>>params.aquif_Sy;
  inp>>params.ini_satu_thickness;
  params.aquif_Kz = kappa*params.aquif_Kr;
}

int main()
{
  int N=0,M=10,gl_M=100;
  double tD,rD,zD,j0[M+1],max_tD,dtD,x[gl_M],w[gl_M];
  DP sD;

  ifstream inp("besJ0zeros.dat");
  ifstream inp2("gausleg.dat");
  ofstream outf("neuman72.dat",ios::out);
  outf<<setiosflags(ios::fixed)
    <<setiosflags(ios::scientific)
    <<setiosflags(ios::showpoint)
    <<setprecision(6);
	
  for(int i=0;i<M+1;i++){//Read zeroes of zeroth order Bessel function
    inp>>j0[i];}

  for(int i=0;i<gl_M;i++){//Read in the abscissas and weights of Gaussian quadrature
    inp2>>x[i]>>w[i];}

  Aquif_params params;
  aquif_parameters(params);
  rD = params.piez_pos_r/params.ini_satu_thickness;
  tD = 1.0e-4;
  dtD = 1.0e-5;
  max_tD = 10.0*tD;
  for(int k=0;k<10;k++){
    cout<<"cycle: "<<k+1<<endl;
    while(tD < max_tD){
      sD = totalIntgrl(func,N,M,tD,j0,params,gl_M,x,w);
      outf<<tD/pow(rD,2.0)<<"\t"<<setw(9)<<abs(sD)<<endl;
      //cout<<tD<<"\t"<<abs(sD)<<endl;
      tD += dtD;}
    dtD *= 10.0;
    max_tD *= 10.0;
  }
  return 0;
}
