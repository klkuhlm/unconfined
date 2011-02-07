#include "../stdafx.h"
#include "../nrecip/nr.h"

using namespace std;

const double Err=1.0e-8;//Relative error for Laplace transform inversion
const double pi = 4.0*atan(1.0);//defining the number "pi"

struct Aquif_params{//Stucture for input parameters
  char pumpwelltype,obswelltype,deriv,inversion;
  int numtimes;
  double Qrate, obs_rw, shapefactor, piez_pos_r, piez_pos_z, screen_top, screen_bottom, anisotropy, aquif_Kz, aquif_Kr, aquif_Ss, aquif_Sy, ini_satu_thickness, decay_const;
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

complex<long double> func_u(long double a, complex<long double> p, Aquif_params params)
{//Theis solution in Laplace & Hankel space incorporating wellbore storage

  double tDb,Cw,b;
  
  b = params.ini_satu_thickness;//Initial saturated thickness
  Cw = pi*pow(params.obs_rw,2.0);//Wellbore storage
  tDb = Cw/(pow(b,2.0)*params.shapefactor*params.aquif_Ss);//Dimensionless wellbore response time
  
  return ((long double)2.0/(p*(p + pow(a,(long double)2.0))*(p*(long double)tDb + (long double)1.0)));
}

complex<long double> delta(long double a, complex<long double> p, Aquif_params params)
{//The function \Delta that appears in solution
  long double kappa,aD,a_r,S_y,betaD,b;
  
  b = (long double)params.ini_satu_thickness;//Initial saturated thickness
  kappa = (long double)(params.aquif_Kz/params.aquif_Kr);//Anisotropy ratio
  S_y = params.aquif_Sy;//Specific yield
  a_r = (long double)(params.aquif_Kr/params.aquif_Ss);//Hydraulic diffusivity
  aD = (long double)(params.aquif_Kz/(params.decay_const*b*S_y));//Measure of vertical diffusivity associated with delayed yield
  betaD = (long double)params.decay_const*pow(b,(long double)2.0)/a_r;//Dimensionless Moench decay constant

  complex<long double> eta = sqrt((p + pow(a,(long double)2.0))/kappa);
  return (aD*eta*sinh(eta) + (p/(p + betaD))*cosh(eta));
}

complex<long double> LH_drawdown(long double a, complex<long double> p, Aquif_params params)
{//Computes drawdown in a piezometer
  complex<long double> eta,vD;
  long double kappa,zD,betaD,a_r,b,Cw,tDb;
  
  b = (long double)params.ini_satu_thickness;//Initial saturated thickness
  zD = (long double)(params.piez_pos_z/b);//Dimensionless vertical position of piezometer
  kappa = (long double)(params.aquif_Kz/params.aquif_Kr);//Anisotropy ratio
  a_r = (long double)(params.aquif_Kr/params.aquif_Ss);//Hydraulic diffusivity
  betaD = (long double)params.decay_const*pow(b,(long double)2.0)/a_r;//Dimensionless Moench decay constant

  Cw = (long double)pi*pow(params.obs_rw,2.0);//Wellbore storage
  tDb = Cw/(pow(b,(long double)2.0)*params.shapefactor*params.aquif_Ss);//Dimensionless wellbore response time
  
  eta = sqrt((p+pow(a,(long double)2.0))/kappa);
  vD = (long double)2.0*cosh(eta*zD)/((p+pow(a,(long double)2.0))*(p+betaD));
  return (func_u(a,p,params) - vD/delta(a,p,params))/(p*(long double)tDb + (long double)1.0);
}

complex<long double> LH_drawdown_slope(long double a, complex<long double> p, Aquif_params params)
{//Computes drawdown derivative in a piezometer
  complex<long double> eta,vD;
  long double kappa,zD,betaD,a_r,b,Cw,tDb;
  
  b = (long double)params.ini_satu_thickness;//Initial saturated thickness
  zD = (long double)(params.piez_pos_z/b);//Dimensionless vertical position of piezometer
  kappa = (long double)(params.aquif_Kz/params.aquif_Kr);//Anisotropy ratio
  a_r = (long double)(params.aquif_Kr/params.aquif_Ss);//Hydraulic diffusivity
  betaD = (long double)params.decay_const*pow(b,(long double)2.0)/a_r;//Dimensionless Moench decay constant

  Cw = (long double)pi*pow(params.obs_rw,2.0);//Wellbore storage
  tDb = Cw/(pow(b,(long double)2.0)*params.shapefactor*params.aquif_Ss);//Dimensionless wellbore response time
  
  eta = sqrt((p+pow(a,(long double)2.0))/kappa);
  vD = (long double)2.0*cosh(eta*zD)/((p+pow(a,(long double)2.0))*(p+betaD));
  return (p/(p*(long double)tDb + (long double)1.0))*(func_u(a,p,params) - vD/delta(a,p,params));
}

complex<long double> LH_mean_drawdown(long double a, complex<long double> p, Aquif_params params)
{//Drawdown in Fully penetrating observation well
  complex<long double> eta,vD;
  long double kappa,a_r,betaD,b,Cw,tDb;

  b = (long double)params.ini_satu_thickness;//Initial saturated thickness
  kappa = (long double)(params.aquif_Kz/params.aquif_Kr);//Anisotropy ratio
  a_r = (long double)(params.aquif_Kr/params.aquif_Ss);//Hydraulic diffusivity
  betaD = params.decay_const*pow(b,(long double)2.0)/a_r;//Dimensionless Moench decay constant
  
  Cw = (long double)pi*pow(params.obs_rw,2.0);//Wellbore storage
  tDb = Cw/(pow(b,(long double)2.0)*params.shapefactor*params.aquif_Ss);//Dimensionless wellbore response time
  
  eta = sqrt((p+pow(a,(long double)2.0))/kappa);
  vD = (long double)2.0*sinh(eta)/((p+pow(a,(long double)2.0))*(p+betaD));
  return (func_u(a,p,params) - vD/(eta*delta(a,p,params)))/(p*(long double)tDb + (long double)1.0);
}

double fixedTalbot(int M, double t, double rho, complex<double> *Lft)
{//Implementation of Fixed Talbot method
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
{//Inversion of Laplace transform using Fixed Talbot method
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
    if(params.obswelltype == 'P')
      v[i] = LH_drawdown(a,p,params);
    else if(params.obswelltype == 'M')
      v[i] = LH_mean_drawdown(a,p,params);
    else
      v[i] = LH_drawdown_slope(a,p,params);
  }
  return fixedTalbot(M,tD,rho,v);
}

DP func(DP x, double t, Aquif_params params)
{//This is the function that is integrated during Hankel transform inversion

  double rD = params.piez_pos_r/params.ini_satu_thickness;
  return x*lap_invert_u(t,x,params)*abs(NR::bessj0(x*rD));
}

DP integral_Ik(DP func(DP, double, Aquif_params), int k, double tD, double *j0, Aquif_params params,int M, double *x, double *w)
{//Returns the integral I_k
  double a,b,s;
  double rD = params.piez_pos_r/params.ini_satu_thickness;
  a = j0[k]/rD;
  b = j0[k+1]/rD;

    return NR::qgausBM(func,a,b,M,x,w,tD,params); 
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
  ifstream inp("parameters.dat");
  if(!inp){
    cerr<<"File could not be opened"<<endl;
    exit(1);}
  const int MAX=80;
  char ch,buff[MAX];
  inp.getline(buff,MAX);
  inp>>params.pumpwelltype;//choose between partial (P) and full (F) penetration for pumping well
  inp>>params.obswelltype;//Choose between full (F) and partial (P) penetration for observation well
  inp>>params.deriv;//choose between drawdown (N) or its time derivative (Y)
  inp>>params.inversion;
  inp>>params.numtimes;
  inp>>params.Qrate;
  inp>>params.obs_rw;
  inp>>params.shapefactor;
  inp>>params.piez_pos_r;
  inp>>params.piez_pos_z;
  inp>>params.screen_top;//top of pumping well screen
  inp>>params.screen_bottom;//bottom of pumping well screen
  inp>>params.aquif_Kr;
  inp>>params.anisotropy;
  inp>>params.aquif_Ss;
  inp>>params.aquif_Sy;
  inp>>params.ini_satu_thickness;
  inp>>params.decay_const;
  params.aquif_Kz = params.anisotropy*params.aquif_Kr;
}

int main()
{
  int N=0,M=10,gl_M=100;
  double *tD,tD2,rD,zD,j0[M+1],max_tD,dtD,x[gl_M],w[gl_M],b,a_r,H,Q;
  DP sD;

  ifstream inp("besJ0zeros.dat");
  ifstream inp2("gausleg.dat");
  ifstream inp3("times.dat");
  ofstream outf("moench.dat",ios::out);
  outf<<setiosflags(ios::fixed)
    <<setiosflags(ios::showpoint)
    <<setprecision(6);
	
  for(int i=0;i<M+1;i++){//Read zeroes of zeroth order Bessel function
    if(i==0)
      j0[0]=0.0;
    else
      inp>>j0[i];}

  for(int i=0;i<gl_M;i++){//Read in the abscissas and weights of Gaussian quadrature
    inp2>>x[i]>>w[i];}

  Aquif_params params;
  aquif_parameters(params);
  
  b = params.ini_satu_thickness;
  a_r = params.aquif_Kr/params.aquif_Ss;
  Q = params.Qrate;

  H = 6.309e-5*Q/(4.0*pi*b*params.aquif_Kr);
  cout<<params.aquif_Kr<<"\t"<<params.anisotropy<<"\t"<<params.aquif_Ss<<"\t"<<params.aquif_Sy<<"\t"<<params.decay_const<<"\t"<<params.shapefactor<<endl;

  tD = new double[params.numtimes];
  for(int n=0;n<params.numtimes;n++){
    inp3>>tD[n];
    tD[n] *= 60.0*a_r/pow(b,2.0); 
  }

  if(params.inversion == 'Y'){//Run following block for parameter estimation
    for(int k=0;k<params.numtimes;k++){
      sD = totalIntgrl(func,N,M,tD[k],j0,params,gl_M,x,w);
      outf<<setw(10)<<H*abs(sD)/0.3048<<endl;
    }
  }else{//run following block for forward simulation
    tD2 = 1.0e-3;
    dtD = 1.0e-4;
    max_tD = 10.0*tD2;
    for(int k=0;k<10;k++){
      cout<<"log-cycle: "<< k+1 << endl;
      while(tD2 < max_tD){
        sD = totalIntgrl(func,N,M,tD2,j0,params,gl_M,x,w);
        outf<<tD2<<"\t"<<setw(9)<<abs(sD)<<endl;
        tD2 += dtD;}
      dtD *= 10.0;
      max_tD *= 10.0;
    }
  }
  
  return 0;
}
