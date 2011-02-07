#include "../../stdafx.h"
#include "../../nrecip/nr.h"

using namespace std;

const double Err=1.0e-8;//Relative error for Laplace transform inversion
const double pi = 4.0*atan(1.0);//defining the number "pi"

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
}

int main()
{
  int N=0,M=100;
  Vec_O_DP x(M), w(M);
  DP x1 = 1.0;
  DP x2 = -1.0;
  ofstream outf("gausleg.dat",ios::out);
  outf<<setiosflags(ios::fixed)
	<<setiosflags(ios::scientific)
	<<setiosflags(ios::showpoint)
	<<setprecision(16);

  NR::gauleg(x1,x2,x,w);
  for(int i=M/2 -1;i>=0;i--){
    outf<<x[i]<<"\t"<<fabs(w[i])<<endl;}
  return 0;
}