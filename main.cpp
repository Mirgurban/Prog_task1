#include <iostream>
#include <math.h>

using namespace std;

#define pi 3.141592653589793238462643383279

double func(int k, int m, double h)
{
  double x = k * h;
  double y = m * h;
  return 2*(x-x*x + y - y*y);
}

double analit(int k, int m, double h)
{
  double x = k * h;
  double y = m * h;
  return x*(1-x)*y*(1-y);
}

double norm(double *u, int n, double h)
{
  double sum = 0;
  double max = 0;
  for (int k = 0; k < n; ++k) {
    for (int m = 0; m < n; ++m) {
      sum += (analit(k,m,h) - u[k * n + m]) * (analit(k,m,h) - u[k * n + m]);
      if ((analit(k,m,h) - u[k * n + m]) > max) {
        max = (analit(k,m,h) - u[k * n + m]);
      }
    }  
  }
//   return sqrt(sum)/n/n;
  return max;
}

void solveMatrix (int n, double *a, double *c, double *b, double *f, double *x)
{
	double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
    {
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
  
  for (int i = 0; i < n; ++i) {
    c[i] = c[0];
  }
}

void difference_scheme(double *a, double *b, double *c, int n, double tau, double h)
{
    for (int i = 0; i < n; ++i) {
        a[i] = tau/2/h/h;
        c[i] = -1 - tau/h/h;
        b[i] = tau/2/h/h;
    }
    a[0] = 0;
    b[n-1] = 0;
}

double TAU(double h)
{
  return h/2/sin(pi*h);
  // return 2*sin(pi*h)/h/h;
  // return 1;
}


double LU(double x, double y, double z, double h)
{
  // -(uk,m−1 − 2uk,m + uk,m+1)/h^2
  //       or
  //−(uk−1,m − 2uk,m + uk+1,m)/h^2
  return -(x-2*y+z)/h/h;
}

void F1(double *f, int m, int n, double tau, double h, double *u)
{
    for  (int i = 0; i < n - 2; ++i){
      f[i] = tau/2*(LU(u[m-1+n*(i + 1)],u[m + n*(i + 1)],u[m+1 + n*(i + 1)],h) 
                        - func(i + 1,m,h)) - u[m + n*(i + 1)];
    }
}

void F2(double *f, int k, int n, double tau, double h, double *u)
{
    for  (int i = 0; i < n - 2; ++i){
      f[i] = tau/2*(LU(u[i + 1 + n*(k - 1)],u[i + 1 + n*(k)],u[i + 1 + n*(k + 1)],h)
                         - func(k,i + 1,h)) - u[i + 1 + n*(k)];
    }
}

int main()
{
  int N = 100;
  double h = 1./(N-1);

  double *a = new double[N];
  double *b = new double[N];
  double *c = new double[N];
  double *f = new double[N];    
  double *u = new double[N*N + 1];
  double *temp = new double[N*N + 1];
  double tau = TAU(h);
  difference_scheme(a,b,c,N-2,tau,h);

  for (int i = 0; i < N * N + 1; ++i) {
    u[i] = 0;
    temp[i] = 0;
  }

  for (int l = 0; l < 12; ++l) {
    for (int m = 1; m < N - 1; ++m) {
        F1(f,m,N,tau,h,u);
        solveMatrix(N-2,a,c,b,f,temp + m*N + 1);
      }

      for (int k = 1; k < N - 1; ++k) {
        F2(f,k,N,tau,h,temp);
        solveMatrix(N-2,a,c,b,f,u + k*N + 1);
      }
    }

  cout << norm(u,N,h) << endl;
  // cout << u[N+1] << ' ' << analit(1,1,h) << ' ' << u[2*N-2] << ' ' <<  analit(1,N-2,h) << endl;
  if (N < 15){
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < N; ++j) {
        cout << u[i*N + j] << ' ';
      }
      cout << endl;
    }
  }
  delete []a;
  delete []b;
  delete []c;
  delete []temp;
  delete []u;
  delete []f;
}
