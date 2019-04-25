# include <ctime>
# include <cmath>
# include <cctype>
# include <string>
# include <iostream>
# include <math.h>
# include <stdio.h>
# include <fstream>
# include <cassert>
# include <sstream>
# include <omp.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_sf_legendre.h>
# include <gsl/gsl_sf_expint.h>
# include <gsl/gsl_integration.h>
# include <gsl/gsl_errno.h>
# include <gsl/gsl_eigen.h>
# include <gsl/gsl_spline.h>
# include <gsl/gsl_sf_bessel.h>
# include <gsl/gsl_sf_gamma.h>
# include <gsl/gsl_sf_coupling.h>
# include <gsl/gsl_bspline.h>
# include <gsl/gsl_multifit.h>
# include <gsl/gsl_statistics.h>
# include <gsl/gsl_linalg.h>
# include <gsl/gsl_sf_erf.h>
# include <gsl/gsl_sf.h>
# include <gsl/gsl_sort_vector.h>
# include <gsl/gsl_rng.h>
# include <gsl/gsl_randist.h>
# include <gsl/gsl_roots.h>
# include <gsl/gsl_math.h>
# include <gsl/gsl_deriv.h>
# include <vector>
# include <numeric>
# include <algorithm>

using namespace std;

// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// This file contains some functions, based on GSL subroutines, used  in the measurements  **
// of the power spectrum                                                                   **
// Andres Balaguera Antolinez                                                              **
// ******************************************************************************************
// ******************************************************************************************
// ******************************************************************************************
// ********************************************************************
double long factorial(int n){
  if(n<=0 || n==1 ) return 1;
  else return tgammal(n+1);
}
// ********************************************************************

double real_sh(int l, int m, double theta, double phi){
  return  gsl_sf_legendre_sphPlm(l,m,cos(theta))*cos(m*phi);
}
// ********************************************************************

double imag_sh(int l, int m, double theta, double phi){
  return gsl_sf_legendre_sphPlm(l,m,cos(theta))*sin(m*phi);
}

// *******************************************************************************************************************************************************
double top_hat(double x){
  double ans;
  if(fabs(x)>0.5)ans=0.0;
  else if(fabs(x)==0.5)ans=0.5;
  else if(fabs(x)<0.5)ans=1.0;
  return ans;
}

// *******************************************************************************************************************************************************
void comp_time(time_t start, long full, long step){
  float fraction=100.0*((double)(step))/((double)full);
  cout.precision(3);
  time_t end;
  time(&end);
  double lapse=difftime(end,start);
  if(lapse<=60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse<<"  secs \r";cout.flush();
  if(lapse>60)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/60.<<"  mins \r";cout.flush();
  if(lapse>3600)cout<<"\r "<<fraction<<" % completed. Time elapsed: "<<lapse/3600.<<"  hrs \r";cout.flush();

}
// ********************************************************************
string to_string (double Number){
  stringstream ss;
  ss << Number;
  return ss.str();
}
// ********************************************************************
double gsl_integration(double (*function)(double, void *) ,void *p,double LowLimit,double UpLimit){
  gsl_integration_glfixed_table *wf =  gsl_integration_glfixed_table_alloc (500);
  gsl_function F;
  F.params   = p;  
  F.function = function;
  double result=gsl_integration_glfixed(&F,LowLimit,UpLimit,wf);
  gsl_integration_glfixed_table_free(wf);
  return result;
}    


// ********************************************************************

double gsl_integration3(int N, double (*function)(double, void *) ,void *p,double LowLimit,double UpLimit){
  gsl_integration_glfixed_table *wf =  gsl_integration_glfixed_table_alloc (N);
  gsl_function F;
  F.params   = p;  
  F.function = function;
  double result=gsl_integration_glfixed(&F,LowLimit,UpLimit,wf);
  gsl_integration_glfixed_table_free(wf);
  return result;
}    

// ********************************************************************

double gsl_integration2(double (*function)(double, void *) ,void *p,vector<double>XX, vector<double>WW){
  double result;
  int nn=WW.size();
  result=0;
  for(int i=0;i<nn;++i)result+=WW[i]*function(XX[i],p);
  return result;

}    


// ********************************************************************

void gsl_get_GL_weights(double LowLimit,double UpLimit, gsl_integration_glfixed_table *wf, vector<double> &XX, vector<double>&WW){
  double xi, wi;
  int nn=XX.size();
  omp_set_num_threads(omp_get_max_threads());

 // If this is already inside a parllelized loop, do not used it till you are
 // sure. It might lead to random errors!
//#pragma omp parallel for
  for(int i=0;i<nn;++i){
    gsl_integration_glfixed_point(LowLimit,UpLimit, i, &xi, &wi, wf);
    XX[i]=xi;
    WW[i]=wi;
  }
}    

// ********************************************************************

double gsl_integration_sin_kernel(double (*function)(double, void *), void *p, double omega, double LowLimit,double UpLimit){
  int NIT=1e3;
  double result, error;  
  double L=UpLimit-LowLimit;
  double relerr=0.1;      
  double abserr=0.01;

  gsl_integration_workspace  *w  =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_workspace  *wc =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,1000);   

  // La subrutina qawo integra una funcion f peseada con un sin(omega x), 
  // en un intervalo definido
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)

  gsl_function F;
  F.params   = p;  
  F.function = function; 
  result=0;
  for(;;){
    if(gsl_integration_qawo(&F,LowLimit,abserr,relerr,NIT,w,wf,&result,&error)){
      abserr*=10;
    }
    else{
      break;
    }
  }
  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (wc);
  gsl_integration_qawo_table_free(wf);
  return result;
}    




// ********************************************************************
// ********************************************************************

double gsl_integration_sin_kernel_lowlimit_inf(double (*function)(double, void *), void *p, double omega, double LowLimit){
  int NIT=1e3;
  double result, error;  
  double L=1;
  double relerr=1.0;      
  gsl_integration_workspace  *w  =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_workspace  *wc =  gsl_integration_workspace_alloc  (NIT); 
  gsl_integration_qawo_table *wf =  gsl_integration_qawo_table_alloc (omega,L,GSL_INTEG_SINE,1000);   

  // La subrutina qawo integra una funcion f peseada con un sin( omega x), 
  // el cual es en nuestro caso el que viene de j0(x)
  // con esto no necesito pasarle parametros al integrando pues en este caso f=kP(k)
  // Integramos en k para poder usar esta subrutina 
  // utiliza w = 1 con eta=kr, eta lo llamamos k arriba en las funcciones
  
  gsl_function F;
  F.params   = p;  
  F.function = function; 
  result=0;
  for(;;){
    if(gsl_integration_qawf(&F,LowLimit,relerr,1000,w,wc,wf,&result,&error)){
      relerr*=10;
    }
    else{
      break;
    }
  }
  gsl_integration_workspace_free (w);
  gsl_integration_workspace_free (wc);
  gsl_integration_qawo_table_free(wf);
  return result;
}    




// ********************************************************************
// ********************************************************************
double gsl_inter_pointers(double *x, double *y, int n, double xx){
  double ans;
  double xa[n];
  double ya[n];
  for(int i=0;i<n;i++)xa[i]=*(x+i);
  for(int i=0;i<n;i++)ya[i]=*(y+i);
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, xa, ya, n);
  ans=(xx<*x? *x : gsl_spline_eval (spline, xx, acc));
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return ans;
}


// ********************************************************************
// ********************************************************************

void gsl_bspline(vector<double>&xx, vector<double>&yy, vector<double>&new_xx, vector<double> &new_yy){
  int n=xx.size();
  int new_n=new_yy.size();
  double xo[n],yo[n];
  for(int i=0;i<n;i++){xo[i]=xx[i];yo[i]=yy[i];}
  double x_ini=xo[0];
  double x_fin=xo[n-1];
  gsl_matrix *X, *cov;
  gsl_vector *c, *w;
  gsl_vector *x, *y;
  int k = 4;   /*cubic spline*/
  int nbreak=new_n+2-k;
  gsl_bspline_workspace *bw;
  gsl_multifit_linear_workspace *mw;
  double chisq,Rsq,dof,tss;
  gsl_vector *B;
  X= gsl_matrix_alloc(n,new_n);
  cov= gsl_matrix_alloc(new_n,new_n);
  x = gsl_vector_alloc(n);
  y = gsl_vector_alloc(n);
  X = gsl_matrix_alloc(n, new_n);
  c = gsl_vector_alloc(new_n);
  w = gsl_vector_alloc(n);

  mw = gsl_multifit_linear_alloc(n, new_n);
  bw=gsl_bspline_alloc(k,nbreak);
  B=gsl_vector_alloc(new_n);

  gsl_bspline_knots_uniform(x_ini, x_fin,bw);

  for(int i=0;i<n;i++)gsl_vector_set(x,i,xo[i]);
  for(int i=0;i<n;i++)gsl_vector_set(y,i,yo[i]);

  
  for(int i=0;i<n;i++){
    gsl_vector_set(w,i,1e7); /**/
    double xi=gsl_vector_get(x,i);
    gsl_bspline_eval(xi,B,bw);
    for(int j=0;j<new_n;j++){
      double Bj=gsl_vector_get(B,j);
      gsl_matrix_set(X,i,j,Bj);
    }
  }
  gsl_multifit_wlinear(X,w,y,c,cov,&chisq,mw);
  {
    double xi, yi, yerr;
    for(int i=0;i<new_n; i++){
      xi=x_ini+i*(x_fin-x_ini)/(new_n-1);
      gsl_bspline_eval(xi, B, bw);
      gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
      new_xx[i]=xi;
      new_yy[i]=yi;
    }
  }
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);
  return; 
}

// ********************************************************************
// ********************************************************************
double gsl_inter(double *x, double *y, int n, double xn){
  double ans;
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  //  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_cspline, n);
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}
// ********************************************************************
// ********************************************************************

double gsl_inter_new(vector<double> &xx, vector<double> &yy, double xn){
  double ans;
  int n=xx.size();
  double x[n],y[n];

//  cout<<xn<<"  "<<xx[0]<<endl;
  for(int i=0;i<n;i++){x[i]=xx[i];y[i]=yy[i];}
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}

// ********************************************************************
// ********************************************************************

double gsl_inter_new2(vector<double> &xx, vector<vector<double> > &yy, int li, double xn){
  double ans;
  int n=xx.size();
  double x[n],y[n];
  for(int i=0;i<n;++i){x[i]=xx[i];y[i]=yy[li][i];}
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline    = gsl_spline_alloc (gsl_interp_linear, n);
  gsl_spline_init (spline, x, y, n);
  ans = gsl_spline_eval (spline, xn, acc);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  return  ans; 
}


// ****************************************************************************************************************************
// ****************************************************************************************************************************
void sort_index(int i, int j, int k, int *ii, int *jj, int *kk){
  vector<double> maa;
  maa.push_back(i);
  maa.push_back(j);
  maa.push_back(k);
  sort(maa.begin(), maa.end());
  *ii=maa.at(0);
  *jj=maa.at(1);
  *kk=maa.at(2);
  return ;
}


void matrix_inversion(vector< vector<double> > &G, vector< vector<double> > &y){

  int n=G[0].size();
  gsl_matrix *A  = gsl_matrix_alloc(n,n);
  gsl_matrix *Ai = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;++i)for(int j=0;j<n;++j)gsl_matrix_set(A, i,j,G[i][j]);
  gsl_permutation * p = gsl_permutation_alloc(n);
  int signum;

  gsl_linalg_LU_decomp(A, p, &signum);
  gsl_linalg_LU_invert(A, p, Ai);
  gsl_permutation_free(p);
  for (int i=0;i<n;++i)for(int j=0;j<n;++j)y[i][j]=gsl_matrix_get(Ai, i,j);
  gsl_matrix_free (A);
  gsl_matrix_free (Ai);    
  return ;
}


void matrix_det(vector< vector<double> > &G, double &det){
  int n=G[0].size();
  gsl_matrix *A  = gsl_matrix_alloc(n,n);
  for (int i=0;i<n;++i)for(int j=0;j<n;++j)gsl_matrix_set(A, i,j,G[i][j]);
  gsl_permutation * p = gsl_permutation_alloc(n);
  int signum;
  gsl_linalg_LU_decomp(A, p, &signum);
  det=gsl_linalg_LU_det(A, signum);
  gsl_permutation_free(p);
  gsl_matrix_free (A);
}


void get_eigen(vector<vector<double>>&icov_masses, vector<double>&masses){
  gsl_matrix *icova;
  gsl_vector *eig_vec, *eig;
  gsl_eigen_symm_workspace *workspace;
  int n_par=masses.size();
  icova = gsl_matrix_alloc(n_par,n_par);
  eig = gsl_vector_alloc(n_par);
  for(int i=0;i<n_par;i++)  for(int j=0;j<n_par;j++)gsl_matrix_set(icova, i,j, icov_masses[i][j]);
  workspace = gsl_eigen_symm_alloc(n_par);
  gsl_eigen_symm(icova, eig, workspace);
  gsl_sort_vector(eig);
  gsl_eigen_symm_free(workspace);
  for(int ip=0;ip<n_par;ip++)masses[ip]=sqrt(gsl_vector_get(eig,ip));
}


// ************************************************************************************
void get_det_matrix(vector<vector<double>>&matriz, double &determinant){
  vector<double>det(matriz.size(),0);
  get_eigen(matriz, det);
  determinant=1;         ;//Para escalar
  for(int i=0;i<matriz.size();i++)determinant*=det[i]; 
}


