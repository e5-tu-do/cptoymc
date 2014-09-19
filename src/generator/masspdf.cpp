// from ROOT
#include "TMath.h"
#include "Math/SpecFunc.h"

namespace cptoymc {
namespace generator {
namespace masspdf {

const double sq2pi = TMath::Sqrt(2.*TMath::ACos(-1.));
const double sq2pi_inv = 1./sq2pi;
const double logsq2pi = TMath::Log(sq2pi);
const double log_de_2 = TMath::Log(2.);

double low_x_BK(double nu,double x){
  return TMath::Gamma(nu)*TMath::Power(2.,nu-1.)*TMath::Power(x,-nu);
}


double low_x_LnBK(double nu, double x){
  return TMath::Log(TMath::Gamma(nu)) + (nu-1.)*log_de_2 - nu * TMath::Log(x);
}

double BK(double ni, double x) {
  double nu = TMath::Abs(ni);
  if ( x < 1.e-06 && nu > 0.) return low_x_BK(nu,x);
  if ( x < 1.e-04 && nu > 0. && nu < 55.) return low_x_BK(nu,x);
  if ( x < 0.1 && nu >= 55.) return low_x_BK(nu,x);

  //return gsl_sf_bessel_Knu(nu, x);
  return ROOT::Math::cyl_bessel_k(nu, x);
}

double LnBK(double ni, double x) {
  double nu = TMath::Abs(ni);
  if ( x < 1.e-06 && nu > 0.) return low_x_LnBK(nu,x);
  if ( x < 1.e-04 && nu > 0. && nu < 55.) return low_x_LnBK(nu,x);
  if ( x < 0.1 && nu >= 55.) return low_x_LnBK(nu,x);

  //return gsl_sf_bessel_lnKnu(nu, x);
  return TMath::Log(ROOT::Math::cyl_bessel_k(nu, x));
}

double LogEval(double d, double l, double alpha, double beta, double delta) {
  //double d = x-mu;
  //double sq2pi = TMath::Sqrt(2*TMath::ACos(-1));
  double gamma = alpha;//TMath::Sqrt(alpha*alpha-beta*beta);
  double dg = delta*gamma;
  double thing = delta*delta + d*d;
  double logno = l*TMath::Log(gamma/delta) - logsq2pi -LnBK(l, dg);
  
  return TMath::Exp(logno + beta*d +(0.5-l)*(TMath::Log(alpha)-0.5*TMath::Log(thing)) + LnBK(l-0.5,alpha*TMath::Sqrt(thing)));// + TMath::Log(TMath::Abs(beta)+0.0001) );
}

double diff_eval(double d, double l, double alpha, double beta, double delta){
  //double sq2pi = TMath::Sqrt(2*TMath::ACos(-1));
  //double cons1 = 1./sq2pi;
  double gamma = alpha;// TMath::Sqrt(alpha*alpha-beta*beta);
  double dg = delta*gamma;
  //double mu_ = mu;// - delta*beta*BK(l+1,dg)/(gamma*BK(l,dg));
  //double d = x-mu;
  double thing = delta*delta + d*d;
  double sqthing = TMath::Sqrt(thing);
  double alphasq = alpha*sqthing;
  double no = TMath::Power(gamma/delta,l)/BK(l,dg)*sq2pi_inv;
  double ns1 = 0.5-l;
  
  return no*TMath::Power(alpha, ns1)*TMath::Power(thing, l/2. - 1.25)*(-d*alphasq*(BK(l - 1.5, alphasq) + BK(l + 0.5, alphasq)) + (2.*(beta*thing + d*l) - d)*BK(ns1, alphasq))*TMath::Exp(beta*d)/2.;
}

double EvalIpatia(double x, double l, double zeta, double fb, double sigma, double mu, double a, double n, double a2, double n2) {
   double d = x-mu;
   double cons0 = TMath::Sqrt(zeta);
   double alpha, beta, delta, delta2, cons1, phi, A, B, k1, k2;
   double asigma = a*sigma;
   double a2sigma = a2*sigma;
   double out = 0.;
   if (zeta!= 0.) {
     phi = BK(l+1.,zeta)/BK(l,zeta); // careful if zeta -> 0. You can implement a function for the ratio, but carefull again that |nu + 1 | != |nu| + 1 so you jave to deal wiht the signs
     cons1 = sigma/TMath::Sqrt(phi);
     alpha  = cons0/cons1;//*TMath::Sqrt((1 - fb*fb));
     beta = fb;//*alpha;
     delta = cons0*cons1;
     
     if (d < -asigma){
       //printf("-_-\n");
       //printf("alpha %e\n",alpha);
       //printf("beta %e\n",beta);
       //printf("delta %e\n",delta);
       
       k1 = LogEval(-asigma,l,alpha,beta,delta);
       k2 = diff_eval(-asigma,l,alpha,beta,delta); 
       B = -asigma + n*k1/k2;
       A = k1*TMath::Power(B+asigma,n);
       //printf("k1 is %e\n",k1);
       //printf("k2 is %e\n",k2);
       //printf("A is%e\n",A);
       //printf("B is%e\n",B);
       out = A*TMath::Power(B-d,-n);
     }
     else if (d>a2sigma) {
       //printf("uoeo\n");
       k1 = LogEval(a2sigma,l,alpha,beta,delta);
       k2 = diff_eval(a2sigma,l,alpha,beta,delta);
       
       B = -a2sigma - n2*k1/k2;
       
       A = k1*TMath::Power(B+a2sigma,n2);
       
       out =  A*TMath::Power(B+d,-n2);
       
     }
     else {
       //printf("HERE\n");
       out = LogEval(d,l,alpha,beta,delta);
     }
     


   }
   else if (l < 0.) {
     beta = fb;
     cons1 = -2.*l;
     //delta = sigma;
     if (l<=-1.0) { delta = sigma *sqrt(-2+cons1);}
     else {
       printf("WARNING: zeta ==0 and l > -1 ==> not defined rms. Changing the meaning of sigma, but I keep fitting anyway\n");
       delta = sigma;
      
     }
     delta2 = delta*delta;
     if (d < -asigma ) {
       cons1 = TMath::Exp(-beta*asigma);
       phi = 1. + asigma*asigma/delta2;
       k1 = cons1*TMath::Power(phi,l-0.5);
       k2 = beta*k1- cons1*(l-0.5)*TMath::Power(phi,l-1.5)*2*asigma/delta2;
       B = -asigma + n*k1/k2;
       A = k1*TMath::Power(B+asigma,n);
       out = A*TMath::Power(B-d,-n);
     }
     else if (d > a2sigma) {
       cons1 = TMath::Exp(beta*a2sigma);
       phi = 1. + a2sigma*a2sigma/delta2;
       k1 = cons1*TMath::Power(phi,l-0.5);
       k2 = beta*k1+ cons1*(l-0.5)*TMath::Power(phi,l-1.5)*2.*a2sigma/delta2;
       B = -a2sigma - n2*k1/k2;
       A = k1*TMath::Power(B+a2sigma,n2);
       out =  A*TMath::Power(B+d,-n2);
       
     }
     else { out = TMath::Exp(beta*d)*TMath::Power(1. + d*d/delta2,l-0.5);}
   }
   else {
     //printf("zeta = 0 only suported for l < 0, while l = %e\n",0);
   }
   //printf("result is %e\n",out);
   return out;
 }

} // namespace masspdf
} // namespace generator
} // namespace cptoymc