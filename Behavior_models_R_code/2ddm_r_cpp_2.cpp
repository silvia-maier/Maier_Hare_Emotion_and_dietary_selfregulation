#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <ctime> 

using namespace Rcpp;
using namespace RcppParallel;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator


///-----------------------------------------------
///-----------------------------------------------

// [[Rcpp::depends(RcppParallel)]]
struct ddm2w : public Worker {
   
   // input vectors/matrices to read from
   const double d_v;
   const double d_h;
   const double thres;
   const double nDT;
   const double tIn;
   const double bias;
   const double vd;
   const double hd;
   GEN gen;
   
   // output vector to write to
   RVector<double> vecOut;
   
   // initialize from Rcpp input and output matrixes/vectors (the RMatrix/RVector class
   // can be automatically converted to from the Rcpp matrix/vector type)
   ddm2w(const double d_v, const double d_h, const double thres, const double nDT, const double tIn, const double bias, const double vd, const double hd, NumericVector vecOut, GEN gen)
      : d_v(d_v), d_h(d_h), thres(thres), nDT(nDT), tIn(tIn), bias(bias), vd(vd), hd(hd), vecOut(vecOut), gen(gen) {}
   
   // function call operator that work for the specified range (begin/end)
   	void operator()(std::size_t begin, std::size_t end) {
   		
   		double T = 5.2, dt = 0.008, lt;
		lt = (int)(T/dt);	
   		
   		std::vector<double> vec_tHealth(lt,1);
   		std::vector<double> vec_tVal(lt,1);
   		int aux = std::abs((int)(tIn/dt));
   		
   		if (tIn > 0) {
   			for (int t=0; t<aux; t++) {
				vec_tHealth[t] = 0;
			}
   		}
   		else if (tIn < 0) {
   			for (int t=0; t<aux; t++) {
				vec_tVal[t] = 0;
			}
   		}
		
////		std::cout << vec_tHealth[0] << "<--->" << vec_tVal[0] << '\n';
   		
   		
      	for (std::size_t i = begin; i < end; i++) {
	
			double X = bias;
			int flag = 0;
			double cont = 0;
			while (flag==0 && cont<lt) {
					
				X = X + (d_v*vd*vec_tVal[cont] + d_h*hd*vec_tHealth[cont])*dt + gen()*sqrt(dt); 
				if (X > thres) {
					flag=1;
					vecOut[i] = nDT + cont*dt;
				}
				else if (X < -thres) {
					flag=1;
					vecOut[i] = -nDT -cont*dt;
				}
				cont++;
			}
			
     	}
  	}
};


// [[Rcpp::export]]
NumericVector ddm2_parallel(double d_v, double d_h, double thres, double nDT, double tIn, double bias, double vd, double hd, unsigned int N) {
  
    const double sd_n = 1.4;
    ENG  eng;
    eng.seed(static_cast<unsigned int>(std::time(0)));
    DIST dist(0,sd_n);
    GEN  gen(eng,dist);
  
	//output vector
   	NumericVector vecOut(N);

   	// create the worker
   	ddm2w ddm2w(d_v, d_h, thres, nDT, tIn, bias, vd, hd, vecOut, gen);
     
   	// call the worker
   	parallelFor(0, N, ddm2w);

	// vector with RTs is returned to R
   	return vecOut;
}


