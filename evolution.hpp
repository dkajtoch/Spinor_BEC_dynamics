#ifndef _EVOLUTION_HPP_
#define _EVOLUTION_HPP_

// read necessary libraries
#include <complex>
#include <cmath>
#include <vector>
#include <array>
#include <fstream>
#include <string>
#include <boost/numeric/odeint.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>

namespace spinor {

typedef std::vector< std::complex<double> > state_type;
constexpr std::complex<double> re {1.0, 0.0};
constexpr std::complex<double> im {0.0, 1.0};
constexpr double PI = boost::math::constants::pi<double>();

/* ---------------------------------
 * Definition of the norm
 * --------------------------------- */
double norm( const state_type& x )
{
   double sum = 0.0;

   for( auto xx : x )
      sum += boost::math::pow<2>( std::abs( xx ) );

   return sum;
}

/* -------------------------------
 * Initial state definition
 * ------------------------------- */
state_type Fock_Mzero( std::size_t N, std::size_t l )
{
   state_type x(N/2+1);
   x.at(l)= {1.0, 0.0};

   return x;
}

/* -----------------------------------------------------
 * Total mean number of particles in modes m_f = +1, -1
 * ----------------------------------------------------- */
double meanNs( const state_type& psi )
{
   std::size_t dim = psi.size();

   double mean  = 0.0;

   for( std::size_t i=0; i<dim; i++ ) 
   {
      mean += boost::math::pow<2>( std::abs(psi[i]) ) * 2.0 * static_cast<double>(i);
   }

   return mean;
}

/* -----------------------------------------------
 * Variance of total particle number in side modes
 * ----------------------------------------------- */
double varNs( const state_type& psi )
{
   std::size_t dim = psi.size();

   double mean  = meanNs(psi);
   double mean2 = 0.0;

   for( std::size_t i=0; i<dim; i++ ) 
   {
      mean2 += boost::math::pow<2>( std::abs(psi[i]) ) * 
         4.0 * boost::math::pow<2>( static_cast<double>(i) );
   }

   return mean2 - mean*mean;
}

/* -------------------------- */
/* Quantum Fisher Information */
/* -------------------------- */
std::array<double,4> quantumFisherInf( const state_type &psi )
{
   std::size_t dim = psi.size();
   std::size_t N   = 2 * (dim-1);

   double Gamma55 = 0.0;
   double meanY   = 0.0;
   double varY    = 0.0;
   std::complex<double> corre {0.0, 0.0};
   double numberTerms = 0.0;

   for( std::size_t i=0; i<dim; i++)
   {
      // covariance matrix element 'Dxy' or 'Qxy'
      Gamma55 += std::pow( std::abs(psi[i]), 2 ) * 
         2.0 * static_cast<double>(i) * static_cast<double>(i+1) ;

      // covariance matrix element 'Y'
      meanY += std::pow( std::abs(psi[i]), 2 ) * 
         4.0/3.0 * ( 3.0 * static_cast<double>(i) - static_cast<double>(N) );

      varY += std::pow( std::abs(psi[i] ), 2 ) * 
            std::pow( 4.0/3.0 * ( 3.0 * static_cast<double>(i) - static_cast<double>(N) ) , 2 );

      // covariance matrix element Jx, Qyz
      numberTerms += std::pow( std::abs(psi[i]), 2 ) * 
               ( static_cast<double>(N - i) + 
                    2.0 * static_cast<double>(i) * static_cast<double>(N-2*i) );

      if( i>0 )
      {
         corre += std::conj( psi[i] ) * psi[i-1] * 
                 static_cast<double>(i) * 
                 std::sqrt( static_cast<double>(N-2*i+2) * static_cast<double>(N-2*i+1) );
      }

   }

   double Gamma11 = numberTerms + 2.0 * std::real( corre );
   double Gamma22 = numberTerms - 2.0 * std::real( corre );
   double Gamma12 = 2.0 * std::imag( corre );

   double Gamma77 = varY- meanY * meanY;

   // return 3 different interferometers
   std::array<double, 4> lambda;
   lambda[0] = 0.5 * (Gamma11 + Gamma22) + 
      0.5 * std::sqrt( std::pow( Gamma11 - Gamma22, 2 ) + 4.0 * std::pow( Gamma12, 2 ) );
   lambda[1] = Gamma22;
   lambda[2] = Gamma55;
   lambda[3] = Gamma77;

   return lambda;
}

/* ----------------------------------
 * Equations of motion
 * ---------------------------------- */

// ---------- 1st Policy ------------
class symmetricDynamics {

public:

   std::size_t _N;   // number of atoms

   // constructor
   symmetricDynamics( std::size_t N )
      : _N{N}
   {}

   void operator() ( const state_type& x, state_type& dxdt, const double /* t */ )
   {

      for( std::size_t i=0; i< (_N/2+1); i++ )
      {
         dxdt[i] = -im * 2.0 * static_cast<double>(i) *
            ( 2.0*static_cast<double>(_N-2*i) - 1.0 ) * x[i];

         if( i<_N/2 )
         {
            dxdt[i] += -im * 2.0 * std::sqrt( static_cast<double>(_N - 2*i) *
                                    static_cast<double>(_N - 2*i - 1) *
                                    static_cast<double>(i + 1) * 
                                    static_cast<double>(i + 1)
                                    ) * x[i+1];
         }

         if( i>0 )
         {
            dxdt[i] += -im * 2.0 * std::sqrt( static_cast<double>(_N - 2*i + 2) *
                                       static_cast<double>(_N - 2*i + 1) *
                                       static_cast<double>(i) * 
                                       static_cast<double>(i)
                                    ) * x[i-1];
         }
      }
   }

};

// ---------- 2nd Policy ------------
class resonanceDynamics {

public:

   std::size_t N;

   // constructor
   resonanceDynamics( std::size_t N ) 
      : N{N}
   {}
   
   void operator() ( const state_type& x, state_type& dxdt, const double /* t*/ )
   {
      for( std::size_t i=0; i< N/2+1; i++ )
      {
         if( i == 0 )
         {
            dxdt[i] = -im * 2.0 * std::sqrt( static_cast<double>(N - 2*i) *
                                    static_cast<double>(N - 2*i - 1) *
                                    static_cast<double>(i + 1) * 
                                    static_cast<double>(i + 1)
                                    ) * x[i+1];
         }

         if( i == N/2 )
         {
            dxdt[i] = -im * 2.0 * std::sqrt( static_cast<double>(N - 2*i + 2) *
                                       static_cast<double>(N - 2*i + 1) *
                                       static_cast<double>(i) * 
                                       static_cast<double>(i)
                                    ) * x[i-1];
         }

         if( i>0 && i<N/2 )
         {
            dxdt[i] = -im * 2.0 * std::sqrt( static_cast<double>(N - 2*i) *
                                    static_cast<double>(N - 2*i - 1) *
                                    static_cast<double>(i + 1) * 
                                    static_cast<double>(i + 1)
                                    ) * x[i+1]

                      -im * 2.0 * std::sqrt( static_cast<double>(N - 2*i + 2) *
                                       static_cast<double>(N - 2*i + 1) *
                                       static_cast<double>(i) * 
                                       static_cast<double>(i)
                                    ) * x[i-1];
         }
      }
   }
};

// specialization
template< class systemType >
class quantumSystem {

public:

   systemType s;

   // constructor
   quantumSystem( std::size_t N )
      : s{ systemType(N) }
   {}

   // evolution of the quantum state from t0 till tk
   void propagate( state_type& x, double t0, double tk, double dt=1.0E-04 )
   {
      using namespace boost::numeric::odeint;

      // numerical integration
      runge_kutta4< state_type > stepper;

      integrate_const( stepper, std::ref(s), x, t0, tk, dt );
   }

   // integrate the system from time 't' until time t = t + n * dt
   // define observer class to collect the data
   template< class F >
   void propagate( state_type& x, double t0, double tk, F Observer, std::size_t n=100, double dt=1.0E-04 )
   {
      using namespace boost::numeric::odeint;

      // numerical integration
      runge_kutta4< state_type > stepper;

      double t = t0;
      while( t <= tk )
      {
         // collect data
         Observer( x, t );

         integrate_n_steps( stepper, std::ref(s), x, t, dt, n );
         t += n * dt;
      }
   }
};

/* ----------------------------------
 * Write all quantities into a file
 * ---------------------------------- */
class writeAll {

public:

   // constructor
   writeAll( std::string file_name ) 
      : file{ std::ofstream( file_name, std::ofstream::out ) }
   {
      file << "__time__meanNs__varNs__optFq2x2__varQyz__varDxy__varY__" << std::endl;
   }


   void operator() ( const state_type& x, const double t )
   {
      // calculate physical quantities
      mean = meanNs(x);
      var  = varNs(x);
      fis  = quantumFisherInf(x);

      file << t << "\t"
           << mean << "\t"
           << var << "\t"
           << fis[0] << "\t"
           << fis[1] << "\t"
           << fis[2] << "\t"
           << fis[3] << std::endl;
   }

private:

   // variables
   std::ofstream file;
   // physical quantities
   double mean;
   double var;
   std::array<double,4> fis;
};

} // namespace spinor
#endif


