#include "evolution.hpp"

using namespace spinor;

int main() {

   const std::size_t N = 200;

   quantumSystem<symmetricDynamics> sys(N);
   state_type x = Fock_Mzero(N,0);

   sys.propagate( x, 0.0, 0.1, writeAll("symmetric_N=200.dat"), 1000, 1.0E-06 );

   return 0;
}
