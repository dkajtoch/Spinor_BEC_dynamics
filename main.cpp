#include "evolution.hpp"

#define PI boost::math::constants::pi<double>()

int main() {

   timeEvolution<20> system("symmetric");
   state_type x = system.propagate(0.0,0.1);

   system.saveAll(0.0, (PI/2.0), 500, "symmetric_N=20.dat");

   return 0;
}
