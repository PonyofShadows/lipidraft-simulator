#include "src/lipidraft.h"


#define TEMP 300
#define A_LENGTH 2

/* variables */
static double test_ca;
static double test_blength;
static double test_occupancy;


int main(){
  sprintf(database_name,"test_lb_plus");
  for (test_ca=0.15;test_ca<=0.15;test_ca+=0.01){
    for (test_blength=1.5;test_blength>=1.0;test_blength-=0.1){
      for (test_occupancy=0;test_occupancy<=0.5;test_occupancy+=0.04){
        oneexp(TEMP, test_ca, test_occupancy, A_LENGTH, test_blength);
      }
    }  
  }
  return 0;
}





