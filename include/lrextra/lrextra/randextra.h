#include <stdlib.h>
#include <time.h>


/*function implementation*/
void srandtime(){
  time_t t;
  srand((int) time(&t));
}

double randporb(){
  return (((double)rand())/((double)(RAND_MAX)+1)); // [0,1)
}

int randint(int min,int max){
  return ((rand() % (max-min))+min);
}

double randpolar(){
  if (randint(0,2) == 0)
    return -1;
  else
    return 1;
}

double randdouble(double min, double max){
  return (min + randporb() * (max-min));
}

double randdiscretedouble(double min, double max, int n){
return (
  ((double)randint(0,n)) / ((double)n) *
  (max-min) + min
    );
}

double randdoublestep(double angle, double min, double max, int n){
  double outputnum = angle + randpolar()* (max - min) / ((double)n);
  if (outputnum < min) outputnum = min;
  else if (outputnum > max) outputnum = max;
  return outputnum;

}
