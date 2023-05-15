/* macros */
//constant
#define PI 3.141592653589793
#define KB 1.38064853 // (*10^(-23) J/K)
/* function implementations*/
double limitedchange(double value, double change,
    double min, double max){
  double changedvalue = value + change;
  if (changedvalue > max)
    changedvalue = max;
  else if (changedvalue < min)
    changedvalue = min;
  return changedvalue;
}

void cyclestepindex(int *stepedvalue, int max){
  if (*stepedvalue == -1) *stepedvalue = max-1;
  else if (*stepedvalue == max) *stepedvalue = 0; 
}
