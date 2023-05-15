
/*fuction implementation*/
void exchange(void *var1, void *var2){
  char temp;
  int varsize = sizeof(*var1);
  for (int i=0; i<varsize; i++){
    temp = *((char *) var1 + i);
    *((char *) var1 + i) = *((char *) var2 +i);
    *((char *) var2 + i) = temp;
  }
}
