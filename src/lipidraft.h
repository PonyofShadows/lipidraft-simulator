#include <stdlib.h>
#include <stdio.h> 
#include <stdbool.h>
#include <math.h>
#include <string.h>

// needs mongodb, libmongoc and libbson. See As:
// https://www.mongodb.com/
// https://github.com/mongodb/mongo-c-driver
#include <libmongoc-1.0/mongoc.h>
#include <libbson-1.0/bson.h>

#include "../include/lrextra/lrextra.h"

/* macros */
#define LATTICE_SIZE 32
#define CHOL_LENGTH 2
#define DB 8
#define INTERATION_MAX 10000000
/* enums */
enum NAME { Ahead, Atail, Bhead, Btail, C, Vacancy};
enum STEP { Begin , End};
/* data structures */
typedef struct{
  double varphi, theta;
} Sphercoord;

typedef struct {
  enum NAME pl;
  Sphercoord sc_pl;
  enum NAME chol0, chol1;
  Sphercoord sc_chol0, sc_chol1;
} Cell;

typedef struct{
  Cell matrix[LATTICE_SIZE][LATTICE_SIZE];
} Lattice;

typedef struct{
  double x;
  double y;
  double z;
} Coord;

typedef struct{
  enum NAME names[DB];
  Coord recoord[DB];
} Pairdb;

typedef void (*FuncPtr)(void);
/* function declarations */
bool test;
static void initlattice();

static void randomoffset();
static double randomchangeangle(double angle, double min, double max);
static void trywalk();

static double energybetweencells();
static double energyofchangedpart();
static double energyofall();

static void trial();

static double handlelj(enum NAME n1, enum NAME n2,
    Coord *c1, Coord *c2 );

static double handlepairdb();

void printlattice();

/* set */
// Ahead, Atail, Bhead, Btail, C, Vacancy
static double epsilon[6][6] = {
  {    12,    10,      8,      8,   18, 0,},
  {    10,    11,      8,      8,   18, 0,},
  {     8,     8,     12,     10,    6, 0,},
  {     8,     8,     10,     11,    6, 0,},
  {    18,    18,      6,      6,   10, 0,},
  {     0,     0,      0,      0,    0, 0,},
};
// Ahead, Atail, Bhead, Btail, C, Vacancy
static double rmin[6][6]= {
  { 0.9,   1,   1, 1.2, 0.9,1000,},
  {   1, 0.8,   1,   1, 0.9,1000,},
  {   1,   1, 0.9,   1,   1,1000,},
  { 1.2,   1,   1, 0.8,   1,1000,},
  { 0.9, 0.9,   1,   1, 0.9,1000,},
  {1000,1000,1000,1000,1000,1000,},
};
/* variables */
static double temp = 300;
static double beta;

static double ca, occupancy;
static double a_length, b_length;
static double head_distance;

static const double varphimin = 0;
static const double varphimax = PI / 2;
static const double thetamin = - PI;
static const double thetamax = PI;
static const double anglestep = PI / 512;

//static double sin30 = 0.5;
static double cos30=cos(PI/6);
static double tan30=tan(PI/6);

static Lattice lattice;
static double energychange=0;
static double initenergy=0;
static double currentenergy=0;

static int choice;
static bool trystate;
static Cell cur;
static Cell another;
static int curx, cury;
static bool offset;
static int anotherx, anothery;

static int targetx, targety;
static int neighborx, neighbory;
static Cell* target=NULL;
static Cell* neighbor=NULL;

static Pairdb pairdb;

// mongodb
static const char *uri_string = "mongodb://localhost:27017";
static mongoc_uri_t *uri = NULL;
static mongoc_client_t *client = NULL;
static mongoc_database_t *database = NULL;
static mongoc_collection_t *collection = NULL;

static char database_name[16] = "test";
static const char *collection_name_head = "exp_00_";
static int collection_count=0;
//bson
static bson_error_t error;
static bson_oid_t oid;
static bson_t *parameters = NULL, *datas = NULL, *iterationdata;

// for save
static char pltype[2];
static uint32_t numpl=0,numchol=0;
static uint32_t x,y;
static double plx, ply;
static double chol0x, chol0y, chol1x, chol1y;
static const char *key,*key1,*key2;
static char keybuf[48],keybuf1[48],keybuf2[48];
static size_t keylen,keylen1,keylen2;

static char stepstr[6];
static bson_t child0, child1, child2;
static bson_t grandchild1[5],grandchild2[4];

static int iteration;
/* function implementations */
static void initlattice(){
   double n = (double) (pow(LATTICE_SIZE, 2));
   double n2 = 2 * n;
   double na = n * ca;
   double nc = n2 * occupancy;
   for (int x=0;x<LATTICE_SIZE;x++){
     for (int y=0;y<LATTICE_SIZE;y++){
       //A & B
       if (randporb()< na/n) {
          lattice.matrix[x][y].pl = Ahead;
          na -= 1;
          n -= 1;
        } else {
          lattice.matrix[x][y].pl = Bhead;
          n -= 1;
        }
       //lattice.matrix[x][y].sc_pl.varphi = randdouble(varphimin,varphimax);
       //lattice.matrix[x][y].sc_pl.theta = randdouble(thetamin,thetamax);
       // C
        if (randporb()<nc/n2){
          lattice.matrix[x][y].chol0 = C;
          nc -= 1;
          n2 -= 1;
          //lattice.matrix[x][y].sc_chol0.varphi = randdouble(varphimin,varphimax);
          //lattice.matrix[x][y].sc_chol0.theta = randdouble(thetamin,thetamax);
        } else {
          lattice.matrix[x][y].chol0 = Vacancy;
          n2 -= 1;
          //lattice.matrix[x][y].sc_chol0.varphi = 0;
          //lattice.matrix[x][y].sc_chol0.theta = 0;
        }
        if (randporb()< nc/n2){
          lattice.matrix[x][y].chol1 = C;
          nc -= 1;
          n2 -= 1;
          //lattice.matrix[x][y].sc_chol1.varphi = randdouble(varphimin,varphimax);
          //lattice.matrix[x][y].sc_chol1.theta = randdouble(thetamin,thetamax); 

        } else {
          lattice.matrix[x][y].chol1 = Vacancy;
          n2 -= 1;
          //lattice.matrix[x][y].sc_chol0.varphi = 0;
          //lattice.matrix[x][y].sc_chol0.theta = 0;
        }
    }
   }
}
static void toneighbor(int r){
  switch(r){
    case 0: neighborx = targetx; neighbory = targety; break;
    case 1: neighborx=targetx; neighbory = targety+1; break;
    case 2: neighborx=targetx+1; neighbory = targety; break;
    case 3: neighborx=targetx-1; neighbory=targety+1; break;
    case 4: neighborx=targetx; neighbory = targety-1; break;
    case 5: neighborx=targetx-1; neighbory = targety; break;
    case 6: neighborx=targetx+1; neighbory=targety-1; break;
  }
  cyclestepindex(&neighborx,LATTICE_SIZE); 
  cyclestepindex(&neighbory,LATTICE_SIZE);
}

static void randomoffset(){
  targetx = curx;
  targety = cury;
  toneighbor(randint(1,7));
  anotherx = neighborx;
  anothery = neighbory;
  another = lattice.matrix[anotherx][anothery];
}
static double randomchangeangle(double angle, double min, double max){
  double outputangle;
  //if (randint(0,2) == 1)
    //outputangle = limitedchange(angle,anglestep,min,max);
    //outputangle = randdoublestep(angle,min,max,64);
  //else
    outputangle = randdouble(min,max);
    //outputangle = randdiscretedouble(min,max,64);
  return outputangle;
}

static void trywalk(){
  redo:;
  curx = randint(0, LATTICE_SIZE);
  cury = randint(0, LATTICE_SIZE);
  cur = lattice.matrix[curx][cury];
  //try
 choice = randint(0,5);
  switch (choice){
    case 0://exchange (pl, sc_pl) with another cell
      randomoffset();
      offset = true;
      // for fast
      if (cur.pl == another.pl) goto redo;
      exchange(&(cur.pl),&(another.pl));
      exchange(&(cur.sc_pl),&(another.sc_pl));
      break;
    case 1://exchange (chol0, sc_chol0) with another cell
      randomoffset();
      offset = true;
      if (cur.chol0==Vacancy && another.chol0==Vacancy) goto redo;
      exchange(&(cur.chol0),&(another.chol0));
      exchange(&(cur.sc_chol0),&(another.sc_chol0));
      break;
    case 2://change sc_pl
      offset = false;
      if (randint(0,2) == 1)
        cur.sc_pl.varphi = randomchangeangle(cur.sc_pl.varphi,varphimin,varphimax);
      else
        cur.sc_pl.theta = randomchangeangle(cur.sc_pl.theta,thetamin,thetamax);
      break;    
    case 3://change sc_chol0
      if (cur.chol0 == Vacancy) goto redo;
      offset = false;
      if (randint(0,2) == 1)
        cur.sc_chol0.varphi = randomchangeangle(cur.sc_chol0.varphi,varphimin,varphimax);
      else
        cur.sc_chol0.theta = randomchangeangle(cur.sc_chol0.theta,thetamin,thetamax);
      break;    
    case 4://exchange chol0 and chol1, so the sc_chols.
      if (cur.chol0==Vacancy && cur.chol1==Vacancy) goto redo;
      offset = false;
      exchange(&(cur.chol0),&(cur.chol1));
      exchange(&(cur.sc_chol0),&(cur.sc_chol1));
      break;
}
}
static double handlelj(enum NAME n1, enum NAME n2,
    Coord *c1, Coord *c2 ){
  double r2;
 // Vacancy
 if (n1==Vacancy || n2==Vacancy)
   return 0;
 else if ((n1==Ahead && n2==Bhead) || (n1==Bhead && n2==Ahead))
    return (- epsilon[0][2]);
 else if (n1==Ahead && n2==Ahead)
    return (- epsilon[0][0]);
 else if (n1==Bhead && n2==Bhead)
    return (- epsilon[2][2]);
 // others 
 r2= pow(c2->x - c1->x,2) +pow(c2->y - c1->y,2) +
   pow(c2->z - c1->z,2);
 if (r2 < 0.01) // don't be close to zero
  r2 = 0.01;
  //printf("#############################\n");
  //printf("c1: x:%lf, y:%lf, z: %lf\n",c1->x,c1->y,c1->z);
  //printf("c2: x:%lf, y:%lf, z: %lf\n",c2->x,c2->y,c2->z);
 //for (int i=0;i<=7;i++)
  //printf("handledb[%d]: x:%lf, y:%lf, z: %lf\n",i,pairdb.recoord[i].x,pairdb.recoord[i].y,pairdb.recoord[i].z);

 /*
  * test
  * */
  //double output = (epsilon[n1][n2]*
    // (pow(pow(rmin[n1][n2],2)/r2,6)-2*pow(pow(rmin[n1][n2],2)/r2,3)));
   //printf("n1=%d,n2=%d,epsilon=%lf,r^2=%lf,energy=%lf\n",
     // n1,n2,epsilon[n1][n2],r2,output);
  //return output;
 //
 return (epsilon[n1][n2]*
     (pow(pow(rmin[n1][n2],2)/r2,6)-2*pow(pow(rmin[n1][n2],2)/r2,3)));
}

static double handlepairdbouter(){
  double energy=0;
  for (int i=0;i<DB-4;i++){
    for (int j=4;j<DB;j++){
      //printf("ij: (%d,%d)",i,j);
     energy += handlelj(pairdb.names[i],pairdb.names[j],
         &(pairdb.recoord[i]),&(pairdb.recoord[j]));
    }
  //printf("\n");
  }
  //printf("handlepairdbenergy:%lf\n",energy);
  return energy;
}

static double handlesinglecell(){
  double energy=0;
  energy += handlelj(pairdb.names[2],pairdb.names[3],&(pairdb.recoord[2]),&(pairdb.recoord[3]));
  energy += handlelj(pairdb.names[0],pairdb.names[2],&(pairdb.recoord[0]),&(pairdb.recoord[2]));
  energy += handlelj(pairdb.names[0],pairdb.names[3],&(pairdb.recoord[0]),&(pairdb.recoord[3]));
  energy += handlelj(pairdb.names[1],pairdb.names[2],&(pairdb.recoord[1]),&(pairdb.recoord[2]));
  energy += handlelj(pairdb.names[1],pairdb.names[3],&(pairdb.recoord[1]),&(pairdb.recoord[3]));
  return energy;
}
static double energybetweencells(){
  double varphi,theta;
  double targetlength,neighborlength;
  /*
   * test
   * */
//  printf("trystate:%d,offset:%d,(%d,%d)(%d,%d)\n",
 //     trystate,offset,targetx,targety,neighborx,neighbory);
  //names
  pairdb.names[0] = target->pl;
  pairdb.names[2] = target->chol0;
  pairdb.names[3] = target->chol1;
  pairdb.names[4] = neighbor->pl;
  pairdb.names[6] = neighbor->chol0;
  pairdb.names[7] = neighbor->chol1;
  if (target->pl == Ahead) {
    pairdb.names[1] = Atail;
    targetlength = a_length;
  }
  else {
    pairdb.names[1] = Btail;
    targetlength = b_length;
  }
  if (neighbor->pl == Ahead){
    pairdb.names[5] = Atail;
    neighborlength = a_length;
  }
  else {
    pairdb.names[5] = Btail;
    neighborlength = b_length;
  }
  //recoord
  //target
  pairdb.recoord[0].x = 0;
  pairdb.recoord[0].y = 0;
  pairdb.recoord[0].z = 0;

  varphi = target->sc_pl.varphi;
  theta = target->sc_pl.theta;
  //printf("###varphi:%lf,theta:%lf\n",varphi,theta);
  pairdb.recoord[1].x = targetlength/2*sin(varphi)*cos(theta);
  pairdb.recoord[1].y = targetlength/2*sin(varphi)*sin(theta);
  pairdb.recoord[1].z = targetlength/2 * cos(varphi);
  //printf("###xyz=(%lf,%lf,%lf)\n",pairdb.recoord[1].x,pairdb.recoord[1].y,pairdb.recoord[1].z);
  if (target->chol0 == C){
  varphi = target->sc_chol0.varphi;
  theta = target->sc_chol0.theta;
  pairdb.recoord[2].x = 0.5*head_distance + CHOL_LENGTH/2*sin(varphi)*cos(theta);
  pairdb.recoord[2].y = 0.5*head_distance * tan30 + CHOL_LENGTH/2*sin(varphi)*sin(theta);
  pairdb.recoord[2].z = CHOL_LENGTH/2 * cos(varphi);
}
  if (target->chol1 == C){ 
  varphi = target->sc_chol1.varphi;
  theta = target->sc_chol1.theta;
  pairdb.recoord[3].x = 0.5*head_distance + CHOL_LENGTH/2*sin(varphi)*cos(theta);
  pairdb.recoord[3].y = -( 0.5*head_distance * tan30) + CHOL_LENGTH/2*sin(varphi)*sin(theta);
  pairdb.recoord[3].z = CHOL_LENGTH/2 * cos(varphi);
  }
  //neighbor
  pairdb.recoord[4].y = cos30 * head_distance*((double)(neighbory - targety));
  pairdb.recoord[4].x = tan30 * pairdb.recoord[4].y + head_distance*((double)(neighborx - targetx));
  pairdb.recoord[4].z = 0;
  
  varphi = neighbor->sc_pl.varphi;
  theta = neighbor->sc_pl.theta;
  pairdb.recoord[5].x = pairdb.recoord[4].x + neighborlength/2*sin(varphi)*cos(theta);
  pairdb.recoord[5].y = pairdb.recoord[4].y + neighborlength/2*sin(varphi)*sin(theta);
  pairdb.recoord[5].z = neighborlength/2*cos(varphi);
  
  if (neighbor->chol0 == C){
  varphi = neighbor->sc_chol0.varphi;
  theta = neighbor->sc_chol0.theta;
  pairdb.recoord[6].x = pairdb.recoord[4].x + head_distance/2 +
      CHOL_LENGTH/2*sin(varphi)*cos(theta);
  pairdb.recoord[6].y = pairdb.recoord[4].y + head_distance/2*tan30+ 
      CHOL_LENGTH/2*sin(varphi)*sin(theta);
  pairdb.recoord[6].z = CHOL_LENGTH/2*cos(varphi);
  }

if (neighbor->chol1 == C){
  varphi = neighbor->sc_chol1.varphi;
  theta = neighbor->sc_chol1.theta;
  pairdb.recoord[7].x = pairdb.recoord[4].x + head_distance/2 +
      CHOL_LENGTH/2*sin(varphi)*cos(theta);
  pairdb.recoord[7].y = pairdb.recoord[4].y - head_distance/2*tan30+ 
      CHOL_LENGTH/2*sin(varphi)*sin(theta);
  pairdb.recoord[7].z = CHOL_LENGTH/2*cos(varphi);
}
  return handlepairdbouter();
} 

static double energyofchangedpart(){
  double energy=0;
  // cur
  targetx = curx;
  targety = cury;

  if (trystate)
    target = &cur;
  else 
    target = &(lattice.matrix[targetx][targety]);
  for (int i=1;i<7;i++){
    toneighbor(i);
    if (trystate && offset && neighborx==anotherx && neighbory==anothery)
      neighbor = &another;
    else 
      neighbor = &(lattice.matrix[neighborx][neighbory]);
   // printf("0nowchangeenergy:%lf\n",energy);
    energy += energybetweencells();
    //printf("11111111111111nowchangeenergy:%lf\n",energy);
  }
    energy += handlesinglecell();
    //printf("2nowchangeenergy:%lf\n",energy);
  // another
  if (offset){
    targetx = anotherx;
    targety = anothery;
    if (trystate)
      target = &another;
    else
      target = &(lattice.matrix[targetx][targety]);
    for (int i=1;i<7;i++){
      toneighbor(i);
      if (neighborx==curx && neighbory == cury) continue;
      neighbor = &(lattice.matrix[neighborx][neighbory]);
      energy += energybetweencells();
    //printf("3nowchangeenergy:%lf\n",energy);
    }
    energy += handlesinglecell();
    //printf("4nowchangeenergy:%lf\n\n",energy);
  }
    //if (trystate) printf("---------------energy1:%lf\n",energy);
    //else printf("---------------energy0:%lf\n",energy);
  return energy;
}
static double energyofall(){
  double energy = 0;
  /*
   * Be careful! init energy maybe too high
   * */
  trystate = false;
  for (targetx=0;targetx<LATTICE_SIZE;targetx++){
    for (targety=0;targety<LATTICE_SIZE;targety++){
      target = &(lattice.matrix[targetx][targety]);
      for (int i=1;i<=3;i++){
        toneighbor(i);
        neighbor = &(lattice.matrix[neighborx][neighbory]);
        energy += energybetweencells();
    }
        energy += handlesinglecell();
  }
}
  return energy;
}

static void trial(){
  double energy0, energy1;
  energychange = 0;
  // energychange
  trystate = false;
  energy0 = energyofchangedpart();
  printf("w");
  trystate = true;
  energy1 = energyofchangedpart();
  printf("\b");
  energychange = energy1 - energy0;
  //trial
  if ( (energychange <= 0) ||
      (randporb() < exp(-energychange))){
    lattice.matrix[curx][cury] = cur;
    if (offset)
      lattice.matrix[anotherx][anothery] = another;
   currentenergy += energychange;
    //printf("choice:%d, energechange: %lf,",choice,energychange);
    //printf("current energy:%lf\n",currentenergy);
    //printf("energyofall:%lf\n",energyofall());
    //usleep(10000);
  } else
        energychange = 0;
  //printf("currentenergy:%lf, choice:%d, energechange: %lf\n",
  //    currentenergy, choice,energychange);
  //current energy
}


/* run */
void init(
    double in_temp,
    double in_ca, double in_occupancy ,
    double in_alength, double in_blength
    ){
  // turn user's input to static variable
  temp = in_temp;
  beta = 1/(KB*temp);
  ca = in_ca;
  occupancy = in_occupancy;
  a_length = in_alength;
  b_length = in_blength;
  head_distance = ca*rmin[0][0] + (1-ca)*rmin[2][2];
  printf("ca:%.2lf,oc:%.2lf,bl:%.2lf\n",ca,occupancy,b_length);
  srandtime();
  initlattice();
  initenergy = energyofall();
  //initenergy = 0;
  currentenergy = initenergy;
  energychange = 0;
}

void stepfunc(unsigned localiteration){
    for (int i=0;i<localiteration;i++){
      //printf("before:%lf\n",energyofall(&lattice));
      trywalk();
      trial();
      //printf("choice:%d, currentenergy:%lf,energychange:%lf\n",choice,currentenergy,energychange);
      //if (currentenergy != energyofall(&lattice)) abort();
      //if (energychange < -100000) abort();
  }
}

void stepuntilbalanced(){
  iteration = 400000;
  double beforeenergy;
  printf("simulation %d processing...\n",collection_count);
  stepfunc(iteration);
  while (iteration < INTERATION_MAX){
    beforeenergy = currentenergy;
    stepfunc(5000);
    iteration+=5000;
    //printf("step for %d,choice:%d,currentenergy:%lf\n",
      //  iteration,choice,currentenergy);
    if (abs((beforeenergy-currentenergy)) <= 0.0000001){
    //if (beforeenergy == currentenergy){
      printf("balanced.\n");
      return;
    }
  }
  printf("unbalanced.\n");
  iteration = -1;
  return;
}

void printlattice(){
  printf("the lattice looks like:\n");
  int local_c=0;
  for (y=0;y<LATTICE_SIZE;y++){
    for (x=0;x<LATTICE_SIZE;x++){
      if (lattice.matrix[x][y].pl == Ahead) printf("A");
      else printf("B");
      if (lattice.matrix[x][y].chol0 == C) local_c +=1;
      if (lattice.matrix[x][y].chol1 == C) local_c +=1;
    }
    printf("\n");
  }
    printf("number of C: %d\n",local_c);
}

/* save to mariadb*/
int initdb(){
  /*
    * Safely create a MongoDB URI object from the given string
    */
  char namebuffer[16];
  char collection_name[24];
  mongoc_init();
  uri = mongoc_uri_new_with_error(uri_string, &error);
   if (!uri) {
      fprintf (stderr,
               "failed to parse URI: %s\n"
               "error message:       %s\n",
               uri_string,
               error.message);
      return EXIT_FAILURE;
   }
   /*
    * Create a new client instance
    */
   client = mongoc_client_new_from_uri (uri);
   if (!client) {
      return EXIT_FAILURE;
   }
   /*
    * Register the application name so we can track it in the profile logs
    * on the server. This can also be done from the URI (see other examples).
    */
   mongoc_client_set_appname (client, "lipidraft-mongodb-connect");
   /*
    * Get a handle on the database "db_name" and collection "coll_name"
    */
   database = mongoc_client_get_database (client, database_name);
   sprintf(namebuffer,"%d", collection_count);
   sprintf(collection_name,"%s",collection_name_head);
   strncat(collection_name,namebuffer,3 * sizeof(int)+1);
   printf("data collection %s created.\n",collection_name);
   collection = mongoc_client_get_collection (
       client,
       database_name,
       collection_name
       );
  
  return EXIT_SUCCESS;
}

void traverselattice(FuncPtr addition){
  for (x=0;x<LATTICE_SIZE;x++){
    for (y=0;y<LATTICE_SIZE;y++){

       //the key or colonm of datas
       keylen = bson_uint32_to_string (numpl++, &key, keybuf, sizeof keybuf);
       if (lattice.matrix[x][y].chol0 == C)
         keylen1 = bson_uint32_to_string (numchol++, &key1, keybuf1, sizeof keybuf);
       if (lattice.matrix[x][y].chol1 == C)
          keylen2 = bson_uint32_to_string (numchol++, &key2, keybuf2, sizeof keybuf);
       //xy
       plx = ((double)x) + ((double)y)* 0.5;
       ply = ((double)y) * cos30;

       // the func block
       addition();
       // clear keybuf
       memset(keybuf, 0, sizeof(keybuf));
       memset(keybuf1, 0, sizeof(keybuf1));
       memset(keybuf2, 0, sizeof(keybuf2));
     }
  }
}
///////////////////////////////////////////////
 void savepltype(){
   if (lattice.matrix[x][y].pl == Ahead)
      sprintf(pltype,"A");
  else
      sprintf(pltype,"B");

   BSON_APPEND_UTF8(grandchild1+4,key,pltype);
          };
 void savex(){ 
          BSON_APPEND_DOUBLE (grandchild1,key, plx);
          };
 void savey(){ 
          BSON_APPEND_DOUBLE (grandchild1+1,key, ply);
          };
 void savevarphi(){ 
        BSON_APPEND_DOUBLE (grandchild1+2,key,lattice.matrix[x][y].sc_pl.varphi);
        };
void savetheta(){
        BSON_APPEND_DOUBLE (grandchild1+3,key, lattice.matrix[x][y].sc_pl.theta);
        };
void savecholx(){
        if (lattice.matrix[x][y].chol0 == C){
          chol0x= plx + head_distance/2;
          BSON_APPEND_DOUBLE (grandchild2,key1, chol0x);
        }
       if (lattice.matrix[x][y].chol1 == C){
          chol1x= plx + head_distance/2;
          BSON_APPEND_DOUBLE (grandchild2,key2, chol1x);
       }
        };
void savecholy(){
        if (lattice.matrix[x][y].chol0 == C){
          chol0y= ply + head_distance/2*tan30;
          BSON_APPEND_DOUBLE (grandchild2+1,key1, chol0y);
        }
        if (lattice.matrix[x][y].chol1 == C){
          chol1y = ply - head_distance/2*tan30;
          BSON_APPEND_DOUBLE (grandchild2+1,key2, chol1y);
        }
        };
void savecholvarphi(){
       if (lattice.matrix[x][y].chol0 == C)
          BSON_APPEND_DOUBLE (grandchild2+2,key1, lattice.matrix[x][y].sc_chol0.varphi);
       if (lattice.matrix[x][y].chol1 == C)
          BSON_APPEND_DOUBLE (grandchild2+2,key2, lattice.matrix[x][y].sc_chol1.varphi); 
       };
void savecholtheta(){
       if (lattice.matrix[x][y].chol0 == C)
          BSON_APPEND_DOUBLE (grandchild2+3,key1, lattice.matrix[x][y].sc_chol0.theta);
        
       if (lattice.matrix[x][y].chol1 == C)
          BSON_APPEND_DOUBLE (grandchild2+3,key2, lattice.matrix[x][y].sc_chol1.theta);
        };
///////////////////////////////////////////////
void savelattice(enum STEP step){
  FuncPtr addition;

  datas = bson_new();
  switch (step){
  case Begin: sprintf(stepstr,"begin");break;
  case End:   sprintf(stepstr,"end");break;
  }
          // save
BSON_APPEND_DOCUMENT_BEGIN(datas, stepstr, &child0);
         //save energy
  BSON_APPEND_DOUBLE(&child0,"energy",currentenergy);
         //save pl
  BSON_APPEND_DOCUMENT_BEGIN(&child0,"phospholipids",&child1);
         //save pltype
        BSON_APPEND_ARRAY_BEGIN (&child1,"phospholipid type",grandchild1+4);
          addition = savepltype;
          traverselattice(addition);
        bson_append_array_end(&child1,grandchild1+4);

        BSON_APPEND_ARRAY_BEGIN (&child1,"x",grandchild1);
          addition = savex;
          traverselattice(addition);
        bson_append_array_end (&child1,grandchild1);

      BSON_APPEND_ARRAY_BEGIN (&child1,"y",grandchild1+1);
        addition = savey;
        traverselattice(addition);
      bson_append_array_end (&child1,grandchild1+1);

      BSON_APPEND_ARRAY_BEGIN (&child1,"varphi",grandchild1+2);
        addition = savevarphi;
        traverselattice(addition);
      bson_append_array_end (&child1,grandchild1+2);

      BSON_APPEND_ARRAY_BEGIN (&child1,"theta",grandchild1+3);
        addition = savetheta;
        traverselattice(addition);
      bson_append_array_end (&child1,grandchild1+3);
    
  bson_append_document_end(&child0, &child1);

          // save chol

  BSON_APPEND_DOCUMENT_BEGIN(&child0,"cholesterol",&child2);
        

      BSON_APPEND_ARRAY_BEGIN (&child2, "x",grandchild2);
        addition =savecholx;
        traverselattice(addition);
      bson_append_array_end (&child2,grandchild2);


      BSON_APPEND_ARRAY_BEGIN (&child2,"y",grandchild2+1);
        addition = savecholy;
        traverselattice(addition);
      bson_append_array_end (&child2,grandchild2+1);
      

      BSON_APPEND_ARRAY_BEGIN (&child2,"varphi",grandchild2+2);
       
       addition = savecholvarphi;
       traverselattice(addition);

      bson_append_array_end (&child2,grandchild2+2);


      BSON_APPEND_ARRAY_BEGIN (&child2,"theta",grandchild2+3);
        
        addition = savecholtheta;
        traverselattice(addition);

      bson_append_array_end (&child2,grandchild2+3);

  bson_append_document_end(&child0, &child2);
bson_append_document_end(datas, &child0);
  //insert
  if (!mongoc_collection_insert_one (
           collection, datas, NULL, NULL, &error)) {
        fprintf (stderr, "%s\n", error.message);
    }
  //
  bson_destroy(datas);
}

int exitdb(){
  /*
    * Release our handles and clean up libmongoc
    */
   mongoc_collection_destroy (collection);
   mongoc_database_destroy (database);
   mongoc_uri_destroy (uri);
   mongoc_client_destroy (client);
   mongoc_cleanup ();
   collection_count += 1;
   return EXIT_SUCCESS;
}
/* one exp */
void oneexp(
    double in_temp,
    double in_ca, double in_occupancy ,
    double in_alength, double in_blength
    ){
   initdb();
  /*
   * save parameters
   * */
   bson_oid_init (&oid, NULL);
   parameters = BCON_NEW(
      "parameters","{",
      //"_id",BCON_OID(&oid),
      "experiment id", BCON_INT32(collection_count),
      "temperature", BCON_DOUBLE(in_temp),
      "molar fraction of A", BCON_DOUBLE(in_ca),
      "ccupancy fraction of C", BCON_DOUBLE(in_occupancy),
      "tail chain length of A", BCON_DOUBLE(in_alength),
      "tail chain length of B", BCON_DOUBLE(in_blength),
      "}"
       ); 
   if (!mongoc_collection_insert_one (
           collection, parameters, NULL, NULL, &error)) {
        fprintf (stderr, "%s\n", error.message);
    }
   bson_destroy (parameters);
   /*
   * save datas
   * */
  init(in_temp,in_ca,in_occupancy,in_alength,in_blength);
  //printlattice();
  savelattice(Begin);
  stepuntilbalanced();
  printf("iteration: %d,endenergy:%lf\n",iteration,currentenergy);
  savelattice(End);
  //printlattice();
  iterationdata = bson_new();
  BSON_APPEND_INT32(iterationdata,"iteration",iteration);
  if (!mongoc_collection_insert_one (
           collection, iterationdata, NULL, NULL, &error)) {
        fprintf (stderr, "%s\n", error.message);
    }
  bson_destroy(iterationdata);
  exitdb();
}
