//Kirstin model, two species with a surfactant, before cleaving
//                 _____
//                |     |
//                |  3  |
//            ____|     |____
//           |               |
//           | 3     2     3 |        red type
//           |____       ____|
//                |     |
//                |  3  |
//                |_____|      

//                 _____
//                |     |
//                |  6  |
//            ____|     |____
//           |               |
//           | 6     5     6 |        green type
//           |____       ____|
//                |     |
//                |  6  |
//                |_____|    

//                 _____
//                |     |
//                | 4R  |
//            ____|     |____
//           |               |
//           | 4R    9    7G |        cross-links, yellow type
//           |____       ____|
//                |     |
//                | 7G  |
//                |_____|    






//model of crosses, a[row][column], max total density 0.137
//L(32),Tr(8) = 20 min
//L(40),Tr(8) = 40 min
//L(50),Tr(8) = 50 min
//L(64),Tr(8) = 1 hr 40 min
//L(100),Tr(8) = 5 hr 35 min
//L(128),Tr(6) = 6 min
//L(400),Tr(6) = 45 min


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          52  

static long long    a[L][L]={0},b[L][L]={0},t,Tr,Tm;


int main()
{clock_t tic = clock();

int    Nr,Ng,Ny,Nyp,y,z,i,j,l,p,q,pp,dd,cv,rhor1,rhog1,rhoy1,igg1,irr1,irg1,ioo1,temp1,count2=0,gap,datap,r=0,s=0,u=0,v=0;
double rhor,rhog,rhoy,c,d,e,f,Ei,Ef,dE,diff,temp,igg,irr,irg=0.0,ioo=1.0,rot,clev;
char   out1[1000],out2[1000];
FILE   *fp,*fq;
init_genrand(11111);


//parameters---------
pp    = 8;
dd    = 3;
rot   = 1;
//------------------- |irr| < |igg|
rhor  = 0.019;
rhog  = 0.019;
rhoy  = 0.019;

irr   =-2.8;
igg   = irr;
temp  = 0.65;
//---------------------


Tr    = pow(10,pp);
datap = pow(10,dd);
gap   = Tr/datap;


rhor1=rhor*1000; rhog1=rhog*1000; rhoy1=rhoy*1000; irr1=irr*10; igg1=igg*10; irg1=irg*10; ioo1=ioo*10; temp1=temp*100;

sprintf(out1,"beforeclvL%dt%dr%dg%dy%drr%dgg%dT%d.xpm",L,pp,rhor1,rhog1,rhoy1,irr1,igg1,temp1); fp=fopen(out1,"w");
sprintf(out2,"__beforeclvfnL%dt%dr%dg%dy%drr%dgg%dT%d",L,pp,rhor1,rhog1,rhoy1,irr1,igg1,temp1); fq=fopen(out2,"w");

//initialize to zero
for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0; b[y][z]=0;}   l=L;

//solid boundary ------ CHANGE L at TOP to (L+two) ------
for(y=0;y<L;y++){z=0;a[y][z]=1;} for(z=0;z<L;z++){y=0;a[y][z]=1;} for(y=0;y<L;y++){z=L-1;a[y][z]=1;} for(z=0;z<L;z++){y=L-1;a[y][z]=1;}       l=L-2;


Nr    = l*l*rhor;
Ng    = l*l*rhog;
Ny    = l*l*rhoy;


//random initialization of red,green,yellow particles
  i=0; while (i<Nr){
  y=genrand_real1()*L; z=genrand_real1()*L ;
  if(a[y][z]==0 && a[(y+1)%L][z]==0 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0 && a[y][(z-1+L)%L]==0){
  a[y][z]=2; a[(y+1)%L][z]=3; a[(y-1+L)%L][z]=3; a[y][(z+1)%L]=3; a[y][(z-1+L)%L]=3;
  i++ ; }
                   }

  i=0; while (i<Ng){
  y=genrand_real1()*L; z=genrand_real1()*L ;
  if(a[y][z]==0 && a[(y+1)%L][z]==0 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0 && a[y][(z-1+L)%L]==0){
  a[y][z]=5; a[(y+1)%L][z]=6; a[(y-1+L)%L][z]=6; a[y][(z+1)%L]=6; a[y][(z-1+L)%L]=6;
  i++ ; }
                   }
  i=0; while (i<Ny){
  y=genrand_real1()*L; z=genrand_real1()*L ;
  if(a[y][z]==0 && a[(y+1)%L][z]==0 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0 && a[y][(z-1+L)%L]==0){ c=genrand_real1();
            if(c < 0.25){
                       a[y][z]=9; 
                       a[(y+1)%L][z]=4; 
                       a[(y-1+L)%L][z]=7; 
                       a[y][(z+1)%L]=4; 
                       a[y][(z-1+L)%L]=7;}
            else if(c > 0.25 && c < 0.50){
                       a[y][z]=9; 
                       a[(y+1)%L][z]=7; 
                       a[(y-1+L)%L][z]=4; 
                       a[y][(z+1)%L]=7; 
                       a[y][(z-1+L)%L]=4;}
            else if(c > 0.50 && c < 0.75){
                       a[y][z]=9; 
                       a[(y+1)%L][z]=4; 
                       a[(y-1+L)%L][z]=7; 
                       a[y][(z+1)%L]=7; 
                       a[y][(z-1+L)%L]=4;}
            else if(c > 0.75){
                       a[y][z]=9; 
                       a[(y+1)%L][z]=7; 
                       a[(y-1+L)%L][z]=4;
                       a[y][(z+1)%L]=4; 
                       a[y][(z-1+L)%L]=7;}
            i++ ; 
                  }
                }



//======================================================================================================================
//relaxation dynamics starts
for(t=0;t<Tr;t++){


//Random symmetric rotation of linkers
for(p=0;p<L;p++){for(q=0;q<L;q++){
           if(a[p][q]==9){ 
              c=genrand_real1(); i=0;
                  if(c<rot){
                      if(genrand_real1()<0.5){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}


if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

//----------


if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-2+L)%L][q]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[p][(q-2+L)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[(p+2)%L][q]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+2)%L][q]==0){Ef+=ioo;}


if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[p][(q+2)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
                            u=a[p][(q+1)%L];
                            a[p][(q+1)%L]=a[(p-1+L)%L][q];
                            a[(p-1+L)%L][q]=a[p][(q-1+L)%L];
                            a[p][(q-1+L)%L]=a[(p+1)%L][q];
                            a[(p+1)%L][q]=u;
                                          }
                                }
          else {
             
//-----------------------

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}


if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

//----------


if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+2)%L][q]==0){Ef+=ioo;}


if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[p][(q-2+L)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-2+L)%L][q]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[p][(q+2)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==11){Ef+=irr;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[p][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
                            u=a[(p+1)%L][q];
                            a[(p+1)%L][q]=a[p][(q-1+L)%L];
                            a[p][(q-1+L)%L]=a[(p-1+L)%L][q];
                            a[(p-1+L)%L][q]=a[p][(q+1)%L];
                            a[p][(q+1)%L]=u;
                                           }
                                 }//if loop ends
                            }

                         }
               }}


//mcs starts
for(i=0;i<L;i++){for(j=0;j<L;j++){


p = genrand_real1()*L; 
q = genrand_real1()*L;
d = genrand_real1();


if(a[p][q]==5){
     //---------x+
    if(d < 0.25){
       if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 &&
          a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 &&
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==36){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}





if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[p][(q+1)%L]=5;
                                         a[p][(q+2)%L]=6;
                                         a[p][q]=6;
                                         a[(p-1+L)%L][(q+1)%L]=6;
                                         a[(p+1)%L][(q+1)%L]=6;
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;}
                                     }}

      //---------y+
   else if(d > 0.25 && d < 0.50){
        if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[(p-3+L)%L][q]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[(p-1+L)%L][(q-2+L)%L]!=5 &&
           a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2 &&
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][q]==0
){
        
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}





if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[(p-1+L)%L][q]=5;
                                         a[(p-2+L)%L][q]=6;
                                         a[p][q]=6;
                                         a[(p-1+L)%L][(q-1+L)%L]=6;
                                         a[(p-1+L)%L][(q+1)%L]=6;
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;}
                                     }}
     //---------x-
  else if(d > 0.50 && d < 0.75){
       if(a[(p-1+L)%L][(q-2+L)%L]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[p][(q-3+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 &&
          a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 &&
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[p][(q-2+L)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==36){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}





if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}
                
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[p][(q-1+L)%L]=5;
                                         a[p][(q-2+L)%L]=6;
                                         a[p][q]=6;
                                         a[(p-1+L)%L][(q-1+L)%L]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6;
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;}
                                    }}
     //---------y-
   else if(d > 0.75){
       if(a[(p+2)%L][(q+1)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+3)%L][q]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 &&
          a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 &&
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][q]==0
){
               
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==36){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}




if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}
 
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[(p+1)%L][q]=5;
                                         a[(p+2)%L][q]=6;
                                         a[p][q]=6;
                                         a[(p+1)%L][(q-1+L)%L]=6;
                                         a[(p+1)%L][(q+1)%L]=6;
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;}
                                     }}

                          }//green particle's movements end



else if(a[p][q]==2){
     //---------x+
    if(d < 0.25){
       if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 &&
          a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 &&
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}



if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[p][(q+1)%L]=2;
                                         a[p][(q+2)%L]=3;
                                         a[p][q]=3;
                                         a[(p-1+L)%L][(q+1)%L]=3;
                                         a[(p+1)%L][(q+1)%L]=3;
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;}
                                     }}

      //---------y+
   else if(d > 0.25 && d < 0.50){
        if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[(p-3+L)%L][q]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[(p-1+L)%L][(q-2+L)%L]!=5 &&
           a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2 &&
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][q]==0
){
        
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}




if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[(p-1+L)%L][q]=2;
                                         a[(p-2+L)%L][q]=3;
                                         a[p][q]=3;
                                         a[(p-1+L)%L][(q-1+L)%L]=3;
                                         a[(p-1+L)%L][(q+1)%L]=3;
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;}
                                     }}
     //---------x-
  else if(d > 0.50 && d < 0.75){
       if(a[(p-1+L)%L][(q-2+L)%L]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[p][(q-3+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 &&
          a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 &&
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[p][(q-2+L)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}



if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}
                
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[p][(q-1+L)%L]=2;
                                         a[p][(q-2+L)%L]=3;
                                         a[p][q]=3;
                                         a[(p-1+L)%L][(q-1+L)%L]=3;
                                         a[(p+1)%L][(q-1+L)%L]=3;
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;}
                                    }}
     //---------y-
   else if(d > 0.75){
       if(a[(p+2)%L][(q+1)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+3)%L][q]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 &&
          a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 &&
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][q]==0
){
               
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}



if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}
 
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[(p+1)%L][q]=2;
                                         a[(p+2)%L][q]=3;
                                         a[p][q]=3;
                                         a[(p+1)%L][(q-1+L)%L]=3;
                                         a[(p+1)%L][(q+1)%L]=3;
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;}
                                     }}

                          }//red particle's movements end


else if(a[p][q]==9){
     //---------x+
    if(d < 0.25){
       if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 &&
          a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 &&
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}





if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
                                         u=a[p][(q+1)%L];
                                         a[p][(q+1)%L]=9;
                                         a[p][(q+2)%L]=u;
                                         a[p][q]=a[p][(q-1+L)%L];
                                         a[(p-1+L)%L][(q+1)%L]=a[(p-1+L)%L][q];
                                         a[(p+1)%L][(q+1)%L]=a[(p+1)%L][q];
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;}
                                     }}

      //---------y+
   else if(d > 0.25 && d < 0.50){
        if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[(p-3+L)%L][q]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[(p-1+L)%L][(q-2+L)%L]!=5 &&
           a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2 &&
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][q]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}





if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
                                         u=a[(p-1+L)%L][q];
                                         a[(p-1+L)%L][q]=9;
                                         a[(p-2+L)%L][q]=u;
                                         a[p][q]=a[(p+1)%L][q];
                                         a[(p-1+L)%L][(q-1+L)%L]=a[p][(q-1+L)%L];
                                         a[(p-1+L)%L][(q+1)%L]=a[p][(q+1)%L];
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0;}
                                     }}
     //---------x-
  else if(d > 0.50 && d < 0.75){
       if(a[(p-1+L)%L][(q-2+L)%L]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[p][(q-3+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 &&
          a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 &&
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[p][(q-2+L)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}





if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
                                         u=a[p][(q-1+L)%L];
                                         a[p][(q-1+L)%L]=9;
                                         a[p][(q-2+L)%L]=u;
                                         a[p][q]=a[p][(q+1)%L];
                                         a[(p-1+L)%L][(q-1+L)%L]=a[(p-1+L)%L][q];
                                         a[(p+1)%L][(q-1+L)%L]=a[(p+1)%L][q];
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0;}
                                    }}
     //---------y-
   else if(d > 0.75){
       if(a[(p+2)%L][(q+1)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+3)%L][q]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 &&
          a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 &&
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][q]==0
){
               
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}





if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}

 
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
                                         u=a[(p+1)%L][q];
                                         a[(p+1)%L][q]=9;
                                         a[(p+2)%L][q]=u;
                                         a[p][q]=a[(p-1+L)%L][q];
                                         a[(p+1)%L][(q-1+L)%L]=a[p][(q-1+L)%L];
                                         a[(p+1)%L][(q+1)%L]=a[p][(q+1)%L];
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0;}
                                     }}

                          }//yellow particle's movements end


            }}//1 mcs ends


//for final off, for Gif on
//if(t==count2){



//for final on, for Gif off
   }//Trelax ends



//print------------------------------------------
fprintf(fp,"/* XPM */");            fprintf(fp,"\n");
fprintf(fp,"static char * MyArray[] = {");fprintf(fp,"\n");
fprintf(fp,"\"%d %d 10 1\",",(int)L,(int)L);   fprintf(fp,"\n");      //x,t,type object
fprintf(fp,"\"0      c #000000\",");fprintf(fp,"\n");    //empty sites grey 686969
fprintf(fp,"\"8      c #000000\",");fprintf(fp,"\n");    //blank
fprintf(fp,"\"6      c #1bcb0b\",");fprintf(fp,"\n");    //for green arms
fprintf(fp,"\"5      c #1bcb0b\",");fprintf(fp,"\n");    //for green center
fprintf(fp,"\"3      c #8ccdfb\",");fprintf(fp,"\n");    //for red arms
fprintf(fp,"\"2      c #8ccdfb\",");fprintf(fp,"\n");    //for red center
fprintf(fp,"\"9      c #fdee5d\",");fprintf(fp,"\n");    //for surf center yellow faf81c
fprintf(fp,"\"7      c #fdee5d\",");fprintf(fp,"\n");    //for surf arms-green end 
fprintf(fp,"\"4      c #fdee5d\",");fprintf(fp,"\n");    //for surf arm-red end yellow faf81c
fprintf(fp,"\"1      c #fd8e96\",");fprintf(fp,"\n");    //for boundary black

for(y=0;y<L;y++){fprintf(fp,"\""); for(z=0;z<L;z++){fprintf(fp,"%d",a[y][z]); } fprintf(fp,"\","); fprintf(fp,"\n");}
//print------------------------------------------


//for final off, for Gif on
//count2 = count2 + gap;} }//Trelax ends


//to print final configuration
for(i=0; i<L; i++) {
      for(j=0;j<L;j++) {
         fprintf(fq,"%d ", a[i][j]);
         if(j==L-1){
            fprintf(fq,"\n");
         }
      }
   } fclose(fq);


fclose(fp); 
clock_t toc = clock(); printf("#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
}//main ends
