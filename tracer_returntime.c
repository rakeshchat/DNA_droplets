//Kirstin model, two species with a surfactant with CLEAVAGE, to 
//model of crosses, a[row][column], max total density 0.16
//                 _____
//                |     |
//                |  3  |
//            ____|     |____
//           |               |
//           | 3     2     3 |        blue/red type
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
//           _____|     |
//          |           |          cleaved blue part of yellow
//          |  4R    4  | 
//          |___________|
//                 ___________
//                |           |
//                |  7    7G  |    cleaved green part of yellow
//                |      _____|
//                |     |
//                |  7G |
//                |_____|    


//                 _____
//                |     |
//                |  8  |
//            ____|     |____
//           |               |
//           |  8    9     8 |     tracer particle, blue/red
//           |____       ____|
//                |     |
//                |  8  |
//                |_____|   


//L(52),Tr(6) = 2.5 hrs


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define  L   52  

static long long    a[L][L]={0},b[L][L]={0},cls[L][L]={0},t,Tm;


int main()
{clock_t tic = clock();

int    Nr,Ng,Ny,Nyp,y,z,i,j,l,p,q,p1,q1,pp,dd,cv,rhog1,rhor1,rhoy1,igg1,irr1,irg1,ioo1,temp1,count2=0,gap,datap;
int    r=0,s=0,u=0,v=0,posx=0,posy=0,e1,ens,tim1=0,tim2=0,pt1=0,pt2=0,pt3=0;
double rhor,rhog,rhoy,c,d,e,f,Ei,Ef,dE,diff,temp,igg,irr,ioo=1.0,irg=0.0,rot=1.0,clev=1.0;
char   out1[1000],out2[1000],out3[1000],out4[1000];
FILE   *fb,*fq3,*fread_a,*fread_b,*fread_c,*fp2,*fq2;


//parameters---------
pp    = 6;
ens   = 2;
//-------------------
rhor  = 0.015;
rhog  = 0.015;
rhoy  = 0.015;

irr   =-2.80;
igg   = irr;
temp  = 0.70;
//---------------------
Tm    = pow(10,pp);
datap = pow(10,dd);
gap   = Tm/datap;


rhog1=rhog*1000; rhor1=rhor*1000; rhoy1=rhoy*1000; irr1=irr*10; igg1=igg*10; irg1=irg*10; ioo1=ioo*10; temp1=temp*100;


sprintf(out2,"return_time1"); fq3=fopen(out2,"w");

//************************************************************************
//           Dynamics AFTER CLEVING and with TRACER particle
//************************************************************************

for(e1=1;e1<=ens;e1++){

//init_genrand(278700);
init_genrand((unsigned long)time(NULL));

pt1=0; pt2=0; pt3=0;

//initialization
for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0; b[y][z]=0;}

fread_a= fopen("__aftercleaving52","r");
for(y=0;y<L;y++){for(z=0;z<L;z++){fscanf(fread_a,"%lld",&a[y][z]);} }

fread_b= fopen("__linkers52","r");
for(y=0;y<L;y++){for(z=0;z<L;z++){fscanf(fread_b,"%lld",&b[y][z]);} }

fread_c= fopen("__cluster52","r");
for(y=0;y<L;y++){for(z=0;z<L;z++){fscanf(fread_c,"%lld",&cls[y][z]);} }


//Initial position of tracer particle
//for(y=0;y<L;y++){for(z=0;z<L;z++){if(a[y][z]==9){posx=y; posy=z; fprintf(fq3,"%d  %d\n",posx,posy); }}}


//whole dynamics starts
for(t=0;t<Tm;t++){
//pt1=0; pt2=0;


//Random symmetric rotation of cleaved linkers
for(p=0;p<L;p++){for(q=0;q<L;q++){
  if(b[p][q]!=0){
  if(genrand_real1()<rot){ 
          if(genrand_real1()<0.5){ d=genrand_real1();

if(d<0.25 && a[(p-1+L)%L][q]==7 && a[p][(q+1)%L]==7 && b[p][q]==1 && a[(p+1)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
               a[(p-1+L)%L][q]=0; a[(p+1)%L][q]=7; b[p][q]=4;
                                           }
              }
          

else if(d>0.25 && d<0.50 && a[(p-1+L)%L][q]==7 && a[p][(q-1+L)%L]==7 && b[p][q]==2 && a[p][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q-1+L)%L]=0; a[p][(q+1)%L]=7; b[p][q]=1;}
                          }


else if(d>0.50 && d<0.75 && a[p][(q-1+L)%L]==7 && a[(p+1)%L][q]==7 && b[p][q]==3 && a[(p-1+L)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p+2)%L][q]*a[(p+1)%L][q]==49){Ei+=igg;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==42){Ei+=igg;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==28){Ei+=irg;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==21){Ei+=irg;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==0){Ei+=ioo;}

if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==49){Ef+=igg;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==42){Ef+=igg;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==28){Ef+=irg;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==21){Ef+=irg;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p+1)%L][q]=0; a[(p-1+L)%L][q]=7; b[p][q]=2;}
                             }
                          

else if(d>0.75 && a[p][(q+1)%L]==7 && a[(p+1)%L][q]==7 && b[p][q]==4 && a[p][(q-1+L)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q+1)%L]=0; a[p][(q-1+L)%L]=7; b[p][q]=3;}
                             }

//--------------

else if(d<0.25 && a[(p-1+L)%L][q]==4 && a[p][(q+1)%L]==4 && b[p][q]==1 && a[(p+1)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p-1+L)%L][q]=0; a[(p+1)%L][q]=4; b[p][q]=4; if(cls[p][(q+2)%L]==0 && cls[(p+2)%L][q]==0 && cls[(p+1)%L][(q+1)%L]==0 && cls[(p-1+L)%L][q]==0 && cls[(p-1+L)%L][(q+1)%L]==0){cls[p][q]=0; cls[(p+1)%L][q]=0; cls[p][(q+1)%L]=0;}
      else {cls[p][q]=1; cls[(p+1)%L][q]=1; cls[p][(q+1)%L]=1;}
                                           }
              }
          

else if(d>0.25 && d<0.50 && a[(p-1+L)%L][q]==4 && a[p][(q-1+L)%L]==4 && b[p][q]==2 && a[p][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q-1+L)%L]=0; a[p][(q+1)%L]=4; b[p][q]=1; if(cls[(p-2+L)%L][q]==0 && cls[p][(q+2)%L]==0 && cls[(p-1+L)%L][(q+1)%L]==0 && cls[(p-1+L)%L][(q-1+L)%L]==0 && cls[p][(q-1+L)%L]==0 && cls[(p+1)%L][q]==0 && cls[(p+1)%L][(q+1)%L]==0){cls[p][q]=0; cls[(p-1+L)%L][q]=0; cls[p][(q+1)%L]=0;}
       else {cls[p][q]=1; cls[(p-1+L)%L][q]=1; cls[p][(q+1)%L]=1;} }
                          }


else if(d>0.50 && d<0.75 && a[p][(q-1+L)%L]==4 && a[(p+1)%L][q]==4 && b[p][q]==3 && a[(p-1+L)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p+2)%L][q]*a[(p+1)%L][q]==12){Ei+=irr;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==16){Ei+=irr;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==24){Ei+=irg;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==28){Ei+=irg;}
else if(a[(p+2)%L][q]*a[(p+1)%L][q]==0){Ei+=ioo;}

if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==12){Ef+=irr;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==16){Ef+=irr;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==24){Ef+=irg;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==28){Ef+=irg;}
else if(a[(p-2+L)%L][q]*a[(p+1)%L][q]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p+1)%L][q]=0; a[(p-1+L)%L][q]=4; b[p][q]=2; if(cls[(p-2+L)%L][q]==0 && cls[p][(q-2+L)%L]==0 && cls[(p-1+L)%L][(q-1+L)%L]==0 && cls[(p-1+L)%L][(q+1)%L]==0 && cls[p][(q+1)%L]==0 && cls[(p+1)%L][(q-1+L)%L]==0 && cls[(p+1)%L][q]==0){cls[p][q]=0; cls[(p-1+L)%L][q]=0; cls[p][(q-1+L)%L]=0;}
      else {cls[p][q]=1; cls[(p-1+L)%L][q]=1; cls[p][(q-1+L)%L]=1;} }
                             }
                          

else if(d>0.75 && a[p][(q+1)%L]==4 && a[(p+1)%L][q]==4 && b[p][q]==4 && a[p][(q-1+L)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q+1)%L]=0; a[p][(q-1+L)%L]=4; b[p][q]=3; if(cls[p][(q-2+L)%L]==0 && cls[(p+2)%L][q]==0 && cls[(p+1)%L][(q-1+L)%L]==0 && cls[p][(q+1)%L]==0 && cls[(p+1)%L][(q+1)%L]==0 && cls[(p-1+L)%L][(q-1+L)%L]==0 && cls[(p-1+L)%L][q]==0){cls[p][q]=0; cls[(p+1)%L][q]=0; cls[p][(q-1+L)%L]=0;}
       else {cls[p][q]=1; cls[(p+1)%L][q]=1; cls[p][(q-1+L)%L]=1;} }
                             }

                   }//0.5 loop 1st

//---------------------------------------       

           else { d=genrand_real1();

if(d<0.25 && a[(p-1+L)%L][q]==7 && a[p][(q+1)%L]==7 && b[p][q]==1 && a[p][(q-1+L)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q+1)%L]=0; a[p][(q-1+L)%L]=7; b[p][q]=2;
                                           }
              }
          

else if(d>0.25 && d<0.50 && a[(p-1+L)%L][q]==7 && a[p][(q-1+L)%L]==7 && b[p][q]==2 && a[(p+1)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p-1+L)%L][q]=0; a[(p+1)%L][q]=7; b[p][q]=3;}
                          }


else if(d>0.50 && d<0.75 && a[p][(q-1+L)%L]==7 && a[(p+1)%L][q]==7 && b[p][q]==3 && a[p][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q-1+L)%L]=0; a[p][(q+1)%L]=7; b[p][q]=4;}
                             }
                          

else if(d>0.75 && a[p][(q+1)%L]==7 && a[(p+1)%L][q]==7 && b[p][q]==4 && a[(p-1+L)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p+1)%L][q]=0; a[(p-1+L)%L][q]=7; b[p][q]=1;}
                             }

//---------------

else if(d<0.25 && a[(p-1+L)%L][q]==4 && a[p][(q+1)%L]==4 && b[p][q]==1 && a[p][(q-1+L)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q-2+L)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q+1)%L]=0; a[p][(q-1+L)%L]=4; b[p][q]=2; if(cls[(p-2+L)%L][q]==0 && cls[p][(q-2+L)%L]==0 && cls[(p-1+L)%L][(q-1+L)%L]==0 && cls[(p-1+L)%L][(q+1)%L]==0 && cls[p][(q+1)%L]==0 && cls[(p+1)%L][(q-1+L)%L]==0 && cls[(p+1)%L][q]==0){cls[p][q]=0; cls[(p-1+L)%L][q]=0; cls[p][(q-1+L)%L]=0;}
      else {cls[p][q]=1; cls[(p-1+L)%L][q]=1; cls[p][(q-1+L)%L]=1;}
                                           }
              }
          

else if(d>0.25 && d<0.50 && a[(p-1+L)%L][q]==4 && a[p][(q-1+L)%L]==4 && b[p][q]==2 && a[(p+1)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p+2)%L][q]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p-1+L)%L][q]=0; a[(p+1)%L][q]=4; b[p][q]=3; if(cls[p][(q-2+L)%L]==0 && cls[(p+2)%L][q]==0 && cls[(p+1)%L][(q-1+L)%L]==0 && cls[p][(q+1)%L]==0 && cls[(p+1)%L][(q+1)%L]==0 && cls[(p-1+L)%L][(q-1+L)%L]==0 && cls[(p-1+L)%L][q]==0){cls[p][q]=0; cls[(p+1)%L][q]=0; cls[p][(q-1+L)%L]=0;}
       else {cls[p][q]=1; cls[(p+1)%L][q]=1; cls[p][(q-1+L)%L]=1;} }
                          }


else if(d>0.50 && d<0.75 && a[p][(q-1+L)%L]==4 && a[(p+1)%L][q]==4 && b[p][q]==3 && a[p][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][(q-1+L)%L]=0; a[p][(q+1)%L]=4; b[p][q]=4; if(cls[p][(q+2)%L]==0 && cls[(p+2)%L][q]==0 && cls[(p+1)%L][(q+1)%L]==0 && cls[(p-1+L)%L][q]==0 && cls[(p-1+L)%L][(q+1)%L]==0){cls[p][q]=0; cls[(p+1)%L][q]=0; cls[p][(q+1)%L]=0;}
      else {cls[p][q]=1; cls[(p+1)%L][q]=1; cls[p][(q+1)%L]=1;} }
                             }
                          

else if(d>0.75 && a[p][(q+1)%L]==4 && a[(p+1)%L][q]==4 && b[p][q]==4 && a[(p-1+L)%L][q]==0){

Ei=0; Ef=0; dE=0; diff=0;

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p-2+L)%L][q]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[(p+1)%L][q]=0; a[(p-1+L)%L][q]=4; b[p][q]=1; if(cls[(p-2+L)%L][q]==0 && cls[p][(q+2)%L]==0 && cls[(p-1+L)%L][(q+1)%L]==0 && cls[(p-1+L)%L][(q-1+L)%L]==0 && cls[p][(q-1+L)%L]==0 && cls[(p+1)%L][q]==0 && cls[(p+1)%L][(q+1)%L]==0){cls[p][q]=0; cls[(p-1+L)%L][q]=0; cls[p][(q+1)%L]=0;}
       else {cls[p][q]=1; cls[(p-1+L)%L][q]=1; cls[p][(q+1)%L]=1;} }
                             }

                      }//0.5 loop 2nd

                   }
               }
          }}//rotation ends

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
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9
          && a[(p-1+L)%L][(q+1)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==36){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==48){Ei+=irg;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==48){Ei+=irg;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==48){Ei+=irg;}



if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==48){Ef+=irg;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==48){Ef+=irg;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==48){Ef+=irg;}

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
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9
           && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p-1+L)%L][(q+1)%L]==0 && a[(p-2+L)%L][q]==0
){
        
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==48){Ei+=irg;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==48){Ei+=irg;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==48){Ei+=irg;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==48){Ef+=irg;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==48){Ef+=irg;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==48){Ef+=irg;}


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
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9
          && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0 && a[p][(q-2+L)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==36){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==48){Ei+=irg;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==48){Ei+=irg;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==48){Ei+=irg;}



if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==48){Ef+=irg;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==48){Ef+=irg;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==48){Ef+=irg;}
                
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
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9
          && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[(p+2)%L][q]==0
){
               
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==36){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==48){Ei+=irg;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==48){Ei+=irg;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==48){Ei+=irg;}



if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==48){Ef+=irg;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==48){Ef+=irg;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==48){Ef+=irg;}
 
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

//blue particle
else if(a[p][q]==2){
     //---------x+
    if(d < 0.25){
       if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 &&
          a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 &&
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9
          && a[(p-1+L)%L][(q+1)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irr;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irr;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irr;}



if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irr;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irr;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irr;}

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
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9
           && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p-1+L)%L][(q+1)%L]==0 && a[(p-2+L)%L][q]==0
){
        
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irr;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irr;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irr;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irr;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irr;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irr;}


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
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9
          && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0 && a[p][(q-2+L)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irr;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irr;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irr;}



if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irr;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irr;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irr;}
                
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
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9
          && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[(p+2)%L][q]==0
){
               
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irr;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irr;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irr;}



if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irr;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irr;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irr;}
 
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
          }//blue particle's movements end



//Tracer particle
else if(a[p][q]==9){
     //---------x+
    if(d < 0.25){
       if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[p][(q+3)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+2)%L][(q+1)%L]!=5 &&
          a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[p][(q+3)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+2)%L][(q+1)%L]!=2 &&
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9
          && a[(p-1+L)%L][(q+1)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==32){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==48){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==56){Ei+=irg;}


if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==32){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==48){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==56){Ei+=irg;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==32){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==48){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==56){Ei+=irg;}




if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==32){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==48){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==56){Ef+=irg;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==32){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==48){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==56){Ef+=irg;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==32){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==48){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==56){Ef+=irg;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[p][(q+1)%L]=9;
                                         a[p][(q+2)%L]=8;
                                         a[p][q]=8;
                                         a[(p-1+L)%L][(q+1)%L]=8;
                                         a[(p+1)%L][(q+1)%L]=8;
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0; posy = (posy+1)%L;
                                           }     
                                     }}
                                    

      //---------y+
   else if(d > 0.25 && d < 0.50){
        if(a[(p-2+L)%L][(q+1)%L]!=5 && a[(p-1+L)%L][(q+2)%L]!=5 && a[(p-3+L)%L][q]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[(p-1+L)%L][(q-2+L)%L]!=5 &&
           a[(p-2+L)%L][(q+1)%L]!=2 && a[(p-1+L)%L][(q+2)%L]!=2 && a[(p-3+L)%L][q]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[(p-1+L)%L][(q-2+L)%L]!=2 &&
           a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[(p-3+L)%L][q]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[(p-1+L)%L][(q-2+L)%L]!=9
           && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p-1+L)%L][(q+1)%L]==0 && a[(p-2+L)%L][q]==0
){
        
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==32){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==48){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==57){Ei+=irg;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==32){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==48){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==56){Ei+=irg;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==32){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==48){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==56){Ei+=irg;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==32){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==48){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==56){Ef+=irg;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==32){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==48){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==56){Ef+=irg;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==32){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==48){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==56){Ef+=irg;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[(p-1+L)%L][q]=9;
                                         a[(p-2+L)%L][q]=8;
                                         a[p][q]=8;
                                         a[(p-1+L)%L][(q-1+L)%L]=8;
                                         a[(p-1+L)%L][(q+1)%L]=8;
                                         a[p][(q+1)%L]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q-1+L)%L]=0; posx = (posx-1+L)%L;
                                           }
                                     }}                                     

     //---------x-
  else if(d > 0.50 && d < 0.75){
       if(a[(p-1+L)%L][(q-2+L)%L]!=5 && a[(p-2+L)%L][(q-1+L)%L]!=5 && a[p][(q-3+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 &&
          a[(p-1+L)%L][(q-2+L)%L]!=2 && a[(p-2+L)%L][(q-1+L)%L]!=2 && a[p][(q-3+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 &&
          a[(p-1+L)%L][(q-2+L)%L]!=9 && a[(p-2+L)%L][(q-1+L)%L]!=9 && a[p][(q-3+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9 && a[(p+2)%L][(q-1+L)%L]!=9
          && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0 && a[p][(q-2+L)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==32){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==48){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==56){Ei+=irg;}


if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==32){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==48){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==56){Ei+=irg;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==32){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==48){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==56){Ei+=irg;}



if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==32){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==48){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==56){Ef+=irg;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==9){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==32){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==48){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==56){Ef+=irg;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==32){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==48){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==56){Ef+=irg;}
                
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[p][(q-1+L)%L]=9;
                                         a[p][(q-2+L)%L]=8;
                                         a[p][q]=8;
                                         a[(p-1+L)%L][(q-1+L)%L]=8;
                                         a[(p+1)%L][(q-1+L)%L]=8;
                                         a[(p-1+L)%L][q]=0;
                                         a[(p+1)%L][q]=0;
                                         a[p][(q+1)%L]=0; posy = (posy-1+L)%L;
                                          }
                                    }}
                                    

     //---------y-
   else if(d > 0.75){
       if(a[(p+2)%L][(q+1)%L]!=5 && a[(p+1)%L][(q+2)%L]!=5 && a[(p+3)%L][q]!=5 && a[(p+2)%L][(q-1+L)%L]!=5 && a[(p+1)%L][(q-2+L)%L]!=5 &&
          a[(p+2)%L][(q+1)%L]!=2 && a[(p+1)%L][(q+2)%L]!=2 && a[(p+3)%L][q]!=2 && a[(p+2)%L][(q-1+L)%L]!=2 && a[(p+1)%L][(q-2+L)%L]!=2 &&
          a[(p+2)%L][(q+1)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+3)%L][q]!=9 && a[(p+2)%L][(q-1+L)%L]!=9 && a[(p+1)%L][(q-2+L)%L]!=9
          && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[(p+2)%L][q]==0
){
               
//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==32){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==48){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==56){Ei+=irg;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==32){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==48){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==56){Ei+=irg;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==32){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==48){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==56){Ei+=irg;}



if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==9){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==32){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==48){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==56){Ef+=irg;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==9){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==32){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==48){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==56){Ef+=irg;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==9){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==32){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==48){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==56){Ef+=irg;}
 
               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1();
               if(dE<0 || (dE>0 && f<diff)){
                                         a[(p+1)%L][q]=9;
                                         a[(p+2)%L][q]=8;
                                         a[p][q]=8;
                                         a[(p+1)%L][(q-1+L)%L]=8;
                                         a[(p+1)%L][(q+1)%L]=8;
                                         a[p][(q-1+L)%L]=0;
                                         a[(p-1+L)%L][q]=0;
                                         a[p][(q+1)%L]=0; posx = (posx+1)%L;
                                           }
                                     }}
                                     

                          }//TRACER particle's movements end

//=========== cleaved particles movements
else if(b[p][q]!=0){

//====....====
if(a[(p-1+L)%L][q]==7 && a[p][(q+1)%L]==7 && b[p][q]==1){c=genrand_real1();

if(c<0.25 && a[(p-1+L)%L][(q+1)%L]==0 && a[p][(q+2)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+3)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p-1+L)%L][q]=0;
               a[(p-1+L)%L][(q+1)%L]=7;
               a[p][(q+2)%L]=7;
               a[p][(q+1)%L]=7; b[p][q]=0; b[p][(q+1)%L]=1;}

                        }


else if(c>0.25 && c<0.50 && a[(p-2+L)%L][q]==0 && a[(p-1+L)%L][(q+1)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[p][(q+1)%L]=0;
               a[(p-1+L)%L][(q+1)%L]=7;
               a[(p-2+L)%L][q]=7;
               a[(p-1+L)%L][q]=7; b[p][q]=0; b[(p-1+L)%L][q]=1;}

                              }

else if(c>0.50 && c<0.75 && a[(p-1+L)%L][(q-1+L)%L]==0 && a[p][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[p][(q+1)%L]=0;
               a[(p-1+L)%L][q]=0;
               a[p][(q-1+L)%L]=7;
               a[(p-1+L)%L][(q-1+L)%L]=7; b[p][q]=0; b[p][(q-1+L)%L]=1; }

                              }

else if(c>0.75 && a[(p+1)%L][q]==0 && a[(p+1)%L][(q+1)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[p][(q+1)%L]=0;
               a[(p-1+L)%L][q]=0;
               a[(p+1)%L][q]=7;
               a[(p+1)%L][(q+1)%L]=7; b[p][q]=0; b[(p+1)%L][q]=1;}

                              }

                  }//1/4th part of 977 ends
//====....====


//====....====
else if(a[(p-1+L)%L][q]==7 && a[p][(q-1+L)%L]==7 && b[p][q]==2){c=genrand_real1();

if(c<0.25 && a[(p-1+L)%L][(q+1)%L]==0 && a[p][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[(p-1+L)%L][q]=0;
               a[(p-1+L)%L][(q+1)%L]=7;
               a[p][(q+1)%L]=7;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[p][(q+1)%L]=2;}

                        }


else if(c>0.25 && c<0.50 && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p-2+L)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[p][(q-1+L)%L]=0;
               a[(p-2+L)%L][q]=7;
               a[(p-1+L)%L][(q-1+L)%L]=7;
               a[(p-1+L)%L][q]=7; b[p][q]=0; b[(p-1+L)%L][q]=2;}

                              }

else if(c>0.50 && c<0.75 && a[(p-1+L)%L][(q-1+L)%L]==0 && a[p][(q-2+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p-1+L)%L][q]=0;
               a[p][(q-2+L)%L]=7;
               a[(p-1+L)%L][(q-1+L)%L]=7;
               a[p][(q-1+L)%L]=7; b[p][q]=0; b[p][(q-1+L)%L]=2;}

                              }

else if(c>0.75 && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+1)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==49){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==42){Ei+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}
Ef+=ioo;

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[(p-1+L)%L][q]=0;
               a[(p+1)%L][q]=7;
               a[(p+1)%L][(q-1+L)%L]=7;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[(p+1)%L][q]=2;}

                              }

                  }//2/4th part of 977 ends
//====....====


//====....====
else if(a[p][(q-1+L)%L]==7 && a[(p+1)%L][q]==7 && b[p][q]==3){c=genrand_real1();

if(c<0.25 && a[p][(q+1)%L]==0 && a[(p+1)%L][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[(p+1)%L][q]=0;
               a[(p+1)%L][(q+1)%L]=7;
               a[p][(q+1)%L]=7;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[p][(q+1)%L]=3;}

                        }


else if(c>0.25 && c<0.50 && a[(p-1+L)%L][q]==0 && a[(p-1+L)%L][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[(p+1)%L][q]=0;
               a[(p-1+L)%L][q]=7;
               a[(p-1+L)%L][(q-1+L)%L]=7;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[(p-1+L)%L][q]=3;}

                              }

else if(c>0.50 && c<0.75 && a[p][(q-2+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=0;
               a[p][(q-1+L)%L]=7;
               a[(p+1)%L][(q-1+L)%L]=7;
               a[p][(q-2+L)%L]=7; b[p][q]=0; b[p][(q-1+L)%L]=3;}

                              }

else if(c>0.75 && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+2)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==49){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p+1)%L][q]*a[(p+3)%L][q]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==49){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==42){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=7;
               a[p][(q-1+L)%L]=0;
               a[(p+1)%L][(q-1+L)%L]=7;
               a[(p+2)%L][q]=7; b[p][q]=0; b[(p+1)%L][q]=3;}

                              }

                  }//3/4th part of 977 ends
//====....====


//====....====
else if(a[p][(q+1)%L]==7 && a[(p+1)%L][q]==7 && b[p][q]==4){c=genrand_real1();

if(c<0.25 && a[p][(q+2)%L]==0 && a[(p+1)%L][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=0;
               a[p][(q+1)%L]=7;
               a[p][(q+2)%L]=7;
               a[(p+1)%L][(q+1)%L]=7; b[p][q]=0; b[p][(q+1)%L]=4;}

                        }

else if(c>0.25 && c<0.50 && a[(p-1+L)%L][q]==0 && a[(p-1+L)%L][(q+1)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[(p+1)%L][q]=0;
               a[(p-1+L)%L][q]=7;
               a[(p-1+L)%L][(q+1)%L]=7;
               a[p][(q+1)%L]=0; b[p][q]=0; b[(p-1+L)%L][q]=4;}

                              }

else if(c>0.50 && c<0.75 && a[p][(q-1+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==49){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=7;
               a[(p+1)%L][q]=0;
               a[p][(q-1+L)%L]=7;
               a[(p+1)%L][(q-1+L)%L]=7;
               a[p][(q+1)%L]=0; b[p][q]=0; b[p][(q-1+L)%L]=4;}

                              }

else if(c>0.75 && a[(p+1)%L][(q+1)%L]==0 && a[(p+2)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==49){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p+1)%L][q]*a[(p+3)%L][q]==49){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==42){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==49){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==42){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=7;
               a[p][(q+1)%L]=0;
               a[(p+1)%L][(q+1)%L]=7;
               a[(p+2)%L][q]=7; b[p][q]=0; b[(p+1)%L][q]=4;}

                              }

                  }//4/4th part of 977 ends
//====....====




//====....====
if(a[(p-1+L)%L][q]==4 && a[p][(q+1)%L]==4 && b[p][q]==1){c=genrand_real1();

if(c<0.25 && a[(p-1+L)%L][(q+1)%L]==0 && a[p][(q+2)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p-1+L)%L][q]=0;
               a[(p-1+L)%L][(q+1)%L]=4;
               a[p][(q+2)%L]=4;
               a[p][(q+1)%L]=4; b[p][q]=0; b[p][(q+1)%L]=1;}

                        }


else if(c>0.25 && c<0.50 && a[(p-2+L)%L][q]==0 && a[(p-1+L)%L][(q+1)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[p][(q+1)%L]=0;
               a[(p-1+L)%L][(q+1)%L]=4;
               a[(p-2+L)%L][q]=4;
               a[(p-1+L)%L][q]=4; b[p][q]=0; b[(p-1+L)%L][q]=1;}

                              }

else if(c>0.50 && c<0.75 && a[(p-1+L)%L][(q-1+L)%L]==0 && a[p][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[p][(q+1)%L]=0;
               a[(p-1+L)%L][q]=0;
               a[p][(q-1+L)%L]=4;
               a[(p-1+L)%L][(q-1+L)%L]=4; b[p][q]=0; b[p][(q-1+L)%L]=1;}

                              }

else if(c>0.75 && a[(p+1)%L][q]==0 && a[(p+1)%L][(q+1)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[p][(q+1)%L]=0;
               a[(p-1+L)%L][q]=0;
               a[(p+1)%L][q]=4;
               a[(p+1)%L][(q+1)%L]=4; b[p][q]=0; b[(p+1)%L][q]=1;}

                              }

                  }//1/4th part of 944 ends
//====....====


//====....====
else if(a[(p-1+L)%L][q]==4 && a[p][(q-1+L)%L]==4 && b[p][q]==2){c=genrand_real1();

if(c<0.25 && a[(p-1+L)%L][(q+1)%L]==0 && a[p][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[(p-1+L)%L][q]=0;
               a[(p-1+L)%L][(q+1)%L]=4;
               a[p][(q+1)%L]=4;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[p][(q+1)%L]=2;}

                        }


else if(c>0.25 && c<0.50 && a[(p-1+L)%L][(q-1+L)%L]==0 && a[(p-2+L)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==0){Ef+=ioo;}

if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][(q-2+L)%L]*a[p][(q-1+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[p][(q-1+L)%L]=0;
               a[(p-2+L)%L][q]=4;
               a[(p-1+L)%L][(q-1+L)%L]=4;
               a[(p-1+L)%L][q]=4; b[p][q]=0; b[(p-1+L)%L][q]=2;}

                              }

else if(c>0.50 && c<0.75 && a[(p-1+L)%L][(q-1+L)%L]==0 && a[p][(q-2+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p-1+L)%L][q]=0;
               a[p][(q-2+L)%L]=4;
               a[(p-1+L)%L][(q-1+L)%L]=4;
               a[p][(q-1+L)%L]=4; b[p][q]=0; b[p][(q-1+L)%L]=2;}

                              }

else if(c>0.75 && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+1)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==16){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==24){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==28){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==0){Ei+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}
Ef+=ioo;

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[(p-1+L)%L][q]=0;
               a[(p+1)%L][q]=4;
               a[(p+1)%L][(q-1+L)%L]=4;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[(p+1)%L][q]=2;}

                              }

                  }//2/4th part of 944 ends
//====....====


//====....====
else if(a[p][(q-1+L)%L]==4 && a[(p+1)%L][q]==4 && b[p][q]==3){c=genrand_real1();

if(c<0.25 && a[p][(q+1)%L]==0 && a[(p+1)%L][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[(p+1)%L][q]=0;
               a[(p+1)%L][(q+1)%L]=4;
               a[p][(q+1)%L]=4;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[p][(q+1)%L]=3;}

                        }


else if(c>0.25 && c<0.50 && a[(p-1+L)%L][q]==0 && a[(p-1+L)%L][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[(p+1)%L][q]=0;
               a[(p-1+L)%L][q]=4;
               a[(p-1+L)%L][(q-1+L)%L]=4;
               a[p][(q-1+L)%L]=0; b[p][q]=0; b[(p-1+L)%L][q]=3;}

                              }

else if(c>0.50 && c<0.75 && a[p][(q-2+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=0;
               a[p][(q-1+L)%L]=4;
               a[(p+1)%L][(q-1+L)%L]=4;
               a[p][(q-2+L)%L]=4; b[p][q]=0; b[p][(q-1+L)%L]=3;}

                              }

else if(c>0.75 && a[(p+1)%L][(q-1+L)%L]==0 && a[(p+2)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==16){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==28){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}

if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==16){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==28){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=4;
               a[p][(q-1+L)%L]=0;
               a[(p+1)%L][(q-1+L)%L]=4;
               a[(p+2)%L][q]=4; b[p][q]=0; b[(p+1)%L][q]=3;}

                              }

                  }//3/4th part of 944 ends
//====....====


//====....====
else if(a[p][(q+1)%L]==4 && a[(p+1)%L][q]==4 && b[p][q]==4){c=genrand_real1();

if(c<0.25 && a[p][(q+2)%L]==0 && a[(p+1)%L][(q+1)%L]==0){

Ei=0; Ef=0; dE=0; diff=0;

Ei+=ioo;
if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==0){Ef+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==0){Ef+=ioo;}

               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=0;
               a[p][(q+1)%L]=4;
               a[p][(q+2)%L]=4;
               a[(p+1)%L][(q+1)%L]=4; b[p][q]=0; b[p][(q+1)%L]=4;}

                        }

else if(c>0.25 && c<0.50 && a[(p-1+L)%L][q]==0 && a[(p-1+L)%L][(q+1)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[(p+1)%L][q]=0;
               a[(p-1+L)%L][q]=4;
               a[(p-1+L)%L][(q+1)%L]=4;
               a[p][(q+1)%L]=0; b[p][q]=0; b[(p-1+L)%L][q]=4;}

                              }

else if(c>0.50 && c<0.75 && a[p][(q-1+L)%L]==0 && a[(p+1)%L][(q-1+L)%L]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==16){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==28){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==0){Ei+=ioo;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==0){Ef+=ioo;}
Ef+=ioo;


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=4;
               a[(p+1)%L][q]=0;
               a[p][(q-1+L)%L]=4;
               a[(p+1)%L][(q-1+L)%L]=4;
               a[p][(q+1)%L]=0; b[p][q]=0; b[p][(q-1+L)%L]=4;}

                              }

else if(c>0.75 && a[(p+1)%L][(q+1)%L]==0 && a[(p+2)%L][q]==0){

   Ei=0; Ef=0; dE=0; diff=0;
if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==16){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==28){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==0){Ei+=ioo;}
Ei+=ioo;


if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==16){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==28){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==0){Ef+=ioo;}

if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==16){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==28){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==0){Ef+=ioo;}


               dE   = Ef-Ei;
               diff = exp(-dE/temp);
               f = genrand_real1(); u=0;
               if(dE<0 || (dE>0 && f<diff)){
               a[p][q]=0;
               a[(p+1)%L][q]=4;
               a[p][(q+1)%L]=0;
               a[(p+1)%L][(q+1)%L]=4;
               a[(p+2)%L][q]=4; b[p][q]=0; b[(p+1)%L][q]=4;}

                              }

                  }//4/4th part of 944 ends
//====....====

           }//977 and 944 ends

      }}//1 mcs ends
      
      
//taking account of cluster particles for blue and tracer crosses
for(p1=0;p1<L;p1++){for(q1=0;q1<L;q1++){
  if(a[p1][q1]==2){
     if(cls[p1][(q1+2)%L]==0 && cls[(p1-2+L)%L][q1]==0 && cls[p1][(q1-2+L)%L]==0 && cls[(p1+2)%L][q1]==0 && 
        cls[(p1-1+L)%L][(q1+1)%L]==0 && cls[(p1-1+L)%L][(q1-1+L)%L]==0 && cls[(p1+1)%L][(q1-1+L)%L]==0 && cls[(p1+1)%L][(q1+1)%L]==0){cls[p1][q1]=0; cls[(p1+1)%L][q1]=0; cls[p1][(q1+1)%L]=0; cls[p1][(q1-1+L)%L]=0; cls[(p1-1+L)%L][q1]=0;}
//     else {cls[p1][q1]=1; cls[(p1+1)%L][q1]=1; cls[p1][(q1+1)%L]=1; cls[p1][(q1-1+L)%L]=1; cls[(p1-1+L)%L][q1]=1;}
                  }
           }}
      

//taking account of cluster particles for linkers
for(p1=0;p1<L;p1++){for(q1=0;q1<L;q1++){
  if(b[p1][q1]==4 && a[p1][q1]==4){
      if(cls[p1][(q1+2)%L]==0 && cls[(p1+2)%L][q1]==0 && cls[(p1+1)%L][(q1+1)%L]==0 && cls[(p1-1+L)%L][q1]==0 && cls[(p1-1+L)%L][(q1+1)%L]==0){cls[p1][q1]=0; cls[(p1+1)%L][q1]=0; cls[p1][(q1+1)%L]=0;}
      else {cls[p1][q1]=1; cls[(p1+1)%L][q1]=1; cls[p1][(q1+1)%L]=1;}
                              }
  else if(b[p1][q1]==3 && a[p1][q1]==4){
       if(cls[p1][(q1-2+L)%L]==0 && cls[(p1+2)%L][q1]==0 && cls[(p1+1)%L][(q1-1+L)%L]==0 && cls[p1][(q1+1)%L]==0 && cls[(p1+1)%L][(q1+1)%L]==0 && cls[(p1-1+L)%L][(q1-1+L)%L]==0 && cls[(p1-1+L)%L][q1]==0){cls[p1][q1]=0; cls[(p1+1)%L][q1]=0; cls[p1][(q1-1+L)%L]=0;}
       else {cls[p1][q1]=1; cls[(p1+1)%L][q1]=1; cls[p1][(q1-1+L)%L]=1;}
                              }
  else if(b[p1][q1]==2 && a[p1][q1]==4){
      if(cls[(p1-2+L)%L][q1]==0 && cls[p1][(q1-2+L)%L]==0 && cls[(p1-1+L)%L][(q1-1+L)%L]==0 && cls[(p1-1+L)%L][(q1+1)%L]==0 && cls[p1][(q1+1)%L]==0 && cls[(p1+1)%L][(q1-1+L)%L]==0 && cls[(p1+1)%L][q1]==0){cls[p1][q1]=0; cls[(p1-1+L)%L][q1]=0; cls[p1][(q1-1+L)%L]=0;}
      else {cls[p1][q1]=1; cls[(p1-1+L)%L][q1]=1; cls[p1][(q1-1+L)%L]=1;}
                              }
  else if(b[p1][q1]==1 && a[p1][q1]==4){
       if(cls[(p1-2+L)%L][q1]==0 && cls[p1][(q1+2)%L]==0 && cls[(p1-1+L)%L][(q1+1)%L]==0 && cls[(p1-1+L)%L][(q1-1+L)%L]==0 && cls[p1][(q1-1+L)%L]==0 && cls[(p1+1)%L][q1]==0 && cls[(p1+1)%L][(q1+1)%L]==0){cls[p1][q1]=0; cls[(p1-1+L)%L][q1]=0; cls[p1][(q1+1)%L]=0;}
       else {cls[p1][q1]=1; cls[(p1-1+L)%L][q1]=1; cls[p1][(q1+1)%L]=1;}
                              }
                      }}
      
      
for(y=0;y<L;y++){for(z=0;z<L;z++){if(a[y][z]==9){posx=y; posy=z; }}}


     if(pt1==0 && cls[posx][(posy+2)%L]==0 && cls[(posx+2)%L][posy]==0 && cls[posx][(posy-2+L)%L]==0 && cls[(posx-2+L)%L][posy]==0 && cls[(posx-1+L)%L][(posy-1+L)%L]==0 && cls[(posx-1+L)%L][(posy+1)%L]==0 && cls[(posx+1)%L][(posy+1)%L]==0 && cls[(posx+1)%L][(posy-1+L)%L]==0 && pt3==0){pt1=1; pt3=1; tim1=t; /*printf("tim1=%d   pt1=%d   pt2=%d   pt3=%d\n",tim1,pt1,pt2,pt3);*/}
                                              
     else if(cls[posx][(posy+2)%L]==1 || cls[(posx+2)%L][posy]==1 || cls[posx][(posy-2+L)%L]==1 || cls[(posx-2+L)%L][posy]==1 || cls[(posx-1+L)%L][(posy-1+L)%L]==1 || cls[(posx-1+L)%L][(posy+1)%L]==1 || cls[(posx+1)%L][(posy+1)%L]==1 || cls[(posx+1)%L][(posy-1+L)%L]==1){if(pt1==1){tim2=t; pt1=0; pt3=0; fprintf(fq3,"%d\n",tim2-tim1);}}
   
              }//Tm ends
              
//for(y=0;y<L;y++){for(z=0;z<L;z++){fprintf(fq3,"%d",cls[y][z]); } fprintf(fq3,"\n");} fprintf(fq3,"%d    %d\n",posx,posy);

          }//ens ends

fclose(fq3);
clock_t toc = clock(); printf("#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
}//main ends
