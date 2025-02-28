//Kirstin model, two species with a surfactant with CLEAVAGE
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


//L(50),Tr(7) = 400 sec
//L(52),Tr(9) = 16 hr 30 min
//L(64),Tr(8) = 1 hr 36 min
//L(128),Tr(7) = 42 min
//L(256),Tr(6) = 45 min


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define       L    50  

static long long    a[L][L]={0},b[L][L]={0},t,Tm;


int main()
{clock_t tic = clock();

int    Nr,Ng,Ny,Nyp,y,z,i,j,l,p,q,pp,dd,e1,num1=0,ens,cv,rhog1,rhor1,rhoy1,igg1,irr1,irg1,ioo1,temp1,count2=0,gap,datap,r=0,s=0,u=0,v=0,posx=0,posy=0;
double rhor,rhog,rhoy,c,d,e,f,Ei,Ef,dE,diff,temp,igg,irr,ioo=1.0,irg=0.0,rot,clev;
char   out1[1000],out2[1000];
FILE   *fq2,*fq3,*fread_a,*fread_b;
//init_genrand((unsigned long)time(NULL));


//parameters---------
pp    = 7;
dd    = 3;
ens   = 1;
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


//************************************************************************
//           Dynamics AFTER CLEVING and with TRACER particle
//************************************************************************
for(e1=1;e1<=ens;e1++){


//--------------------------
//initialization
for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0; b[y][z]=0;}

num1=5; init_genrand(11000);

//change this for different tracer particles
fread_a= fopen("__5afterclvfnL50t8r15g15y15T70rr-28rg0","r");
for(y=0;y<L;y++){for(z=0;z<L;z++){fscanf(fread_a,"%lld",&a[y][z]);} }
//-------------------------------

sprintf(out2,"pathL%dt%dr%dg%dy%dT%drr%de_%03d",L,pp,rhor1,rhog1,rhoy1,temp1,irr1,num1); fq3=fopen(out2,"w");


fread_b= fopen("__1afterclvBnosL50t8r15g15y15T70rr-28rg0","r");
for(y=0;y<L;y++){for(z=0;z<L;z++){fscanf(fread_b,"%lld",&b[y][z]);} }

//Initial position of tracer particle
for(y=0;y<L;y++){for(z=0;z<L;z++){if(a[y][z]==9){posx=y; posy=z; fprintf(fq3,"%d  %d\n",posx,posy); }}}


//dynamics with cleaved state starts
for(t=0;t<Tm;t++){

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
               a[(p-1+L)%L][q]=0; a[(p+1)%L][q]=4; b[p][q]=4;
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
               a[p][(q-1+L)%L]=0; a[p][(q+1)%L]=4; b[p][q]=1;}
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
               a[(p+1)%L][q]=0; a[(p-1+L)%L][q]=4; b[p][q]=2;}
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
               a[p][(q+1)%L]=0; a[p][(q-1+L)%L]=4; b[p][q]=3;}
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
               a[p][(q+1)%L]=0; a[p][(q-1+L)%L]=4; b[p][q]=2;
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
               a[(p-1+L)%L][q]=0; a[(p+1)%L][q]=4; b[p][q]=3;}
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
               a[p][(q-1+L)%L]=0; a[p][(q+1)%L]=4; b[p][q]=4;}
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
               a[(p+1)%L][q]=0; a[(p-1+L)%L][q]=4; b[p][q]=1;}
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
                                         a[p][(q-1+L)%L]=0; posy = (posy+1)%L; fprintf(fq3,"%d  %d\n",posx,posy);}     
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
                                         a[p][(q-1+L)%L]=0; posx = (posx-1+L)%L; fprintf(fq3,"%d  %d\n",posx,posy);}
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
                                         a[p][(q+1)%L]=0; posy = (posy-1+L)%L; fprintf(fq3,"%d  %d\n",posx,posy);}
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
                                         a[p][(q+1)%L]=0; posx = (posx+1)%L; fprintf(fq3,"%d  %d\n",posx,posy);}
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


//for final off, for Gif on
//if(t==count2){


//for final on, for Gif off
   }//Tm ends


// //print------------------------------------------
// fprintf(fq2,"/* XPM */");            fprintf(fq2,"\n");
// fprintf(fq2,"static char * MyArray[] = {");fprintf(fq2,"\n");
// fprintf(fq2,"\"%d %d 10 1\",",(int)L,(int)L);   fprintf(fq2,"\n");      //x,t,type object
// fprintf(fq2,"\"0      c #000000\",");fprintf(fq2,"\n");    //empty sites
// fprintf(fq2,"\"6      c #1bcb0b\",");fprintf(fq2,"\n");    //for green arms
// fprintf(fq2,"\"5      c #1bcb0b\",");fprintf(fq2,"\n");    //for green center
// fprintf(fq2,"\"3      c #8ccdfb\",");fprintf(fq2,"\n");    //for blue arms
// fprintf(fq2,"\"2      c #8ccdfb\",");fprintf(fq2,"\n");    //for blue center
// fprintf(fq2,"\"9      c #fd6758\",");fprintf(fq2,"\n");    //for tracer particle
// fprintf(fq2,"\"8      c #fd6758\",");fprintf(fq2,"\n");    //for tracer particle
// fprintf(fq2,"\"7      c #1bcb0b\",");fprintf(fq2,"\n");    //for surf arms-green
// fprintf(fq2,"\"4      c #8ccdfb\",");fprintf(fq2,"\n");    //for surf arm-blue
// fprintf(fq2,"\"1      c #fd8e96\",");fprintf(fq2,"\n");    //for boundary red

// for(y=0;y<L;y++){fprintf(fq2,"\""); for(z=0;z<L;z++){fprintf(fq2,"%d",a[y][z]); } fprintf(fq2,"\","); fprintf(fq2,"\n");}
// //print------------------------------------------

//for final off, for Gif on
//count2 = count2 + gap;} }//Tm ends

//fclose(fq2); 


fclose(fq3);

     }//ens ends


clock_t toc = clock(); printf("#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
}//main ends
