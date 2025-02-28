//Kirstin model, two species with a surfactant with CLEAVAGE
//To simulate the average number of opposite arms as a nearest neighbour
//of a cross.
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

//                 _____
//                |     |
//                | 4R  |
//           _____|     |
//          |           |          cleaved red part of yellow
//          |  4R    4  | 
//          |___________|
//                 ___________
//                |           |
//                |  7    7G  |    cleaved green part of yellow
//                |      _____|
//                |     |
//                |  7G |
//                |_____|    


//model of crosses, a[row][column], max total density 0.16
//L(50,52),3*Tr(8),ens(1) = 7.5hrs(L50), (L52)


# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include "mt19937ar.c"


# define L          52  

static long long    a[L][L]={0},b[L][L]={0},t,Tr,Tm;


int main()
{clock_t tic = clock();

int    Nr,Ng,Ny,Nyp,N1,y,z,i,j,l,p,q,pp,qq,ee,dd,cv,rhog1,rhor1,rhoy1,igg1,irr1,irg1,ioo1,temp1,count2=0,gap,datap,r=0,s=0,u=0,v=0,c1,r1,r1a,c2,r2,r2a,e1,ens;
double rhor,rhog,rhoy,c,d,e,f,Ei,Ef,dE,diff,temp,igg,irr,irg=0.0,ioo,rot=1,clev=1,ds1,ds2;
char   out1[1000],out2[1000];
//init_genrand(101010);
init_genrand((unsigned long)time(NULL));
FILE   *fq,*fr,*fread;


//parameters---------
qq    = 3;
pp    = 8;
ens   = 3;

dd    = 5;  //data points
//-------------------
rhor  = 0.019;
rhog  = 0.019;
rhoy  = 0.019;

irr   =-2.8;
igg   =-2.9;
ioo   = 1.0;

temp  = 0.70;
//---------------------

Tr    = qq*pow(10,pp);
datap = pow(10,dd);
gap   = Tr/datap;

double *corr1; corr1 = (double*)calloc(Tr, sizeof(double)); for(t=0;t<Tr;t++){corr1[t]=0;}
double *corr2; corr2 = (double*)calloc(Tr, sizeof(double)); for(t=0;t<Tr;t++){corr2[t]=0;}

rhog1=rhog*1000; rhor1=rhor*1000; rhoy1=rhoy*1000; irr1=irr*10; igg1=igg*10; irg1=irg*10; ioo1=ioo*10; temp1=temp*100;


sprintf(out1,"finalcorrP2L%d_%dt%de%dr%dg%dy%drr%dgg%doo%dT%d",L,qq,pp,ens,rhor1,rhog1,rhoy1,irr1,igg1,ioo1,temp1); fq=fopen(out1,"w");


//ensemble average loop starts
for(e1=0;e1<ens;e1++){

//initialize to zero
for(y=0;y<L;y++)for(z=0;z<L;z++){a[y][z]=0; b[y][z]=0;}   l=L;

//boundary particles
for(y=0;y<L;y++){z=0;a[y][z]=8;} for(z=0;z<L;z++){y=0;a[y][z]=8;} for(y=0;y<L;y++){z=L-1;a[y][z]=8;} for(z=0;z<L;z++){y=L-1;a[y][z]=8;}       l=L-2;

//Reading from last configuration file
fread= fopen("__10fnL52t8r19g19y19rr-28gg-29oo10T65","r");
for(y=0;y<L;y++){for(z=0;z<L;z++){fscanf(fread,"%d",&a[y][z]);} }


Nr    = l*l*rhor;
Ng    = l*l*rhog;
Ny    = l*l*rhoy;
Nyp   = (Ny*2)/3;  //number of broken linkers to be added after cleaving



//Cleaving starts ================================
for(p=0;p<L;p++){for(q=0;q<L;q++){
           if(a[p][q]==9){
               if(genrand_real1()<clev){
                    if(genrand_real1()<0.5){r=0; s=0; u=0; v=0;
                           if(a[p][(q+1)%L]==7){a[p][(q+1)%L]=0;r=1;}
                           if(a[(p-1+L)%L][q]==7){a[(p-1+L)%L][q]=0;s=1;}
                           if(a[p][(q-1+L)%L]==7){a[p][(q-1+L)%L]=0;u=1;}
                           if(a[(p+1)%L][q]==7){a[(p+1)%L][q]=0;v=1;}
                       if(u==1 && v==1){b[p][q]=1; a[p][q]=4;}
                       else if(r==1 && v==1){b[p][q]=2; a[p][q]=4;}
                       else if(r==1 && s==1){b[p][q]=3; a[p][q]=4;}
                       else if(s==1 && u==1){b[p][q]=4; a[p][q]=4;}
                                            }
                     else               {r=0; s=0; u=0; v=0;
                           if(a[p][(q+1)%L]==4){a[p][(q+1)%L]=0;r=1;}
                           if(a[(p-1+L)%L][q]==4){a[(p-1+L)%L][q]=0;s=1;}
                           if(a[p][(q-1+L)%L]==4){a[p][(q-1+L)%L]=0;u=1;}
                           if(a[(p+1)%L][q]==4){a[(p+1)%L][q]=0;v=1;}
                       if(u==1 && v==1){b[p][q]=1; a[p][q]=7;}
                       else if(r==1 && v==1){b[p][q]=2; a[p][q]=7;}
                       else if(r==1 && s==1){b[p][q]=3; a[p][q]=7;}
                       else if(s==1 && u==1){b[p][q]=4; a[p][q]=7;}
                                         }
               	                      }
                         }
            }}


 i=0; while (i<Nyp){
  y=genrand_real1()*L; z=genrand_real1()*L ;
  if(a[y][z]==0){ c=genrand_real1();
            if(c<0.25 && a[(y-1+L)%L][z]==0 && a[y][(z+1)%L]==0){
            	 if(genrand_real1()<0.5){
                       a[y][z]=4; 
                       a[(y-1+L)%L][z]=4; 
                       a[y][(z+1)%L]=4; 
                       b[y][z]=1;
                       i++;
                                        }
                   else {
                   	   a[y][z]=7; 
                       a[(y-1+L)%L][z]=7; 
                       a[y][(z+1)%L]=7; 
                       b[y][z]=1;
                       i++;
                        }
                     }
             else if(c>0.25 && c<0.50 && a[(y-1+L)%L][z]==0 && a[y][(z-1+L)%L]==0){
            	 if(genrand_real1()<0.5){
                       a[y][z]=4; 
                       a[(y-1+L)%L][z]=4; 
                       a[y][(z-1+L)%L]=4; b[y][z]=2; i++;
                                        }
                   else {
                   	   a[y][z]=7; 
                       a[(y-1+L)%L][z]=7; 
                       a[y][(z-1+L)%L]=7; b[y][z]=2; i++;
                        }
                     }
              else if(c>0.50 && c<0.75 && a[(y+1)%L][z]==0 && a[y][(z-1+L)%L]==0){
            	 if(genrand_real1()<0.5){
                       a[y][z]=4; 
                       a[(y+1)%L][z]=4; 
                       a[y][(z-1+L)%L]=4; b[y][z]=3; i++;
                                        }
                   else {
                   	   a[y][z]=7; 
                       a[(y+1)%L][z]=7; 
                       a[y][(z-1+L)%L]=7; b[y][z]=3; i++;
                        }
                     }
              else if(c>0.75 && a[(y+1)%L][z]==0 && a[y][(z+1)%L]==0){
            	 if(genrand_real1()<0.5){
                       a[y][z]=4; 
                       a[(y+1)%L][z]=4; 
                       a[y][(z+1)%L]=4; b[y][z]=4; i++;
                                        }
                   else {
                   	   a[y][z]=7; 
                       a[(y+1)%L][z]=7; 
                       a[y][(z+1)%L]=7; b[y][z]=4; i++;
                        }
                     } 
                  }
             }
//=============================cleaving ends


//dynamics with cleaved state starts
for(t=0;t<Tr;t++){

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

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}



if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==42){Ef+=igg;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==42){Ef+=igg;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==42){Ef+=igg;}

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

if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==42){Ef+=igg;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==42){Ef+=igg;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==42){Ef+=igg;}


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


if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==24){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==36){Ei+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==42){Ei+=igg;}



if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==42){Ef+=igg;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==36){Ef+=igg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==42){Ef+=igg;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==42){Ef+=igg;}
       	        
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

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==24){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==36){Ei+=igg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==42){Ei+=igg;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==24){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==36){Ei+=igg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==42){Ei+=igg;}



if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==24){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==36){Ef+=igg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==42){Ef+=igg;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==24){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==36){Ef+=igg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==42){Ef+=igg;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==24){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==36){Ef+=igg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==42){Ef+=igg;}
 
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
          a[(p-2+L)%L][(q+1)%L]!=9 && a[(p-1+L)%L][(q+2)%L]!=9 && a[p][(q+3)%L]!=9 && a[(p+1)%L][(q+2)%L]!=9 && a[(p+2)%L][(q+1)%L]!=9
          && a[(p-1+L)%L][(q+1)%L]==0 && a[(p+1)%L][(q+1)%L]==0 && a[p][(q+2)%L]==0
){

//measuring interaction energy
Ei=0; Ef=0; dE=0; diff=0;

if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==21){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==18){Ei+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==12){Ei+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][q]==9){Ei+=irr;}

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}



if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q+1)%L]==9){Ef+=irr;}


if(a[p][(q+1)%L]*a[p][(q+3)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+3)%L]==9){Ef+=irr;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q+1)%L]==9){Ef+=irr;}

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

if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}


if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-3+L)%L][q]==9){Ef+=irr;}

if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p-1+L)%L][(q-2+L)%L]==9){Ef+=irr;}

if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p-1+L)%L][(q+2)%L]==9){Ef+=irr;}


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


if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}


if(a[(p+1)%L][q]*a[(p+2)%L][q]==21){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==18){Ei+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==12){Ei+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][q]==9){Ei+=irr;}



if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-3+L)%L]==9){Ef+=irr;}


if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p-1+L)%L][q]*a[(p-2+L)%L][(q-1+L)%L]==9){Ef+=irr;}


if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+2)%L][(q-1+L)%L]==9){Ef+=irr;}
                
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

if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==21){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==18){Ei+=irg;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==12){Ei+=irr;}
else if(a[p][(q-1+L)%L]*a[p][(q-2+L)%L]==9){Ei+=irr;}

if(a[p][(q+1)%L]*a[p][(q+2)%L]==21){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==18){Ei+=irg;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==12){Ei+=irr;}
else if(a[p][(q+1)%L]*a[p][(q+2)%L]==9){Ei+=irr;}



if(a[(p+1)%L][q]*a[(p+3)%L][q]==21){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==18){Ef+=irg;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==12){Ef+=irr;}
else if(a[(p+1)%L][q]*a[(p+3)%L][q]==9){Ef+=irr;}


if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==21){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==18){Ef+=irg;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==12){Ef+=irr;}
else if(a[p][(q-1+L)%L]*a[(p+1)%L][(q-2+L)%L]==9){Ef+=irr;}


if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==21){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==18){Ef+=irg;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==12){Ef+=irr;}
else if(a[p][(q+1)%L]*a[(p+1)%L][(q+2)%L]==9){Ef+=irr;}
 
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


//Calculation of Pearson correlation coefficient
ds1=0.0; ds2=0.0; N1=0;
c1 = 2; 
r1=3;  r1a=4;   //ds1
r2=6;  r2a=7;   //ds2


//if(c1==2){r1=6;} else if(c1==5){r1=3;}
//if(c1==2){r1=3;} else if(c1==5){r1=6;}

for(p=0;p<L;p++){for(q=0;q<L;q++){

if(a[p][q]==c1){
          if(a[p][(q+2)%L]==r1 || a[p][(q+2)%L]==r1a){ds1++;}
          if(a[p][(q-2+L)%L]==r1 || a[p][(q-2+L)%L]==r1a){ds1++;}
          if(a[(p-2+L)%L][q]==r1 || a[(p-2+L)%L][q]==r1a){ds1++;}
          if(a[(p+2)%L][q]==r1 || a[(p+2)%L][q]==r1a){ds1++;}
          if(a[(p-1+L)%L][(q+1)%L]==r1 || a[(p-1+L)%L][(q+1)%L]==r1a){ds1++;}
          if(a[(p-1+L)%L][(q-1+L)%L]==r1 || a[(p-1+L)%L][(q-1+L)%L]==r1a){ds1++;}
          if(a[(p+1)%L][(q+1)%L]==r1 || a[(p+1)%L][(q+1)%L]==r1a){ds1++;}
          if(a[(p+1)%L][(q-1+L)%L]==r1 || a[(p+1)%L][(q-1+L)%L]==r1a){ds1++;}
               

          if(a[p][(q+2)%L]==r2 || a[p][(q+2)%L]==r2a){ds2++;}
          if(a[p][(q-2+L)%L]==r2 || a[p][(q-2+L)%L]==r2a){ds2++;}
          if(a[(p-2+L)%L][q]==r2 || a[(p-2+L)%L][q]==r2a){ds2++;}
          if(a[(p+2)%L][q]==r2 || a[(p+2)%L][q]==r2a){ds2++;}
          if(a[(p-1+L)%L][(q+1)%L]==r2 || a[(p-1+L)%L][(q+1)%L]==r2a){ds2++;}
          if(a[(p-1+L)%L][(q-1+L)%L]==r2 || a[(p-1+L)%L][(q-1+L)%L]==r2a){ds2++;}
          if(a[(p+1)%L][(q+1)%L]==r2 || a[(p+1)%L][(q+1)%L]==r2a){ds2++;}
          if(a[(p+1)%L][(q-1+L)%L]==r2 || a[(p+1)%L][(q-1+L)%L]==r2a){ds2++;}
               }
         }} //corr loop ends

         corr1[t] += ds1;
         corr2[t] += ds2;

     }//Tmax ends

  }//ens ends


for(t=0;t<Tr;t+=gap){fprintf(fq,"%d    %lf   %lf\n",t,corr1[t]/ens,corr2[t]/ens);}


fclose(fq); free(corr1); free(corr2);
clock_t toc = clock(); printf("#CPU RUNTIME: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);
}//main ends
