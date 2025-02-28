
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{

long int Tm;
int      L,tmax,counter=0,count=0,s,step1,i,j,k,time;
double   ens;
char     c;
char     filename[2000];

FILE *fp,*fo,*ft;

    printf("Enter file name: "); 
    scanf("%s", filename); 


//------
ft = fopen("convert.sh","w");

for(time=1;time<=60;time++){
fprintf(ft,"convert %s-%d.xpm %s-%d.png\n",filename,time,filename,time);
                          }
fclose(ft);
system("chmod +x *.sh");
system("sleep 5");

//system("./convert.sh");
//system("rm convert.sh");


}//main ends
