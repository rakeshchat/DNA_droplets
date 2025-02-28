//for spliting files to make gif or to make png files from different initial config1.
//this program creates a file with all the splitted file names.

//============================================
//to generate gif:
//convert -delay 4 -loop 1 image-*.xpm animated.gif
//(4 means 4/100 second delay between each frame, 1 means no loop)
//============================================

//to generate png files
//after executing splitfiles.c
//ececute the generated file dosplit.sh
//ls -1 *.xpm | xargs -n 1 bash -c 'convert "$0" "${0%.*}.png"'
//rm *.xpm
//now run the matlab peak code
//to print every 7th line of a file
//awk 'NR % 7 == 0' input > output
//gifsicle -i input.gif --optimize=3 -o output.gif
//convert input.gif -fuzz 10% -layers Optimize output.gif
//ffmpeg -i input.ogv -vcodec libx264 "output.mp4"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main()
{

long int Tm,time;
int      L,tmax,counter=0,count=0,s,step1,i,j,k;
double   ens;
char     c;
char     filename[2000];

FILE *fp,*fo;

    printf("Enter file name: "); 
    scanf("%s", filename);


    printf("Enter system size L= ");
    scanf("%d", &L);

    // Open the file
    fp = fopen(filename, "r"); 
  
  
    // Extract characters from file and store in character c 
    for (c = getc(fp); c != EOF; c = getc(fp)) 
        if (c == '\n') // Increment count if this character is newline 
            count = count + 1; 

//=============

//L= ;

//=============

step1 = 12+L;
ens  = count/step1;
Tm   = step1*ens;

fo = fopen("dosplit.sh","w");

//=============
counter=0;


for(time=0;time<Tm;time+=step1){
fprintf(fo,"awk \'NR>=%ld && NR<=%d\' %s > mat2imagel102-%06d.xpm\n",time+1, step1+step1*counter, filename, counter+1);
fprintf(fo,"convert mat2imagel102-%06d.xpm mat2imagel102-%06d.png\n",counter+1, counter+1);
counter++;
                              }


fclose(fo); fclose(fp);

//printf("-> %s has no. of lines -> %d, no. of snapshots => %d\n", filename, count, ens);
printf("-> %s has snapshots ==> %.1f\n", filename, ens);


system("chmod +x *.sh");
system("sleep 5");
system("./dosplit.sh");

//system("sleep 5"); 

//system("convert -delay 14 -loop 1 matimage-*.xpm c3.40.gif"); system("sleep 5"); 

//system("rm mat2image*.xpm");


}//main ends
