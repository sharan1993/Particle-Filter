//Code Sample: Implementation of region segmentation on the basis of range.
//Conversion of range values to cartesian coordinates were already provided.
//Implemented a segmentation on the basis of calculated surface normals and similarity in
//the direction of these normals.
//By Sharan Rajendran
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int ROWS	=128;
int COLS	=128;
int MAX_QUEUE=128*128;

double SN[3][128*128];
/*
**	This routine converts the data in an Odetics range image into 3D
**	cartesian coordinate data.  The range image is 8-bit, and comes
**	already separated from the intensity image.
*/

//My code
//Function to write an image
void write_image(char filename[],int rows,int cols,unsigned char *image)
{
FILE *fpt;
fpt=fopen(filename,"w");
fprintf(fpt,"P5 %d %d 255\n",cols,rows);
fwrite(image,cols*rows,1,fpt);
fclose(fpt);
}

//Thresholding function
void threshold(unsigned char *in,int value,unsigned char **out)
{
(*out)=(unsigned char *)calloc(ROWS*COLS,sizeof(unsigned char*));
int i,j;
for(i=0;i<ROWS;i++)
for(j=0;j<COLS;j++)
if(in[i*COLS+j]>value)
(*out)[i*COLS+j]=255;
else
(*out)[i*COLS+j]=0;
write_image("thres.ppm",ROWS,COLS,(*out));
}

//Region grow function using queue
void RangeSegFill(unsigned char *t,unsigned char *image,int  ROWS,int COLS,int r,int c,int paint_value,int new_label,int *indices,int *count)
{
double average[]={0,0,0};
double total[]={0,0,0};
int r2,c2;

//queue initialization
int queue[128*128],qh,qt;
*count=0;
if (image[r*COLS+c] !=0)
  return;
image[r*COLS+c]=new_label;
if (indices != NULL)
  indices[0]=r*COLS+c;
queue[0]=r*COLS+c;
qh=1;	/* queue head */
qt=0;	/* queue tail */
*count=1;

double cosine,a,b;
//initialization of average vector components
average[0]=SN[0][r*COLS+c];
average[1]=SN[1][r*COLS+c];
average[2]=SN[2][r*COLS+c];

double dot;
while (qt != qh)
{
if(*count!=1)
{
	average[0]=total[0]/(*count);average[1]=total[1]/(*count);average[2]=total[2]/(*count);}
	for (r2=-1; r2<=1; r2++)
		for (c2=-1; c2<=1; c2++)
		{
			if (r2==0&&c2==0)
				continue;
			if ((queue[qt]/COLS+r2)<0||(queue[qt]/COLS+r2)>=ROWS||(queue[qt]%COLS+c2)<0||(queue[qt]%COLS+c2)>=COLS)
				continue;
			if (image[(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2] != 0)
				continue;
			a=sqrt(pow(average[0],2)+pow(average[1],2)+pow(average[2],2));
			average[0]=average[0]/a;average[1]=average[1]/a;average[2]=average[2]/a;
			b=sqrt(pow(SN[0][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2],2)+pow(SN[1][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2],2)+pow(SN[2][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2],2));
			cosine=(SN[0][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2]*average[0]+SN[1][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2]*average[1]+SN[2][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2]*average[2]);
	
			if (abs(cosine)<0.6)
				continue;
			image[(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2]=new_label;
			if (indices != NULL)
				indices[*count]=(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2;
			total[0]+=SN[0][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2];
			total[1]+=SN[1][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2];
			total[2]+=SN[2][(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2];
			(*count)++;
		        queue[qh]=(queue[qt]/COLS+r2)*COLS+queue[qt]%COLS+c2;
		        qh=(qh+1)%MAX_QUEUE;
  		        if (qh == qt)
        		{
				printf("Max queue size exceeded\n");
				write_image("max.ppm",ROWS,COLS,image);
			        exit(0);
		        }
	      }
	      qt=(qt+1)%MAX_QUEUE;
}
}
//End of my code

//Main function contains the conversion to cartesian coordinates
main(argc,argv)

int	argc;
char	*argv[];

{
int	r,c;
double	cp[7];
double	xangle,yangle,dist;
double	ScanDirectionFlag,SlantCorrection;
unsigned char	RangeImage[128*128],*thres;
double		P[3][128*128];
int             ImageTypeFlag;
char	Filename[160],Outfile[160];

FILE	*fpt;

//printf("Enter range image file name:");
//scanf("%s",Filename);

if ((fpt=fopen("chair-range.ppm","r")) == NULL)
  {
  printf("Couldn't open %s\n",Filename);
  return 0;
  }
fread(RangeImage,1,128*128,fpt);
fclose(fpt);

threshold(RangeImage,160,&thres);
printf("Up(-1), Down(1) or Neither(0)? ");
scanf("%d",&ImageTypeFlag);


cp[0]=1220.7;		/* horizontal mirror angular velocity in rpm */
cp[1]=32.0;		/* scan time per single pixel in microseconds */
cp[2]=(COLS/2)-0.5;		/* middle value of columns */
cp[3]=1220.7/192.0;	/* vertical mirror angular velocity in rpm */
cp[4]=6.14;		/* scan time (with retrace) per line in milliseconds */
cp[5]=(ROWS/2)-0.5;		/* middle value of rows */
cp[6]=10.0;		/* standoff distance in range units (3.66cm per r.u.) */

cp[0]=cp[0]*3.1415927/30.0;	/* convert rpm to rad/sec */
cp[3]=cp[3]*3.1415927/30.0;	/* convert rpm to rad/sec */
cp[0]=2.0*cp[0];		/* beam ang. vel. is twice mirror ang. vel. */
cp[3]=2.0*cp[3];		/* beam ang. vel. is twice mirror ang. vel. */
cp[1]/=1000000.0;		/* units are microseconds : 10^-6 */
cp[4]/=1000.0;			/* units are milliseconds : 10^-3 */

switch(ImageTypeFlag)
  {
  case 1:		/* Odetics image -- scan direction upward */
    ScanDirectionFlag=-1;
    break;
  case 0:		/* Odetics image -- scan direction downward */
    ScanDirectionFlag=1;
    break;
  default:		/* in case we want to do this on synthetic model */
    ScanDirectionFlag=0;
    break;
  }

	/* start with semi-spherical coordinates from laser-range-finder: */
	/*			(r,c,RangeImage[r*COLS+c])		  */
	/* convert those to axis-independant spherical coordinates:	  */
	/*			(xangle,yangle,dist)			  */
	/* then convert the spherical coordinates to cartesian:           */
	/*			(P => X[] Y[] Z[])			  */

if (ImageTypeFlag != 3)
  {
  for (r=0; r<ROWS; r++)
    {
    for (c=0; c<COLS; c++)
      {
      SlantCorrection=cp[3]*cp[1]*((double)c-cp[2]);
      xangle=cp[0]*cp[1]*((double)c-cp[2]);
      yangle=(cp[3]*cp[4]*(cp[5]-(double)r))+	/* Standard Transform Part */
	SlantCorrection*ScanDirectionFlag;	/*  + slant correction */
      dist=(double)RangeImage[r*COLS+c]+cp[6];
      P[2][r*COLS+c]=sqrt((dist*dist)/(1.0+(tan(xangle)*tan(xangle))
	+(tan(yangle)*tan(yangle))));
      P[0][r*COLS+c]=tan(xangle)*P[2][r*COLS+c];
      P[1][r*COLS+c]=tan(yangle)*P[2][r*COLS+c];
      }
    }
  }



//My Code

//P0-x cooridnate, p1-y coordinate, p2-z coordinate
double v1i,v2i,v1j,v2j,v1k,v2k;

//Distance between pixel of interest and pixels used for vector calculation is 5.
int vl=4;

for(r=0;r<ROWS;r++)
{
	for(c=0;c<COLS;c++)
	{
		v1i=P[0][(r+vl)*COLS+c]-P[0][r*COLS+c];
		v2i=P[0][(r)*COLS+c+vl]-P[0][r*COLS+c];
		v1j=P[1][(r+vl)*COLS+c]-P[1][r*COLS+c];
		v2j=P[1][(r)*COLS+c+vl]-P[1][r*COLS+c];
		v1k=P[2][(r+vl)*COLS+c]-P[2][r*COLS+c];
		v2k=P[2][(r)*COLS+c+vl]-P[2][r*COLS+c];
		
		//Calculation of Surface Normals
		SN[0][r*COLS+c]=v1j*v2k-v1k*v2j;
		SN[1][r*COLS+c]=v1k*v2i-v1i*v2k;
		SN[2][r*COLS+c]=v1i*v2j-v1j*v2i;
		SN[0][r*COLS+c]=SN[0][r*COLS+c]/sqrt(pow(SN[0][r*COLS+c],2)+pow(SN[1][r*COLS+c],2)+pow(SN[2][r*COLS+c],2));
		SN[1][r*COLS+c]=SN[1][r*COLS+c]/sqrt(pow(SN[0][r*COLS+c],2)+pow(SN[1][r*COLS+c],2)+pow(SN[2][r*COLS+c],2));
		SN[2][r*COLS+c]=SN[2][r*COLS+c]/sqrt(pow(SN[0][r*COLS+c],2)+pow(SN[1][r*COLS+c],2)+pow(SN[2][r*COLS+c],2));
	}
}

//Unsigned label image creation.
unsigned char *label;
label=(unsigned char*)calloc(ROWS*COLS,sizeof(unsigned char*));
for(r=0;r<ROWS;r++)
	for(c=0;c<COLS;c++)
		label[r*COLS+c]=0;

int RegionSize;
int region_no=1;

//unlabelled ->0
for(r=0;r<ROWS;r++)
{
	for(c=0;c<COLS;c++)
	{
		//Pixels in the threshold image that are white are assumed to form parts of the wall and is included in one region
		if(thres[r*COLS+c]==255)
			label[r*COLS+c]=1;
	}
}

//Creating other regions
int flag=0;int r1,c1;
//Ignores border pixels
for(r=2;r<ROWS-2;r++)
{
	for(c=2;c<COLS-2;c++)
	{
		flag=0;
		if(label[r*COLS+c]==0)
		{
			//Creates a 5x5 window to look for labelled pixels
			for(r1=-2;r1<=2;r1++)
			{
				for(c1=-2;c1<=2;c1++)
				{
					if(label[(r+r1)*COLS+c+c1]!=0)
						flag=1; //if labelled pixel is found, Region No is not incremented	
				}
			}
			if(flag==0)
			{
				region_no++; //if no labelled pixels are found, Region No is incremented.
				if(region_no==255) //if number of regions exceed 255, process in terminated with 255 regions
				{
					printf("Max label\n");
					break;
				}
			RangeSegFill(thres,label,ROWS,COLS,r,c,0,region_no,NULL,&RegionSize); //Call region grow function
			printf("Region labeled %d is %d in size\n",region_no,RegionSize);
			}
		}
	}
}

//Creates an image where regions have different grayscale colors
int color_label_inc=255/region_no;
for(r=0;r<ROWS;r++)
{
	for(c=0;c<COLS;c++)
	{
		label[r*COLS+c]=label[r*COLS+c]*color_label_inc;
	}
}

write_image("label.ppm",ROWS,COLS,(label));
return 0;
}



