#include <stdio.h>
#include <stdlib.h>
#include <math.h>
typedef struct{
	unsigned char red;
	unsigned char green;
	unsigned char blue;
}pixel;
typedef struct{
unsigned char valoare;	
float val_corelatie;
}tablou;
typedef struct{
int I,J;
unsigned char VALOARE;
float VAL_CORELATIE;
unsigned char adevar;
}detect;
void citireHEADER2(FILE *f,int *header)
{
	fread(&header[0],2,1,f);//signature
    fread(&header[1],4,1,f);//size
    fread(&header[2],2,1,f);//rezzerved1
    fread(&header[3],2,1,f);//rezerved2
    fread(&header[4],4,1,f);//offset
    fread(&header[5],4,1,f);//bitmapinfo
    fread(&header[6],4,1,f);// AICI E LATIMEA==WIDTH[header]
    fread(&header[7],4,1,f);//AICI E INTALTIMEA==HEIGHT[header]
    fread(&header[8],2,1,f);//planes
    fread(&header[9],2,1,f);//bitsPERpixel
    fread(&header[10],4,1,f);///compression
    fread(&header[11],4,1,f);//imageSIZE
    fread(&header[12],4,1,f);//horizontalREZOLUTION
    fread(&header[13],4,1,f);//verticalREZOLUTION
    fread(&header[14],4,1,f);//numCOLORS
    fread(&header[15],4,1,f);//importantCOLORS
}
void copiereHEADER2(FILE *g,int *header)
{
	fwrite(&header[0],2,1,g);//signature
    fwrite(&header[1],4,1,g);//size
    fwrite(&header[2],2,1,g);//rezzerved1
    fwrite(&header[3],2,1,g);//rezerved2
    fwrite(&header[4],4,1,g);//offset
    fwrite(&header[5],4,1,g);//bitmapinfo
    fwrite(&header[6],4,1,g);// AICI E LATIMEA==WIDTH[header]
    fwrite(&header[7],4,1,g);//AICI E INTALTIMEA==HEIGHT[header]
    fwrite(&header[8],2,1,g);//planes
    fwrite(&header[9],2,1,g);//bitsPERpixel
    fwrite(&header[10],4,1,g);///compression
    fwrite(&header[11],4,1,g);//imageSIZE
    fwrite(&header[12],4,1,g);//horizontalREZOLUTION
    fwrite(&header[13],4,1,g);//verticalREZOLUTION
    fwrite(&header[14],4,1,g);//numCOLORS
    fwrite(&header[15],4,1,g);//importantCOLORS
}
pixel* citire_gri_liniarizare(FILE *imag,int *header)
{
	citireHEADER2(imag,header);
	pixel *p=(pixel*)malloc(header[6]*header[7]*sizeof(pixel));
	unsigned int i;
	for(i=0;i<header[6]*header[7];i++)
	fread(&p[i],3,1,imag);

	for(i=0;i<header[6]*header[7];i++)	
	{	unsigned char x=0.299*p[i].red+0.587*p[i].green+0.114*p[i].blue;
		p[i].red=x;
		p[i].green=x;
		p[i].blue=x;
	}
	return p;
}
pixel* urcare(FILE *f,int *header)

{
	citireHEADER2(f,header);
	pixel *p,padding;
	p=(pixel*)malloc(sizeof(pixel)*header[6]*header[7]);
	int i,j;
	for(i=0;i<header[7];i++)
		for(j=0;j<header[6];j++)
		{	fread(&p[i*header[6]+j],3,1,f);
			if(header[6]%4!=0)
			fread(&padding,header[6]%4,1,f);
		}
	return p;
}
pixel* liniarizare_cifra(char nume[20])
{
	FILE *f;
	f=fopen(nume,"rb");
	pixel* x=(pixel*)malloc(11*15*sizeof(pixel));
	unsigned i,j;
    int *header1;
	header1=(int*)malloc(20*sizeof(int));
	fread(&header1[0],2,1,f);//signature
    fread(&header1[1],4,1,f);//size
    fread(&header1[2],2,1,f);//rezzerved1
    fread(&header1[3],2,1,f);//rezerved2
    fread(&header1[4],4,1,f);//offset
    fread(&header1[5],4,1,f);//bitmapinfo
    fread(&header1[6],4,1,f);// AICI E LATIMEA==WIDTH[header]
    fread(&header1[7],4,1,f);//AICI E INTALTIMEA==HEIGHT[header]
    fread(&header1[8],2,1,f);//planes
    fread(&header1[9],2,1,f);//bitsPERpixel
    fread(&header1[10],4,1,f);///compression
    fread(&header1[11],4,1,f);//imageSIZE
    fread(&header1[12],4,1,f);//horizontalREZOLUTION
    fread(&header1[13],4,1,f);//verticalREZOLUTION
    fread(&header1[14],4,1,f);//numCOLORS
    fread(&header1[15],4,1,f);//importantCOLORS
	
	
    for(i=0;i<15;i++)
      {
       for(j=0;j<11; j++)
        fread(&x[i*11+j],3,1,f);
		
		fseek(f,3,SEEK_CUR);			//padding=3
      }
	  	for(i=0;i<11*15;i++)	
	{	float Z=0.299*x[i].red+0.587*x[i].green+0.114*x[i].blue;
		if(Z<=255)
		{x[i].red=(unsigned char)Z;	x[i].green=(unsigned char)Z; x[i].blue=(unsigned char)Z;}
		else if(Z<0)
		{x[i].red=0;	x[i].green=0;		x[i].blue=0;}
		else
		{x[i].red=255;		x[i].green=255;		x[i].blue=255;}
	}
	fclose(f);
	return x;
}
float  corelatie(unsigned x, unsigned y,pixel *imag,pixel *sablon,int *header)
{
	float S1=0,S2=0,S3=0,S4=0,S5=0,S6=0,SUM=0;
	unsigned i,j;
	
	for(i=0;i<165;i++)
			S1+=(float)sablon[i].red;
	S1/=(float)165; 
		
	for(i=0;i<165;i++)
			S2+=((float)(sablon[i].red-S1))*((float)(sablon[i].red-S1));
	S2/=(float)164;
	S2=sqrt(S2);	
	
	for(i=x;i<x+15;i++)
		for(j=y;j<y+11;j++)
		 	S3+=(float)imag[header[6]*i+j].red;
	S3/=(float)165;

	for(i=x;i<x+15;i++)
		for(j=y;j<y+11;j++)
		S4+=((float)(imag[header[6]*i+j].red-S3))*((float)(imag[header[6]*i+j].red-S3));	
	S4/=(float)164;
	S4=sqrt(S4);
	
	for(i=x;i<x+15;i++)
		for(j=y;j<y+11;j++)
		{
			S5=(float)(imag[header[6]*i+j].red-S3);
			S6=(float)(sablon[(i-x)*11+(j-y)].red-S1);
			SUM+= ((float)(S5*S6))/((float)(S4*S2));
		}
		
	SUM/=(float)165; 
	return SUM;
	
}
void desenculoare(pixel *imag, unsigned i, unsigned j,pixel culoare_chenar)
{
	unsigned x,y;
	//desenare:
	for(x=0;x<11;x++)
		imag[i*500+j+x]=culoare_chenar;
	for(x=0;x<11;x++)
		imag[(i+15)*500+j+x]=culoare_chenar;
	for(x=0;x<15;x++)
		imag[(i+x)*500+j]=culoare_chenar;
	for(x=0;x<15;x++)
		imag[(i+x)*500+j+11]=culoare_chenar;
}
float template_matching(pixel *imag,unsigned i, unsigned j,float  ps,char numesablon[20],int *header)
{

	pixel *Sablon; 
	unsigned k,l;
	Sablon=liniarizare_cifra(numesablon);
		
	//inversare sablon:
	for(k=0;k<15/2;k++)
		for(l=0;l<11;l++)
		{pixel aux=Sablon[k*11+l];
		Sablon[k*11+l]=Sablon[(15-k-1)*11+l];
		Sablon[(15-k-1)*11+l]=aux;}
		
	float  corr=corelatie(i,j,imag,Sablon,header);
	if(corr>=ps)
	return corr;
	return 0;
	
}
int comparare(const void *a, const void *b)
{
	tablou x=*(tablou*)a;
	tablou y=*(tablou*)b;
	if(x.val_corelatie>y.val_corelatie)
		return -1;
	return 1;
}
void eliminare_non_maxime_si_desenare(unsigned ABC,detect *detectie,pixel *imagine)
{	unsigned i,j;
	
	for(i=0;i<ABC-1;i++)
		for(j=i+1;j<ABC;j++)
			if((detectie[i].adevar+detectie[j].adevar)==2)
			{
				float intersectiei=0,intersectiej=0;
				intersectiei=(int)(detectie[i].I-detectie[j].I);
				if(intersectiei<0)
					intersectiei=intersectiei*(-1);
				intersectiej=(int)(detectie[i].J-detectie[j].J);
				if(intersectiej<0)
					intersectiej=intersectiej*(-1);
				
				if(intersectiei<15&&intersectiej<11)
				{
					float intersectie=(15-intersectiei)*(11-intersectiej);
					float suprapunere;
					suprapunere=(float)(intersectie/(165+165-intersectie));
					if(suprapunere>0.2)
					{
						if(detectie[i].VAL_CORELATIE>detectie[j].VAL_CORELATIE)
							detectie[j].adevar=0;
						else detectie[i].adevar=0;
					}
				}
			}
	for(i=0;i<ABC;i++)
		if(detectie[i].adevar==1)
		{		pixel culoare_chenar;
				switch(detectie[i].VALOARE){
						case(0):culoare_chenar.red=0; culoare_chenar.green=0; culoare_chenar.blue=255; break;
						case(1):culoare_chenar.red=0; culoare_chenar.green=255; culoare_chenar.blue=255; break;
						case(2):culoare_chenar.red=0; culoare_chenar.green=255; culoare_chenar.blue=0; break;
						case(3):culoare_chenar.red=255; culoare_chenar.green=255; culoare_chenar.blue=0; break;
						case(4):culoare_chenar.red=255; culoare_chenar.green=0; culoare_chenar.blue=255; break;
						case(5):culoare_chenar.red=255; culoare_chenar.green=0; culoare_chenar.blue=0; break;
						case(6):culoare_chenar.red=192; culoare_chenar.green=192; culoare_chenar.blue=192; break;
						case(7):culoare_chenar.red=0; culoare_chenar.green=140; culoare_chenar.blue=255; break;
						case(8):culoare_chenar.red=128; culoare_chenar.green=0; culoare_chenar.blue=128; break;
						case(9):culoare_chenar.red=0; culoare_chenar.green=0; culoare_chenar.blue=128; break; 	}		
				desenculoare(imagine,detectie[i].I,detectie[i].J,culoare_chenar);
			
		}
	
	
}
int main()
{

	printf("File name:");
	char nume_fisier[100];
	scanf("%s",nume_fisier);


	FILE *fout,*IMAGINE;
	IMAGINE=fopen(nume_fisier,"rb");
	if(IMAGINE==NULL)
		{printf("File not found\n");
return 0;}
else 
	printf("Wait...\n");
	fout=fopen("smart_file.bmp","wb");
	unsigned i,j,k,l,ABC;
	int *header=(int*)malloc(20*sizeof(int));
	char numesablon[20]="cifrax.bmp";
	unsigned cif;
	
	pixel *imag=citire_gri_liniarizare(IMAGINE,header);
	fclose(IMAGINE);
	IMAGINE=fopen(nume_fisier,"rb");
	pixel *imagine=urcare(IMAGINE,header);
	fclose(IMAGINE);
	detect *detectie=(detect*)malloc(header[6]*header[7]*sizeof(detect));
	
	for(ABC=0;ABC<header[6]*header[7];ABC++)
		detectie[ABC].adevar=1;
		//inversare imag:
	for(k=0;k<header[7]/2;k++)
		for(l=0;l<header[6];l++)
		{pixel aux=imag[k*header[6]+l]; imag[k*header[6]+l]=imag[(header[7]-k-1)*header[6]+l];imag[(header[7]-k-1)*header[6]+l]=aux;}
		//inversare imagine:
	for(k=0;k<header[7]/2;k++)
		for(l=0;l<header[6];l++)
		{pixel aux=imagine[k*header[6]+l];	imagine[k*header[6]+l]=imagine[(header[7]-k-1)*header[6]+l];imagine[(header[7]-k-1)*header[6]+l]=aux;}
	

	ABC=0;
	
		for(i=0;i<header[7]-15;i++)
			for(j=0;j<header[6]-11;j++)	
				for(cif=0;cif<=9;cif++)	
			{	
				tablou *v=(tablou*)malloc(10*sizeof(tablou));
				numesablon[5]=(char)(cif+'0');
				v[cif].val_corelatie=template_matching(imag,i,j,0.5,numesablon,header);
				v[cif].valoare=cif;
				qsort(v,10,sizeof(tablou),comparare);
				
				if(v[0].val_corelatie>=0.5){
				detectie[ABC].I=i;
				detectie[ABC].J=j;
				detectie[ABC].VALOARE=v[0].valoare;
				detectie[ABC].VAL_CORELATIE=v[0].val_corelatie;
				ABC++;}
				
				free(v);
				}
	
	eliminare_non_maxime_si_desenare(ABC,detectie,imagine);
	
	
//reinversare imagine:
	for(k=0;k<header[7]/2;k++)
		for(l=0;l<header[6];l++)
		{pixel aux=imagine[k*header[6]+l];
		imagine[k*header[6]+l]=imagine[(header[7]-k-1)*header[6]+l];
		imagine[(header[7]-k-1)*header[6]+l]=aux;}
		
		//afisare:
		copiereHEADER2(fout,header);
		for(i=0;i<header[6]*header[7];i++)
			fwrite(&imagine[i],3,1,fout);
	free(header);
	free(imag);
	free(imagine);    
	free(detectie);
	fclose(fout);
	fclose(IMAGINE);
	return 0;
}