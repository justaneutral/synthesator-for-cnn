#include "factorizator.h"

int *malloc1i(unsigned int size)
{
	int *arrayPtr=(int*)malloc(sizeof *arrayPtr*size);
	if(arrayPtr)
	{
		unsigned int i;
		for(i=0;i<size;i++)
		{
			arrayPtr[i]=0;
		}
	}
	return arrayPtr;
}

double *malloc1d(unsigned int size)
{
	double *arrayPtr=(double*)malloc(sizeof *arrayPtr*size);
	if(arrayPtr)
	{
		unsigned int i;
		for(i=0;i<size;i++)
		{
			arrayPtr[i]=0.0;
		}
	}
	return arrayPtr;
}

int **malloc2i(unsigned int rows,unsigned int cols)
{
	int **arrayPtr=(int**)malloc(sizeof *arrayPtr*rows);
	if(arrayPtr)
	{
		unsigned int i;
		for(i=0;i<rows;i++)
		{
			arrayPtr[i]=(int*)malloc(sizeof *arrayPtr[i]*cols);
			if(arrayPtr[i])
			{
				unsigned int j;
				for(j=0;j<cols;j++)
				{
					arrayPtr[i][j]=0;
				}
			}
		}
	}
	return arrayPtr;
}

void free2i(int **arrayPtr,unsigned rows)
{
	if(arrayPtr)
	{
		unsigned i;
		for(i=0;i<rows;i++)
		{
			free(arrayPtr[i]);
		}
	}
	free(arrayPtr);
	arrayPtr = NULL;
}

double **malloc2d(unsigned rows,unsigned cols)
{
	double **arrayPtr=(double**)malloc(sizeof *arrayPtr*rows);
	if(arrayPtr)
	{
		unsigned int i;
		for(i=0;i<rows;i++)
		{
			arrayPtr[i]=(double*)malloc(sizeof *arrayPtr[i]*cols);
			if(arrayPtr[i])
			{
				unsigned int j;
				for(j=0;j<cols;j++)
				{
					arrayPtr[i][j]=0.0;
				}
			}
		}
	}
	return arrayPtr;
}

void free2d(double **arrayPtr,unsigned rows)
{
	//puts("free2d\n");
	if(arrayPtr)
	{
		//printf("root pointer address = %d, number of 1-D pointers = %d\n",long(arrayPtr),rows);
		unsigned i;
		for(i=0;i<rows;i++)
		{
			//printf("node pointer # %d address = %d\n",i,long(arrayPtr[i])); 
			free(arrayPtr[i]);
		}
	}
	free(arrayPtr);
	arrayPtr = NULL;
}
void print_matrix_i(char *name,int **matr,unsigned hights,unsigned length)
{
	unsigned i,j;
	printf("%s: %d x %d\n",name,hights,length);
	if(logfile)
			fprintf(logfile,"%s: %d x %d\n",name,hights,length);
	if(matr && hights>0 && length>0)
	{
		for(i=0;i<hights;i++)
		{
			for(j=0;j<length-1;j++)
			{
				printf("%+d ",matr[i][j]);
				if(logfile)
					fprintf(logfile,"%+d ",matr[i][j]);
			}
			printf("%+d\n",matr[i][length-1]);
			if(logfile)
				fprintf(logfile,"%+d\n",matr[i][length-1]);
		}
		puts("\n");
		if(logfile)
			fprintf(logfile,"\n");
	}
}

void print_matrix_d(char *name,double **matr,unsigned hights,unsigned length)
{
	unsigned i,j;
	printf("%s: %d x %d\n",name,hights,length);
	if(logfile)
			fprintf(logfile,"%s: %d x %d\n",name,hights,length);
	if(matr && hights>0 && length>0)
	{
		for(i=0;i<hights;i++)
		{
			for(j=0;j<length-1;j++)
			{
				printf("%+f ",matr[i][j]);
				if(logfile)
				fprintf(logfile,"%+f ",matr[i][j]);
			}
			printf("%+f\n",matr[i][length-1]);
			if(logfile)
				fprintf(logfile,"%+f\n",matr[i][length-1]);
		}
		puts("\n");
		if(logfile)
			fprintf(logfile,"\n");
	}
}

void print_matrix_d_for_matlab(char *name,double **matr,unsigned hights,unsigned length)
{
	unsigned i,j;
	printf("%s=[...\n",name);
	if(logfile)
			fprintf(logfile,"%s=[...\n",name);
	if(matr)
	{
		for(i=0;i<hights;i++)
		{
			for(j=0;j<length-1;j++)
			{
				printf("%+f ",matr[i][j]);
				if(logfile)
				fprintf(logfile,"%+f ",matr[i][j]);
			}
			if(i<hights-1)
			{
				printf("%+f;...\n",matr[i][length-1]);
				if(logfile)
					fprintf(logfile,"%+f;...\n",matr[i][length-1]);
			}
			else
			{
				printf("%+f]\n",matr[i][length-1]);
				if(logfile)
					fprintf(logfile,"%+f]\n",matr[i][length-1]);			
			}
		}
		puts("\n");
		if(logfile)
			fprintf(logfile,"\n");
	}
}



void print_vector_d_offset(char *name,double *vect,unsigned length,unsigned offset)
{
        unsigned i;
        printf("%s: %d, {%d}\n",name,length,offset);
        if(logfile)
                        fprintf(logfile,"%s: %d, {%d}\n",name,length,offset);
        if(vect && length>0)
        {
                for(i=0;i<length-1;i++)
                {
                        printf("%+f ",vect[(i+offset)%length]);
                        if(logfile)
                                fprintf(logfile,"%+f ",vect[(i+offset)%length]);
                }
                printf("%+f\n\n",vect[(length-1+offset)%length]);
                if(logfile)
                        fprintf(logfile,"%+f\n\n",vect[(length-1+offset)%length]);
        }
}



void print_vector_d(char *name,double *vect,unsigned length)
{
	unsigned i;
	printf("%s: %d\n",name,length);
	if(logfile)
			fprintf(logfile,"%s: %d\n",name,length);
	if(vect && length>0)
	{
		for(i=0;i<length-1;i++)
		{
			printf("%+f ",vect[i]);
			if(logfile)
				fprintf(logfile,"%+f ",vect[i]);
		}
		printf("%+f\n\n",vect[length-1]);
		if(logfile)
			fprintf(logfile,"%+f\n\n",vect[length-1]);
	}
}

void print_vector_i(char *name,int *vect,unsigned length)
{
	unsigned i;
	printf("%s: %d\n",name,length);
	if(logfile)
			fprintf(logfile,"%s: %d\n",name,length);
	if(vect && length>0)
	{
		for(i=0;i<length-1;i++)
		{
			printf("%+d ",vect[i]);
			if(logfile)
				fprintf(logfile,"%+d ",vect[i]);
		}
		printf("%+d\n\n",vect[length-1]);
		if(logfile)
			fprintf(logfile,"%+d\n\n",vect[length-1]);
	}
}

int powi(int arg,int pow)
{
	int i,r=1;
	if(arg==0)
		if(pow==0)
			return 1;
		else
			return 0;
	for(i=0;i<pow;i++)
		r*=arg;
	return r;
}

unsigned powu(unsigned arg,unsigned pow)
{
	unsigned i,r=1;
	if(arg==0)
		if(pow==0)
			return 1;
		else
			return 0;
	for(i=0;i<pow;i++)
		r*=arg;
	return r;
}



unsigned factorizeC(double **matr,double *kernel,unsigned kernel_size_max,int **commutator_image,unsigned m_width,unsigned m_length)
{
	unsigned kernel_size,v,h,v1,h1,vh;
	double **matrix=NULL;
	kernel_size=0;
	matrix=malloc2d(m_width,m_length);
	if(matrix)
	{
		for(v=0;v<m_width;v++)
			for(h=0;h<m_length;h++)
				matrix[v][h]=matr[v][h];
		for(v=0;(v<m_width)&&(kernel_size_max>kernel_size) ;v++)
		{
			for(h=0;h<m_length;h++)
			{	
				if(matrix[v][h])
				{
					if(kernel_size_max>kernel_size)
					{
						kernel[kernel_size]=matrix[v][h];
						for(vh=v*m_length+h;vh<m_width*m_length;vh++)
						{
							v1=vh/m_length;
							h1=vh-v1*m_length;
							if(matrix[v1][h1]==kernel[kernel_size])
							{
								matrix[v1][h1]=0.0;
								commutator_image[v1][h1]=kernel_size+1;
							}
						}
					}
					else
					{
						//printf("Kernel too long, exceeding %d elements\n",kernel_size_max);
						//return 0;
					}
					kernel_size++;
				}
			}
		}
		free2d(matrix,m_width);
	}
	return kernel_size;
}


unsigned countpatternoccurence(int pattern1,int pattern2,int distancebetweenpatternelements,int **matrix,unsigned numberofvectors,unsigned vectorlength)
{
	int tmp[64];
	unsigned vectornumber,position,numberofpatternoccurences;
	numberofpatternoccurences=0;
	for(vectornumber=0;vectornumber<numberofvectors;vectornumber++)
	{	
		for(position=0;position<vectorlength-distancebetweenpatternelements;position++)
		{
			tmp[position]=matrix[vectornumber][position];
		}
		for(position=0;position<vectorlength-distancebetweenpatternelements;position++)
			if(tmp[position]==pattern1)
				if(matrix[vectornumber][position+distancebetweenpatternelements]==pattern2)
				{
					tmp[position]=0;
					numberofpatternoccurences++;
				}
	}
	return numberofpatternoccurences;
}

unsigned counttodesiredpatternoccurence(int pattern1,int pattern2,int distancebetweenpatternelements,int **matrix,unsigned numberofvectors,unsigned vectorlength,unsigned desired_count,int *tmp)
{
	unsigned numberofpatternoccurences = 0;
	if(tmp)
	{
#ifdef __USE_OMP__
#pragma omp parallel for
		for (int vectornumber = 0;vectornumber<numberofvectors;vectornumber++)
		{
			tmp[vectornumber*(vectorlength + 1) + vectorlength] = 0;
			memcpy(&(tmp[vectornumber*(vectorlength+1)]), matrix[vectornumber], sizeof(int)*vectorlength);
			for (unsigned position = vectornumber*(vectorlength + 1);position < vectornumber*(vectorlength + 1) + vectorlength - distancebetweenpatternelements;position++)
				if (tmp[position] == pattern1)
					if (tmp[position + distancebetweenpatternelements] == pattern2)
					{
						tmp[position] = 0;
						tmp[position + distancebetweenpatternelements] = 0;
						tmp[vectornumber*(vectorlength + 1) + vectorlength]++;
					}
		}
		for (int vectornumber = 0;vectornumber<numberofvectors;vectornumber++)
			if((numberofpatternoccurences += tmp[vectornumber*(vectorlength + 1) + vectorlength])>=desired_count)
				return numberofpatternoccurences;
#else
		for(unsigned vectornumber=0;vectornumber<numberofvectors;vectornumber++)
		{
			memcpy(tmp,matrix[vectornumber],sizeof(int)*vectorlength);
			for(unsigned position=0;position<vectorlength-distancebetweenpatternelements;position++)
				if(tmp[position]==pattern1)
					if(tmp[position+distancebetweenpatternelements]==pattern2)
					{
						tmp[position]=0;
						tmp[position+distancebetweenpatternelements]=0;
						if(desired_count<=++numberofpatternoccurences)
							return numberofpatternoccurences;
					}
		}
#endif
	}
	return numberofpatternoccurences;
}


unsigned getmostcommonpatternC(int **matrix,unsigned numberofvectors,unsigned vectorlength,int *mcp)
{
	unsigned numberofterms,vectornumber,position,patternelement1,patternelement2,distancebetweenpatternelements,max_repeat_num,repeat_num;
	numberofterms=0;
	for(vectornumber=0;vectornumber<numberofvectors;vectornumber++)
	{
		for(position=0;position<vectorlength;position++)
		{
			if (numberofterms < unsigned(matrix[vectornumber][position])) numberofterms = unsigned(matrix[vectornumber][position]);
		}
	}

	max_repeat_num=0;
	for(vectornumber=0;vectornumber<numberofvectors;vectornumber++)
	{
		//find the first pattern term
		for(position=0;position<vectorlength;position++)
		{
			if(matrix[vectornumber][position])
			{
				patternelement1=matrix[vectornumber][position];
				//find the second term
				for(distancebetweenpatternelements=1;distancebetweenpatternelements<vectorlength-position;distancebetweenpatternelements++)
				{
					if(matrix[vectornumber][position+distancebetweenpatternelements])
					{
						patternelement2=matrix[vectornumber][position+distancebetweenpatternelements];
						repeat_num=countpatternoccurence(patternelement1,patternelement2,distancebetweenpatternelements,matrix,numberofvectors,vectorlength);
						if(repeat_num>max_repeat_num)
						{
							max_repeat_num=repeat_num;
							mcp[0]=numberofterms+1;
							mcp[1]=patternelement1;
							mcp[2]=patternelement2;
							mcp[3]=distancebetweenpatternelements;
						}
					}
				}
				patternelement1=0;
			}
		}
	}
	//no new patterns.
	return max_repeat_num;
}


unsigned getnextpatternC(int **matrix,unsigned numberofvectors,unsigned vectorlength,int *mcp,unsigned desired_repeat_num,int *tmp)
{
	unsigned numberofterms,vectornumber,position,patternelement1,patternelement2,distancebetweenpatternelements,max_repeat_num,repeat_num;
	numberofterms=0;
	for(vectornumber=0;vectornumber<numberofvectors;vectornumber++)
	{
		for(position=0;position<vectorlength;position++)
		{
			if (numberofterms < unsigned(matrix[vectornumber][position])) numberofterms = unsigned(matrix[vectornumber][position]);
		}
	}

	max_repeat_num=0;
	for(vectornumber=0;vectornumber<numberofvectors;vectornumber++)
	{
		//find the first pattern term
		for(position=0;position<vectorlength;position++)
		{
			if(matrix[vectornumber][position])
			{
				printf("loc = {%d of %d, %d of %d}\n",vectornumber,numberofvectors,position,vectorlength);
				patternelement1=matrix[vectornumber][position];
				//find the second term
				for(distancebetweenpatternelements=1;distancebetweenpatternelements<vectorlength-position;distancebetweenpatternelements++)
				{
					if(matrix[vectornumber][position+distancebetweenpatternelements])
					{
						patternelement2=matrix[vectornumber][position+distancebetweenpatternelements];
						repeat_num=1;
						if(desired_repeat_num>1)
							repeat_num=counttodesiredpatternoccurence(patternelement1,patternelement2,distancebetweenpatternelements,matrix,numberofvectors,vectorlength,desired_repeat_num,tmp);
						if(repeat_num>max_repeat_num)
						{
							max_repeat_num=repeat_num;
							mcp[0]=numberofterms+1;
							mcp[1]=patternelement1;
							mcp[2]=patternelement2;
							mcp[3]=distancebetweenpatternelements;
							if(repeat_num>=desired_repeat_num)
								return max_repeat_num;
						}
					}
				}
				patternelement1=0;
			}
		}
	}
	//no new patterns.
	return max_repeat_num;
}

unsigned replacepatternC(int *pattern,int **matrix,unsigned numberofvectors,unsigned vectorlength)
{
	//unsigned
	int vectornumber,position,numberofpatternoccurences;
	int distancebetweenpatternelements;
	numberofpatternoccurences=0;
	distancebetweenpatternelements=pattern[3];
#ifdef __USE_OMP__
#pragma omp parallel for
#endif
	for(vectornumber=0;vectornumber<numberofvectors;vectornumber++)
		for(position=0;position<vectorlength-distancebetweenpatternelements;position++)
			if(matrix[vectornumber][position]==pattern[1])
				if(matrix[vectornumber][position+distancebetweenpatternelements]==pattern[2])
				{
					matrix[vectornumber][position]=0;
					matrix[vectornumber][position+distancebetweenpatternelements]=pattern[0];
					numberofpatternoccurences++;
				}
	return numberofpatternoccurences;
}

unsigned reducetermsC(int **matrix,unsigned numberofvectors,unsigned vectorlength,int **patterns,unsigned pattern_size_max,int *occurences,unsigned dropthreshold)
{
	// pattern format: razpatnumber patternelementlate patternelementearly distancebetweenpatternelements
	unsigned pattern_size, numberofpatternoccurences, desired_occurence;//, dropthreshold;
	int *tmp=NULL;
	desired_occurence=100000;
	//dropthreshold=8;
	pattern_size=0;
#ifdef __USE_OMP__
	tmp=(int *)malloc(sizeof(int)*(vectorlength+1)*numberofvectors);
#else
	tmp = (int *)malloc(sizeof(int)*vectorlength);
#endif
	if(tmp)
	{
		numberofpatternoccurences=1;
		while(numberofpatternoccurences>0)
		{
			numberofpatternoccurences=getnextpatternC(matrix,numberofvectors,vectorlength,patterns[pattern_size],desired_occurence,tmp);
			if(numberofpatternoccurences>0)
			{
				if(desired_occurence)
				{
					desired_occurence=numberofpatternoccurences;
					if(numberofpatternoccurences<=dropthreshold)
						desired_occurence=0;
				}
#if (1)
				print_vector_i("Pattern to replace",patterns[pattern_size],5);
				//print_matrix_i("Commutator image befor the pattern is replaced",matrix,numberofvectors,vectorlength);
#endif
				occurences[pattern_size]=replacepatternC(patterns[pattern_size],matrix,numberofvectors,vectorlength);
#if (1)
				printf("Pattern # %d encountered %d times\n",1+pattern_size,occurences[pattern_size]);
				//print_matrix_i("Commutator image with replaced pattern",matrix,numberofvectors,vectorlength);
#endif
				pattern_size++;
				if(pattern_size_max<=pattern_size)
				{
					printf("too many patterns %d\n",pattern_size);
					return 0;
				}
			}
		}
		free(tmp);
	}
	return pattern_size;
}

int makemultichannelpatterns(int **patterns,unsigned pattern_size,unsigned number_of_channels)
{
	unsigned i,j;
	if(!patterns || pattern_size==0)
		return -2;
	if(number_of_channels>1)
	{
		for(i=0;i<pattern_size;i++)
			for(j=3;j<=4;j++)
				patterns[i][j]*=number_of_channels;
		return 0;
	}
	return -1;
}

//.//.//.//.//.//.//.//.//.//.//.//.//.//.//
int retardarguments(int **patterns,unsigned pattern_size,unsigned pipelinedelay,int *time_aligning_delays)
{
	unsigned i,j;
	int k,l;
	l=0;
#if (0)
	print_matrix_i("initial patterns",patterns,pattern_size,5);
#endif
	for(i=1;i<=pattern_size;i++)
	{
		k=0;
		for(j=2;j<=3;j++)
			if(time_aligning_delays[patterns[i-1][j-1]-1]>k)
			   k=time_aligning_delays[patterns[i-1][j-1]-1];
		time_aligning_delays[i-1+patterns[0][0]-1]=k+pipelinedelay;
		for(j=1;j<=2;j++)
			patterns[i-1][j-1+3]+=k-time_aligning_delays[patterns[i-1][1+j-1]-1];
		j= (patterns[i-1][3] < patterns[i-1][4]) ? patterns[i-1][3] : patterns[i-1][4];
		time_aligning_delays[i-1+patterns[0][0]-1]-=j;
		if (l < time_aligning_delays[i-1+patterns[0][0]-1]) l = time_aligning_delays[i-1+patterns[0][0]-1];
		for(k=3;k<=4;k++)
			patterns[i-1][k]-=j;
	}
#if (0)
	print_matrix_i("patterns for components with delay",patterns,pattern_size,5);
#endif
	for(i=0;i<unsigned(patterns[pattern_size-1][0]);i++)
		time_aligning_delays[i]=l-time_aligning_delays[i];
#if (0)
	print_vector_i("time_aligning_delays",time_aligning_delays,patterns[pattern_size-1][0]);
#endif
	return l;
}
//.//.//.//.//.//.//.//.//.//.//.//.//.//.//
int calculatetotaldelay(int pattern,int **patterns,unsigned pattern_size,unsigned pipelinedelay,int *time_aligning_delays)
{
	int i,j,b,d,p,delay[2]={0,0};
	b=patterns[0][0];
	j=pattern-b;
	if(j>=0)
	{
		for(i=0;i<2;i++)
		{
			p=patterns[j][i+1];
			d=patterns[j][i+3];
			if(p>=b)
				delay[i]=d+calculatetotaldelay(p,patterns,pattern_size,pipelinedelay,time_aligning_delays);
		}
		if(delay[0]<delay[1])
			delay[0]=delay[1];
		return pipelinedelay+delay[0];
	}
	return 0;
}
///..///..///..///..///..///..///..///..////

void update_pattern_timescale(int timescale, int ** patterns, unsigned pattern_size, int *spectrumtermtaps, unsigned width)
{
	if (timescale <= 0)
	{
		if (patterns && pattern_size)
		{
			for (unsigned i = 0;i < pattern_size; i++)
				for (unsigned j = 3;j < 5;j++)
				{
					patterns[i][j] = 0;
				}
		}
		if (width && spectrumtermtaps)
		{
			for (unsigned i = 0;i < width;i++)
				spectrumtermtaps[i] = 0;
		}
	}
}

unsigned serialize_commutatorC(int **commutator_image,unsigned numberofvectors,unsigned vectorlength,unsigned pipelinedelay,unsigned number_of_channels,int **patterns,unsigned pattern_size_max,int *patternoccurences,int *spectrumterms,int *spectrumtermtaps,int *extention, unsigned dropthreshold)
{
	unsigned pattern_size,m,n,conveyor_size,max_conveyor_size;
	int min;
	int *time_aligning_delays=NULL;
	pattern_size=reducetermsC(commutator_image,numberofvectors,vectorlength,patterns,pattern_size_max,patternoccurences,dropthreshold);
#if (0)
    print_matrix_i("Preliminary pattern set",patterns,pattern_size,5);
#endif
	makemultichannelpatterns(patterns,pattern_size,number_of_channels);
#if (0)
	print_matrix_i("Preliminary pattern set with correct number of channels",patterns,pattern_size,5);
#endif
	time_aligning_delays=(pattern_size>0)?malloc1i(patterns[pattern_size-1][0]):NULL;
	//if(time_aligning_delays)
	{
		max_conveyor_size= time_aligning_delays==NULL?0:retardarguments(patterns,pattern_size,pipelinedelay,time_aligning_delays);
	    print_matrix_i("Updated pattern set with retarded delays",patterns,pattern_size,5);
		
		min=vectorlength+max_conveyor_size;
		for(m=0;m<numberofvectors;m++)
		{
			spectrumterms[m]=0;
			for(n=0;n<vectorlength;n++)
				if(commutator_image[m][n])
				{
					spectrumterms[m]=commutator_image[m][n];
					spectrumtermtaps[m] = (time_aligning_delays==NULL)?n:(n-time_aligning_delays[spectrumterms[m]-1]);
					if (min > spectrumtermtaps[m]) min = spectrumtermtaps[m];
					conveyor_size= (time_aligning_delays==NULL)?0:time_aligning_delays[spectrumterms[m]-1]+calculatetotaldelay(spectrumterms[m],patterns,pattern_size,pipelinedelay,time_aligning_delays);
#if (0)
					printf("conveyor_size = %d\n",conveyor_size);
#endif
					if (max_conveyor_size < conveyor_size) max_conveyor_size = conveyor_size;
				}
		}
		if (time_aligning_delays!=(int*)NULL)
		{
			free(time_aligning_delays);
			time_aligning_delays = (int*)NULL;
		}
	}
#if (0)
	printf("max conveyor size = %d\n",max_conveyor_size);
#endif
	*extention=max_conveyor_size;
	for(m=0;m<numberofvectors;m++)
		spectrumtermtaps[m]-=min;
#if (0)
	print_vector_i("Spectrum terms",spectrumterms,numberofvectors);
	print_vector_i("Spectrum term taps",spectrumtermtaps,numberofvectors);
#endif
    return pattern_size;
}
