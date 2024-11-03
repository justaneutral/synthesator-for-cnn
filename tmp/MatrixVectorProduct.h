#include "factorizator.h"

//#define __AUTONOMOUS_INDEXING__

class MatrixVectorProduct
{
public:
int initialized;
double **matrix;
unsigned matrix_width;
unsigned matrix_length;
double **terms;
unsigned conveyerlength;
unsigned conveyerwidth;
int conveyerindex;
double *kernel;
unsigned kernel_size;
int *spectrumterms;
int *spectrumtermtaps;
int **patterns;
unsigned pattern_size;
double *spectrum;
double *testspectrum;
double *vector;
unsigned vectorlength;
double *testvector;
unsigned testvectorlength;
int offset;
double currentelement;
unsigned extention;
int numberofiterations;


MatrixVectorProduct(void)
{
	initialized = (int)0;
	matrix = (double **)NULL;
	matrix_width = (unsigned)0;
	matrix_length = (unsigned)0;
	terms = (double **)NULL;
	conveyerlength = (unsigned)0;
	conveyerwidth = (unsigned)0;
		conveyerindex = (int)0;
	kernel = (double *)NULL;
	kernel_size = (unsigned)0;
	spectrumterms = (int *)NULL;
	spectrumtermtaps = (int *)NULL;
	patterns = (int**)NULL;
	pattern_size = (unsigned)0;
	spectrum = (double *)NULL;
	testspectrum = (double *)NULL;
	vector = (double *)NULL;
	vectorlength = (unsigned)0;
	testvector = (double *)NULL;
	testvectorlength = (unsigned)0;
		offset = (int)0;
	currentelement = (double)0.0;
	extention = (unsigned)0;
	numberofiterations = (int)0;
}




int initialize(double **matrix, unsigned matrix_width, unsigned matrix_length, double *kernel, unsigned kernel_size, unsigned **patterns, unsigned pattern_size, double **terms, int *spectrumtermtaps, unsigned conveyerlength, unsigned conveyerwidth, unsigned extention) 
{
	if(!initialized)
	{
		this->matrix = matrix;
		this->matrix_width = matrix_width;
		this->matrix_length = matrix_length;

		this->kernel = kernel;
		this->kernel_size = kernel_size;

		this->patterns = patterns;
		this->pattern_size = pattern_size;

		this->terms = terms;
		

		this->conveyerlength = conveyerlength;
		this->spectrum = spectrum;
		this->spectrumterms = spectrumterms;
		this->spectrumtermtaps = spectrumtermtaps;		
		this->matrix_width = matrix_width;
		this->conveyerwidth = conveyerwidth;
		this->conveyerlength = conveyerlength;
		this->extention = extention;

	this->extention = extention;
		this->spectrumtermtaps = spectrumtermtaps;
		this->kernel_size = kernel_size;
		this->pattern_size = pattern_size;
		this->conveyerlength=this->matrix_length+this->extention;
#if (1)
        	print_vector_i("taps",this->spectrumtermtaps,this->matrix_width);
#endif
        	this->conveyerwidth=this->kernel_size+this->pattern_size;
        	terms=malloc2d(conveyerwidth,conveyerlength);
        	if(terms)
        	{
			initialized = 1;
		}
	}
	return initialized;
}


int uninitialize(void)
{
	int rc = initialized;
	if(initialized)
	{
	//print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
	}
	//free buffers that are not used
	if(terms && conveyerwidth > 0)
	{
		free2d(terms,conveyerwidth);
        }
	initialized = 0;
	return rc;
}



int matrix_multiplication_conventional(double **matrix,unsigned matrix_width,unsigned matrix_length,double *spectrum,double *testvector,unsigned testvectorlength,int offset)
{
	unsigned ind1,ind2;
	int ind3;
	for(ind1=0;ind1<matrix_width;ind1++)
	{
		spectrum[ind1]= 0.0;
	}
	for (ind2 = 0;ind2 < matrix_length;ind2++)
	{
		for(ind1=0;ind1<matrix_width;ind1++)
		{
			ind3 = offset + ind2;
			if (ind3 >= 0 && ind3<int(testvectorlength))
			{
				spectrum[ind1] += matrix[ind1][ind2] * testvector[ind3];
			}
			else
			{
				print_vector_i("Error - test vector preambula is too short",NULL,ind3);
				return -1;
			}
		}
	}
	return 0;
}

int matrix_multiplication_conventional_for_offset_search(double **matrix,unsigned matrix_width,unsigned matrix_length,double *spectrum)
{
	unsigned ind1,ind2;
	ind2=0;
	for(ind1=0;ind1<matrix_width;ind1++)
	{
		spectrum[ind1]=matrix[ind1][matrix_length-1];
		if(spectrum[ind1]!=0.0)
			ind2++;
	}
	return ind2;
}

void matrix_multiplication_iterationC(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement,int conveyerindex)
{
#ifdef __AUTONOMOUS_INDEXING__
        static int x = conveyerlength-1; //stronger test
        x = x >= int(conveyerlength)-1 ? 0 : x+1;
#else //__AUTONOMOUS_INDEXING__
        int x = (conveyerindex+1)%conveyerlength;
#endif //__AUTONOMOUS_INDEXING__
	int insert_position;
	int i,j,termsrow,patternsrow,p1,p1delay,p1pos,p2,p2delay,p2pos;
	double arg1,arg2;
        insert_position = conveyerlength-1 + x;
        if (insert_position < 0) insert_position += conveyerlength;
        else if (insert_position >= conveyerlength) insert_position -= conveyerlength;
	//populate products
    for(i=0;i<kernel_size;i++)
		terms[i][insert_position]=currentelement*kernel[i];
    for(i=kernel_size;i<int(conveyerwidth);i++)
		terms[i][insert_position]=0;
    //here we update terms one by one if they are.
	if (pattern_size > 0)
	{
		for (termsrow = patterns[0][0] - 1;termsrow < patterns[pattern_size - 1][0];termsrow++)
		{
			patternsrow = termsrow - patterns[0][0] + 1;
			p1 = patterns[patternsrow][1] - 1;
			p1delay = patterns[patternsrow][3];
			p1pos = conveyerlength - p1delay - 1;
			p1pos = p1pos + x;
			if (p1pos < 0) p1pos += conveyerlength;
			else if (p1pos >= conveyerlength) p1pos -= conveyerlength;
			arg1 = terms[p1][p1pos];
			p2 = patterns[patternsrow][2] - 1;
			p2delay = patterns[patternsrow][4];
			p2pos = conveyerlength - p2delay - 1;
			p2pos = p2pos + x;
			if (p2pos < 0) p2pos += conveyerlength;
			else if (p2pos >= conveyerlength) p2pos -= conveyerlength;
			arg2 = terms[p2][p2pos];
			terms[termsrow][insert_position] = arg1 + arg2;
		}
	}
}

void extract_spectrumC(double *spectrum,double **terms,int *spectrumterms,int *spectrumtermtaps,unsigned matrix_width,unsigned conveyerlength, int conveyerindex)
{
#ifdef __AUTONOMOUS_INDEXING__
        static int x = conveyerlength-1; //stronger test
        x = x >= int(conveyerlength)-1 ? 0 : x+1;
	printf("extract spectrum: x=%d, conveyerindex=%d\n", x, (conveyerindex+1)%conveyerlength);
#else //__AUTONOMOUS_INDEXING__
	int x = (conveyerindex+1)%conveyerlength;
#endif //__AUTONOMOUS_INDEXING__
	unsigned spectrumindex,spectrumtermindex;
	for(spectrumindex=0;spectrumindex<matrix_width;spectrumindex++)
	{
		if (spectrumterms[spectrumindex] > 0)
		{
			spectrumtermindex = (spectrumterms[spectrumindex] - 1);
			int spectrumtermtap = spectrumtermtaps[spectrumindex] + x;
			if (spectrumtermtap < 0) spectrumtermtap += conveyerlength;
                        else if (spectrumtermtap >= conveyerlength) spectrumtermtap -= conveyerlength;
			spectrum[spectrumindex] = terms[spectrumtermindex][spectrumtermtap];
		}
		else
		{
			spectrum[spectrumindex] = 0;
			//printf(".");
		}
	}
}



#define __SHOW_RESULTS_FOR_EACH_SHIFT__
void matrix_multiplicationConveyor(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned matrix_width,unsigned matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,int numberofiterations)
{
	unsigned conveyerlength,conveyerwidth,vectorindex,maxindex;
	double currentelement;
	double **terms=NULL;
	conveyerlength=matrix_length+extention;
	maxindex = numberofiterations + extention >= vectorlength ? vectorlength : numberofiterations + extention;
	//maxindex = numberofiterations + matrix_length - 1 >= vectorlength ? vectorlength : numberofiterations + matrix_length - 1;
    	conveyerwidth=kernel_size+pattern_size;
    	terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
		for(vectorindex=0;vectorindex<maxindex;vectorindex++)
		{
			currentelement=vector[vectorindex];
			matrix_multiplication_iterationC(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,vectorindex);
			//print_matrix_d("Terms",terms,conveyerwidth,conveyerlength);
#ifdef __SHOW_RESULTS_FOR_EACH_SHIFT__
			printf("%d:",vectorindex);
			extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,matrix_width,conveyerlength,vectorindex);
			print_vector_d("cur spec conveyor:",spectrum,matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
		}
		//print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
#ifndef __SHOW_RESULTS_FOR_EACH_SHIFT__
		extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,matrix_width,conveyerlength,maxindex-1);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
		print_vector_d("cur spec conveyor:",spectrum,matrix_width);
		free2d(terms,conveyerwidth);
	}
}

void matrix_multiplicationOneVector(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned matrix_width,unsigned matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,int numberofiterations)
{
/*
	unsigned conveyerlength,conveyerwidth,i,j,maxindex;
	double currentelement;
	double **terms=NULL;
	conveyerlength=matrix_length+extention;
*/
	unsigned maxindex = numberofiterations + extention >= vectorlength ? vectorlength : numberofiterations + extention;
	//maxindex = (numberofiterations + matrix_length - 1) >= vectorlength ? vectorlength : (numberofiterations + matrix_length - 1);

#if (1)
	print_vector_i("taps",spectrumtermtaps,matrix_width);
#endif
/*
	conveyerwidth=kernel_size+pattern_size;
        terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
*/
		//feed all vector elements to the conveyor and calculate intermittant terms
		for(int i=0;i<maxindex;i++)
		{
			currentelement=vector[i]; //can be any or from different channels mixed in the input stream
			matrix_multiplication_iterationC(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,i);
#ifdef __SHOW_RESULTS_FOR_EACH_SHIFT__
			printf("%d:",i);
			extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,matrix_width,conveyerlength,i);
			print_vector_d("cur spec:",spectrum,matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
			//print_matrix_d("Terms",terms,conveyerwidth,conveyerlength);
		}
/*
		//print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
#ifndef __SHOW_RESULTS_FOR_EACH_SHIFT__
		extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,matrix_width,conveyerlength,maxindex-1);
		print_vector_d("cur spec:",spectrum,matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
		//free buffers that are not used
		free2d(terms,conveyerwidth);
	}
*/
}

double calculate_energy(double *vect0,double *vect1, unsigned len)
{
	double energy;
	unsigned i;
	energy=0.0;
	for (i = 0;i < len;i++)
	{
		energy += vect0[i] * vect1[i];
		if (energy < 0.0)
		{
			printf("energy < 0\n");
			print_vector_d("vect0", vect0, len);
			print_vector_d("vect1", vect1, len);
		}
	}
	return energy;
}

unsigned matrix_multiplication_find_offset(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned matrix_width,unsigned matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,double *testspectrum)
{
	unsigned conveyerlength,conveyerwidth,i,taillength,totaloffset;
	double currentelement,energy0,energy1;
	double **terms=NULL;
	totaloffset=0;
	conveyerlength=matrix_length+extention;
	taillength=matrix_length*matrix_width*conveyerlength;//6 for 2 ch 0pd.;//3 for 1 ch. 0pd.
    
	conveyerwidth=kernel_size+pattern_size;
        terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
		//calculate criteria for stopping
		energy0=calculate_energy(testspectrum,testspectrum,matrix_width);
		energy1=0.0;
		//feed all vector elements to the conveyor and calculate intermittent terms
		//continue cycling until all terms travel to the terminals
		for(i=0;energy0>=0 && i<vectorlength+taillength;i++)
		{
			totaloffset++;
			currentelement = (i == 0) ? 1.0 : 0.0;
			matrix_multiplication_iterationC(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,i);
			//if(i>=conveyerlength)
			{
				extract_spectrumC(spectrum, terms, spectrumterms, spectrumtermtaps, matrix_width, conveyerlength, i);
				energy1 = calculate_energy(testspectrum, spectrum, matrix_width);
				if (energy1 || energy1*1.01 > energy0)
				{
					energy0 = -1.0;
				}
			}
		}
		//free buffers that are not used
		free2d(terms,conveyerwidth);
	}
	return totaloffset-1;
}


}
