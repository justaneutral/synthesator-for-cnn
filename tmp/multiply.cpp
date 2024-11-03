#include "factorizator.h"




int matrix_multiplication_conventional_step(double **matrix,unsigned matrix_width,unsigned matrix_length,double *originalspectrum,double *testvector,int *offset, int GF2, int currentvalue)
{
	unsigned ind0,ind1,ind2,ind3;
	for(ind1=0;ind1<matrix_width;ind1++)
	{
		originalspectrum[ind1]= (GF2 < 0)? 1.0 : 0.0;
	}

	ind0 = ((*offset)+matrix_length-1)%matrix_length;
	*offset=ind0;
	testvector[ind0] = currentvalue;

	for (ind2 = 0;ind2 < matrix_length;ind2++)
	{
		for(ind1=0;ind1<matrix_width;ind1++)
		{
			ind3 = (ind0 + ind2) % matrix_length;
			if (GF2 == 0)
			{
				originalspectrum[ind1] += matrix[ind1][ind2] * testvector[ind3];
			}
			else
			{
				if (GF2 < 0)
				{
					if (matrix[ind1][ind2])
					{
						originalspectrum[ind1] *= testvector[ind3];
					}
				}
				else
				{
					unsigned long a = originalspectrum[ind1], b = matrix[ind1][ind2], c = testvector[ind3];
					a ^= (b&c);
					originalspectrum[ind1] = a;
				}
			}
		}
	}
	return 0;
}



int matrix_multiplication_conventional(double **matrix,unsigned matrix_width,unsigned matrix_length,double *originalspectrum,double *testvector,unsigned testvectorlength,int offset, int GF2)
{
	unsigned ind1,ind2;
	int ind3;
	for(ind1=0;ind1<matrix_width;ind1++)
	{
		originalspectrum[ind1]= (GF2 < 0)? 1.0 : 0.0;
	}

	for (ind2 = 0;ind2 < matrix_length;ind2++)
	{
		for(ind1=0;ind1<matrix_width;ind1++)
		{
			ind3 = offset + ind2;
			if (ind3 >= 0 && ind3<int(testvectorlength))
			{
				if (GF2 == 0)
				{
					originalspectrum[ind1] += matrix[ind1][ind2] * testvector[ind3];
				}
				else
				{
					if (GF2 < 0)
					{
						if (matrix[ind1][ind2])
						{
							originalspectrum[ind1] *= testvector[ind3];
						}
					}
					else
					{
						unsigned long a = originalspectrum[ind1], b = matrix[ind1][ind2], c = testvector[ind3];
						a ^= (b&c);
						originalspectrum[ind1] = a;
					}
				}
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

int matrix_multiplication_conventional_for_offset_search(double **matrix,unsigned matrix_width,unsigned matrix_length,double *originalspectrum)
{
	unsigned ind1,ind2;
	ind2=0;
	for(ind1=0;ind1<matrix_width;ind1++)
	{
		originalspectrum[ind1]=matrix[ind1][matrix_length-1];
		if(originalspectrum[ind1]!=0.0)
			ind2++;
	}
	return ind2;
}

void matrix_multiplication_iterationC(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement,double **rez1,unsigned pipelinedelay,int GF2)
{
        //static int x = 0;
        static int x = conveyerlength-1; //stronger test
	int insert_position;
	int i,j,termsrow,patternsrow,p1,p1delay,p1pos,p2,p2delay,p2pos;
	double arg1,arg2,rez;
	
	//propagate
        //x = x == 0 ? int(conveyerlength)-1 : x-1;
        x = x >= int(conveyerlength)-1 ? 0 : x+1;
        insert_position = conveyerlength-1 + x;
        if (insert_position < 0) insert_position += conveyerlength;
        else if (insert_position >= conveyerlength) insert_position -= conveyerlength;
	/*for(i=0;i<int(conveyerwidth);i++)
	{
		for(j=0;j<int(conveyerlength)-1;j++)
			terms[i][j]=terms[i][j+1];
		terms[i][conveyerlength-1]=0;
	}*/
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

			///.///good for pipeline delay 
			if (pipelinedelay > 0)
				for (i = pipelinedelay - 1;i >= 0;i--)
					rez1[termsrow - patterns[0][0] + 1][i + 1] = rez1[termsrow - patterns[0][0] + 1][i];
			if (GF2 == 0)
			{
				rez1[termsrow - patterns[0][0] + 1][0] = arg1 + arg2;
			}
			else
			{
				if (GF2 < 0)
				{
					rez1[termsrow - patterns[0][0] + 1][0] = arg1 * arg2;
				}
				else
				{
					unsigned long a, b=arg1, c=arg2;
					a = b^c;
					rez1[termsrow - patterns[0][0] + 1][0] = a;
				}
			}
			///.///rez=arg1+arg2;	
			rez = rez1[termsrow - patterns[0][0] + 1][pipelinedelay];
			terms[termsrow][insert_position] = rez;
		}
	}
}


void matrix_multiplication_iterationCI(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement,double **rez1,unsigned pipelinedelay,int GF2,int conveyerindex)
{
        int i,termsrow,patternsrow,p1,p1delay,p1pos,p2,p2delay,p2pos;
        double arg1,arg2,rez;

#ifdef __AUTONOMOUS_INDEXING__
        static int x = conveyerlength-1; //stronger test
        x = x >= int(conveyerlength)-1 ? 0 : x+1;
#else //__AUTONOMOUS_INDEXING__
        int x = (conveyerindex+1)%conveyerlength;
#endif //__AUTONOMOUS_INDEXING__

        int insert_position = (conveyerlength-1 + x)%conveyerlength;

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

			///.///good for pipeline delay 
			if (pipelinedelay > 0)
				for (i = pipelinedelay - 1;i >= 0;i--)
					rez1[termsrow - patterns[0][0] + 1][i + 1] = rez1[termsrow - patterns[0][0] + 1][i];
			if (GF2 == 0)
			{
				rez1[termsrow - patterns[0][0] + 1][0] = arg1 + arg2;
			}
			else
			{
				if (GF2 < 0)
				{
					rez1[termsrow - patterns[0][0] + 1][0] = arg1 * arg2;
				}
				else
				{
					unsigned long a, b=arg1, c=arg2;
					a = b^c;
					rez1[termsrow - patterns[0][0] + 1][0] = a;
				}
			}
			///.///rez=arg1+arg2;	
			rez = rez1[termsrow - patterns[0][0] + 1][pipelinedelay];
			terms[termsrow][insert_position] = rez;
		}
	}
}


void matrix_multiplication_iteration(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement, int conveyerindex)
{
        int i,termsrow,patternsrow,p1,p1delay,p1pos,p2,p2delay,p2pos;
        double arg1,arg2,rez;

#ifdef __AUTONOMOUS_INDEXING__
        static int x = conveyerlength-1; //stronger test
        x = x >= int(conveyerlength)-1 ? 0 : x+1;
#else //__AUTONOMOUS_INDEXING__
        int x = (conveyerindex+1)%conveyerlength;
#endif //__AUTONOMOUS_INDEXING__

        int insert_position = (conveyerlength-1 + x)%conveyerlength;
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





void extract_spectrumC(double *spectrum,double **terms,int *spectrumterms,int *spectrumtermtaps,unsigned numberofvectors,unsigned conveyerlength, int conveyerindex)
{
#ifdef __AUTONOMOUS_INDEXING__
        static int x = conveyerlength-1; //stronger test
        x = x >= int(conveyerlength)-1 ? 0 : x+1;
	printf("extract spectrum: x=%d, conveyerindex=%d\n", x, (conveyerindex+1)%conveyerlength);
#else //__AUTONOMOUS_INDEXING__
	int x = (conveyerindex+1)%conveyerlength;
#endif //__AUTONOMOUS_INDEXING__
	unsigned spectrumindex,spectrumtermindex;
	for(spectrumindex=0;spectrumindex<numberofvectors;spectrumindex++)
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



void matrix_multiplicationConveyorC(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,int GF2,int numberofiterations)
{
	unsigned conveyerlength,conveyerwidth,vectorindex,maxindex;
	double currentelement;
	double **terms=NULL;
	double **rez1=NULL;
	conveyerlength=original_matrix_length+extention;
	//maxindex = numberofiterations + extention >= vectorlength ? vectorlength : numberofiterations + extention;
	maxindex = numberofiterations + original_matrix_length - 1 >= vectorlength ? vectorlength : numberofiterations + original_matrix_length - 1;
    	conveyerwidth=kernel_size+pattern_size;
    	terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
		if(rez1=malloc2d(pattern_size,pipelinedelay+1))
		{
			for(vectorindex=0;vectorindex<maxindex;vectorindex++)
			{
				currentelement=vector[vectorindex];
				matrix_multiplication_iterationC(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,rez1,pipelinedelay,GF2);
				//print_matrix_d("Terms",terms,conveyerwidth,conveyerlength);
#ifdef __SHOW_RESULTS_FOR_EACH_SHIFT__
				printf("%d:",vectorindex);
				extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,vectorindex);
				print_vector_d("cur spec conveyor:",spectrum,original_matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
			}
			//print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
#ifndef __SHOW_RESULTS_FOR_EACH_SHIFT__
			extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,maxindex-1);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
			print_vector_d("cur spec conveyor:",spectrum,original_matrix_width);
			free2d(rez1,pattern_size);
		}
		free2d(terms,conveyerwidth);
	}
}

void matrix_multiplicationConveyorCI(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,int GF2,int numberofiterations)
{
	unsigned conveyerlength,conveyerwidth,vectorindex,maxindex;
	double currentelement;
	double **terms=NULL;
	double **rez1=NULL;
	conveyerlength=original_matrix_length+extention;
	//maxindex = numberofiterations + extention >= vectorlength ? vectorlength : numberofiterations + extention;
	maxindex = numberofiterations + original_matrix_length - 1 >= vectorlength ? vectorlength : numberofiterations + original_matrix_length - 1;
    	conveyerwidth=kernel_size+pattern_size;
    	terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
		if(rez1=malloc2d(pattern_size,pipelinedelay+1))
		{
			for(vectorindex=0;vectorindex<maxindex;vectorindex++)
			{
				currentelement=vector[vectorindex];
				matrix_multiplication_iterationCI(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,rez1,pipelinedelay,GF2,vectorindex);
				//print_matrix_d("Terms",terms,conveyerwidth,conveyerlength);
#ifdef __SHOW_RESULTS_FOR_EACH_SHIFT__
				printf("%d:",vectorindex);
				extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,vectorindex);
				print_vector_d("cur spec conveyor:",spectrum,original_matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
			}
			//print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
#ifndef __SHOW_RESULTS_FOR_EACH_SHIFT__
			extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,maxindex-1);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
			print_vector_d("cur spec conveyor:",spectrum,original_matrix_width);
			free2d(rez1,pattern_size);
		}
		free2d(terms,conveyerwidth);
	}
}


void matrix_multiplicationConveyor(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,int numberofiterations)
{
        unsigned conveyerlength,conveyerwidth,vectorindex,maxindex;
        double currentelement;
        double **terms=NULL;
        conveyerlength=original_matrix_length+extention;
        maxindex = numberofiterations + original_matrix_length - 1 >= vectorlength ? vectorlength : numberofiterations + original_matrix_length - 1;
        conveyerwidth=kernel_size+pattern_size;
        terms=malloc2d(conveyerwidth,conveyerlength);
        if(terms)
        {
                for(vectorindex=0;vectorindex<maxindex;vectorindex++)
                {
                        currentelement=vector[vectorindex];
                        matrix_multiplication_iteration(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,vectorindex);
                        //print_matrix_d("Terms",terms,conveyerwidth,conveyerlength);
#ifdef __SHOW_RESULTS_FOR_EACH_SHIFT__
                        printf("%d:",vectorindex);
                        extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,vectorindex);
                        print_vector_d("cur spec conveyor:",spectrum,original_matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
                }
                //print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
#ifndef __SHOW_RESULTS_FOR_EACH_SHIFT__
                extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,maxindex-1);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
                print_vector_d("cur spec conveyor:",spectrum,original_matrix_width);
		free2d(terms,conveyerwidth);
        }
}


void matrix_multiplicationOneVector(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,unsigned number_of_channels,int GF2,int numberofiterations)
{
	unsigned conveyerlength,conveyerwidth,vectorindex,maxindex;
	double currentelement;
	double **terms=NULL;
	double **rez1=NULL;
	conveyerlength=(original_matrix_length)*number_of_channels+extention;
	maxindex = number_of_channels*(numberofiterations + original_matrix_length - 1) >= vectorlength ? vectorlength : number_of_channels*(numberofiterations + original_matrix_length - 1);
#if (1)
	print_vector_i("taps",spectrumtermtaps,original_matrix_width);
#endif
	conveyerwidth=kernel_size+pattern_size;
        terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
		if(rez1=malloc2d(pattern_size,pipelinedelay+1))
		{
			//feed all vector elements to the conveyor and calculate intermittant terms
			for(vectorindex = 0; vectorindex < maxindex; vectorindex++)
			{
				currentelement=vector[vectorindex]; //can be any or from different channels mixed in the input stream
				matrix_multiplication_iterationC(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,rez1,pipelinedelay,GF2);
#ifdef __SHOW_RESULTS_FOR_EACH_SHIFT__
				printf("%d:%d:",vectorindex/number_of_channels,vectorindex%number_of_channels);
				extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,vectorindex);
				print_vector_d("cur spec:",spectrum,original_matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
				//print_matrix_d("Terms",terms,conveyerwidth,conveyerlength);
			}
			//free buffers that are not used
			free2d(rez1,pattern_size);
			//print_matrix_d("Terms, final stage",terms,conveyerwidth,conveyerlength);
#ifndef __SHOW_RESULTS_FOR_EACH_SHIFT__
			extract_spectrumC(spectrum,terms,spectrumterms,spectrumtermtaps,original_matrix_width,conveyerlength,maxindex-1);
			print_vector_d("cur spec:",spectrum,original_matrix_width);
#endif // __SHOW_RESULTS_FOR_EACH_SHIFT__
		}
		free2d(terms,conveyerwidth);
	}
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

unsigned matrix_multiplication_find_offset(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,unsigned number_of_channels,double *testspectrum)
{
	unsigned conveyerlength,conveyerwidth,i,taillength,totaloffset;
	double currentelement,energy0,energy1;
	double **terms=NULL;
	double **rez1=NULL;
	totaloffset=0;
	conveyerlength=(original_matrix_length)*number_of_channels+extention;
	taillength=original_matrix_length*original_matrix_width*conveyerlength;//6 for 2 ch 0pd.;//3 for 1 ch. 0pd.
    
	conveyerwidth=kernel_size+pattern_size;
    terms=malloc2d(conveyerwidth,conveyerlength);
	if(terms)
	{
		if(rez1=malloc2d(pattern_size,pipelinedelay+1))
		{
			//calculate criteria for stopping
			energy0=calculate_energy(testspectrum,testspectrum,original_matrix_width);
			energy1=0.0;
			//feed all vector elements to the conveyor and calculate intermittent terms
			//continue cycling until all terms travel to the terminals
			for(i=0;energy0>=0 && i<vectorlength*number_of_channels+taillength;i++)
			{
				totaloffset++;
				currentelement = (i == 0) ? 1.0 : 0.0;
				matrix_multiplication_iterationC(terms,conveyerwidth,patterns,pattern_size,kernel,kernel_size,conveyerlength,currentelement,rez1,pipelinedelay,0);
				//if(i>=conveyerlength)
				{
					extract_spectrumC(spectrum, terms, spectrumterms, spectrumtermtaps, original_matrix_width, conveyerlength, i);
					energy1 = calculate_energy(testspectrum, spectrum, original_matrix_width);
					if (energy1 || energy1*1.01 > energy0)
					{
						energy0 = -1.0;
					}
				}
			}
			//free buffers that are not used
			free2d(rez1,pattern_size);
		}
		free2d(terms,conveyerwidth);
	}
	return totaloffset-1;
}
