//2018/01/20 - multiple formats, handles matrices with rows having a single element. 

/*How to set command line parameters in Visual Studio
On the main IDE window's menu bar choose Project -> {project name} Properties. That
opens the project's property pages. You can specify command line arguments under
Configuration Properties->Debugging->Command Arguments.Just type in exactly what
you would type at the command line after the program name.
If the command line must look like this:
synthesator.exe ..\..\..\step_3_synthesize_matrix_vector_multiplication_schematics_gl1\H_comp 0 1
put this in the command arguments field :
..\..\..\step_3_synthesize_matrix_vector_multiplication_schematics_gl1\H_comp 0 1


the calling format is:
synthesator.exe [<system name>='matrix' [<GF2>='0' [<pipelinedelay>='0' <timescale>='1']]]
here:
<system name> - char* path and name of a *.txt file with original matrix without .txt extention
<GF2> - int, 
	= 0 - arithmetic matrix (*,+)
	< 0 - geometric matrix (^(0,1),*) - needs hand corrections if kernel blocking is used
	> 0 - logical matrix (&,xor) - - needs hand corrections if kernel blocking is used
<piprlinedalay> - >=0, delays in sum or mult operators if GF2<=0
<timescale> - int,
	<= 0 - all balansing delays are 0, idempotent operation, inputs must be separated.
	> 0	- serial input, time depenency not distorted
*/

#include "stdafx.h"
#include "factorizator.h"
#include <malloc.h>
#include <iostream>
using namespace std;


class CnnFastLayer
{
public:

	int rc;
	tfr tf;
	int synthesised;

	int allocatematrix(tfr *ptf, unsigned precision, unsigned width, unsigned length)
	{
		int rc = -1;
		synthesised = 0;
		ptf->precision = precision;
		ptf->width = width;
		ptf->length = length;
		ptf->matrix = malloc2d(ptf->width, ptf->length);
		double preccoef = 0.5/RAND_MAX;
		for(int i=0;i<precision;i++) preccoef *= 10.0;

		if(ptf->matrix)
		{
			rc = 0;
			int v;
			for (unsigned i = 0;i < ptf->width;i++)
			{
				for (unsigned j = 0;j < ptf->length;j++)
				{
					//v = ((i+j)<(ptf->length/2) && (i+j)>=(ptf->length)) ? -1 : 1;
					v = (int)(((rand()<<2) - RAND_MAX) * preccoef);
					ptf->matrix[i][j] = v;
				}
			}
		}
		return rc;
	}


	CnnFastLayer(unsigned precision, unsigned sizex, unsigned sizey)
	{
		tf.spectrum = NULL;
		tf.vector = NULL;
		tf.terms = NULL;
		rc = allocatematrix(&tf, precision, sizey, sizex);
		if(rc==0)
		{
			tf.spectrum=malloc1d(tf.width);
			if (tf.spectrum)
			{
				tf.vector=malloc1d(tf.length);
			}
		}
		if((NULL == tf.spectrum) || (NULL == tf.vector)) 
		{
			rc |= 0x5;
		} 
	}

	~CnnFastLayer(void)
	{
		if(tf.spectrum)
		{
			std::free(tf.spectrum);
			tf.spectrum = NULL;
		}

		if(tf.vector)
		{
			std::free(tf.vector);
			tf.vector = NULL;
		}

		if(tf.terms)
		{
			free2d(tf.terms,tf.conveyerwidth);
			tf.terms = NULL;
		}

		cleanall();
	}


	void cleansynthesised(void)
	{
		if(synthesised)
		{
			free2i(tf.commutator_image,tf.size);
			free(tf.kernel);
			free(tf.patternoccurences);
			free2i(tf.patterns,tf.number_of_patterns);
			free(tf.spectrumterms);
			free(tf.spectrumtermtaps);
			
			free2d(tf.terms, tf.conveyerwidth);
			synthesised = 0;
			tf.truevectorindex = 0;
			tf.conveyerlength=0;//+tf.extention;
			tf.conveyerwidth=tf.kernel_size+tf.number_of_patterns;
			tf.conveyerindex = 0;
		}
	}




	int synthesise(void)
	{
		rc=0;
		unsigned i,j;
		unsigned ker_num;
		int **pat=NULL;
		double *kern=NULL;
		////p a r a m e t e r s ! ! !

		if(synthesised)
		{
			cleansynthesised();
		}

		if(tf.matrix && tf.length && tf.width) //matrix exists
		{
			tf.GF2 = 0; //0 -arithmetic, <0 - multiplication of terms, >0 - logical
			tf.pipelinedelay = 0;
			tf.timescale = 1;
			tf.numberofchannels = 1;
			tf.nummux = 1;
			tf.dropthreshold = 0;
			tf.offset = 1000000;
			if (make_idempotent(&tf))
			{
				factorizator(&tf);
#ifdef __SHOW_SYNTHESYS_PROCESS__
				print_matrix_i("patterns", tf.patterns, tf.number_of_patterns, 5);
#endif
#ifdef __SHOW_SYNTHESYS_PROCESS__
				print_vector_d("kernel", tf.kernel, tf.kernel_size);
#endif
				//find delays necessary for taps to output ports.
				int minval=tf.spectrumtermtaps[0];
				for(int i=0;i<tf.size;i++)
				{
					if(minval< tf.spectrumtermtaps[i]) minval = tf.spectrumtermtaps[i];
				}
				//adjust delays necessary for taps to output ports.
				for(int i=0;i<tf.size;i++)
				{
					tf.spectrumtermtaps[i]=minval-tf.spectrumtermtaps[i];
				}
#ifdef __SHOW_SYNTHESYS_PROCESS__
				printf("totaloffset=%d\n", tf.totaloffset);
				printf("length=%d\n", tf.length);
				printf("numberofchannels=%d\n", tf.numberofchannels);
				printf("pipelinedelay=%d\n", tf.pipelinedelay);
				printf("maxtap=%d\n", tf.maxtap);
				printf("extention=%d\n", tf.extention);
#endif
			}
			else
			{
				puts("wrong time scale for idempotent operation\n");
				rc = -2;
			}
		}
		else
		{
			puts("Matrix not present\n");
			rc = -1;
		}
		synthesised = (rc==0) ? 1 : 0;

		if(synthesised)
		{
			tf.truevectorindex = tf.length-1;
			tf.conveyerlength=tf.length;//+tf.extention;
			tf.conveyerwidth=tf.kernel_size+tf.number_of_patterns;
			tf.terms=malloc2d(tf.conveyerwidth,tf.conveyerlength);
			tf.conveyerindex = 0;
			for(int i=0; i<tf.length; i++) tf.vector[i] = 0;
		}

		return rc;
	}



	void cleanall(void)
	{
		cleanup_tfr(&tf);
	}

	int generate_and_test(unsigned long *prods_fast, unsigned long *adds_fast, unsigned long *prods_conv, unsigned long *adds_conv)
	{
		int numberofiterations = 40; //tf.length; //tf.extention;
		int ind;
		double sum = 0;
		double *buffer = (double*) malloc(sizeof(double)*tf.width);
		if (buffer)
		{ 
			rc = synthesise();
			*prods_fast += tf.kernel_size * numberofiterations;
			*adds_fast += (tf.patterns[tf.number_of_patterns-1][0]-tf.patterns[0][0])*numberofiterations;
			*prods_conv += tf.length*tf.width*numberofiterations;
			*adds_conv += (tf.length-1)*tf.width*numberofiterations;
			///allocate memory for conventional multiplier
			for(int i=0; (i<numberofiterations) && (sum == 0); i++)
			{
				matrix_multiplication_conventional_step(tf.matrix,tf.width,tf.length,tf.spectrum,tf.vector,&tf.truevectorindex,tf.GF2,i+5);
				print_vector_d_offset("Vector of the conventional method",tf.vector,tf.length,tf.truevectorindex);
				//print_vector_d("Spectrum of the conventional method",tf.spectrum,tf.width);
				ind = 1; for(int x=0; x<=tf.width; x++)
				{
					buffer[x]=tf.spectrum[x];
				}
				matrix_multiplication_iteration(tf.terms,tf.conveyerwidth,tf.patterns,tf.number_of_patterns,tf.kernel,tf.kernel_size,tf.conveyerlength,i+5,tf.conveyerindex);
				extract_spectrumC(tf.spectrum,tf.terms,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.conveyerlength,tf.conveyerindex);
				ind = 0;
				sum = 0.0;
				for(int x=0; x<=tf.width; x++) 
				{ 
					sum+=fabs(tf.spectrum[x]-buffer[x]);
					if(sum>1.0)
					{
						sum = 0.0;
						ind++;
					}
					//tf.spectrum[x] -= buffer[x];
				}
				//print_vector_d("Spectrum difference beween fast and conventional methods:",tf.spectrum,tf.width);
				print_vector_d_difference("Spectrum conv method: red - above, green - below fast.", buffer, tf.spectrum, 0.5, tf.width);
				print_vector_d_difference("Spectrum fast method: red - above, green - below conv.", tf.spectrum, buffer, 0.5, tf.width);
				//printf("conveyerindex = %d:",tf.conveyerindex);
				tf.conveyerindex = (tf.conveyerindex + 1) % tf.conveyerlength;
				//print_matrix_d("Terms",tf.terms,tf.conveyerwidth,tf.conveyerlength);
				printf("sum = %g, index = %d\n", sum, ind);
			}
			//if(sum > 1.0) rc |= 0x7;
		}
		else
		{
			rc = 0x9;
		}
		return rc;
	}


	int mr(unsigned x, unsigned y)
	{
		int rv = 0;
		if(tf.matrix && tf.width > y && tf.length > x)
		{
			rv = tf.matrix[y][x];
		}
		return rv;
	}

	int mw(int v, unsigned x, unsigned y)
	{
		int rv = 0;
		if(tf.matrix && tf.width > y && tf.length > x)
		{
			tf.matrix[y][x] = v;
		}
		else
		{
			rv = -1;
		}
		return rv;
	}

	int mx(void)
	{
		int rv = 0;
		if(tf.matrix && tf.width && tf.length)
		{
			rv = tf.length;
		}
	}

	int my(void)
	{
		int rv = 0;
		if(tf.matrix && tf.width && tf.length)
		{
			rv = tf.width;
		}
	}

};




int usage_example_and_test(void)
{
	int precision = 4, sizex = 8,sizey =1024;
	CnnFastLayer *cnnfl =  NULL;
						//	new CnnFastLayer(precision, sizex, sizey);
	int rc = 0;
	int x,y,u = 0;
	unsigned long prods_fast = 0;
	unsigned long adds_fast = 0;
	unsigned long prods_conv = 0;
	unsigned long adds_conv = 0;
	for(x = 0; x < sizex && rc==0; x++)
	{
		for(y = 0; y < sizey && rc==0; y++)
		{
			cnnfl = new CnnFastLayer(precision, sizex, sizey);
			if(cnnfl)
			{
				//u = (u+y+x+cnnfl->mr(x,y))%10; 
				//cnnfl->mw(u,x,y);
				rc |= cnnfl->generate_and_test( &prods_fast, &adds_fast, &prods_conv, &adds_conv);
				delete cnnfl;
				cout << "x,y = " << x << ", " << y << endl;
				cout << "prods_fast " << prods_fast << ", prods_conv " << prods_conv << ", adds_fast " << adds_fast << ", adds_conv " << adds_conv << ", ratios: prod " << ((double)prods_conv/(double)prods_fast) << ", add " << ((double)adds_conv/(double)adds_fast) << endl;
			}
			else
			{
				cout << "last x,y = " << x << ", " << y << endl;
				for(;;);
			}
		}
	}
	//delete cnnfl;
	cout << "x,y = " << x << ", " << y << endl;
	return rc;
}


int main(int argc, char *argv[])
{
	return usage_example_and_test();
}