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

#include <iostream>
using namespace std;

int rc;
tfr tf;

int synthesise(int argc, char* argv[])
{
	rc=0;
	unsigned i,j;
	char system_name[256];
	//if (argc > 1)
	//	strcpy(system_name, (char*)argv[1]);
	//else
		strcpy(system_name, "matrix");
	char orig_name[256], file_name[256], log_name[256], terms_file_name[256];
	strcpy(orig_name, system_name);
	strcpy(file_name, system_name);
	strcpy(log_name, system_name);
	strcpy(terms_file_name, system_name);
	strcat(orig_name, ".txt");
	strcat(file_name, ".mdl");
	strcat(log_name, "_log.txt");
	strcat(terms_file_name, "_terms.txt");
	unsigned ker_num;
	int **pat=NULL;
	double *kern=NULL;
	//double *matrix=NULL;
	//tfr tf;
	////p a r a m e t e r s ! ! !
	tf.precision = 2;
	tf.width = 3;
	tf.length = 4;
	tf.matrix = malloc2d(tf.width, tf.length);
	if (tf.matrix)
	{
		int v = 1;
		for (unsigned i = 0;i < tf.width;i++)
		{
			for (unsigned j = 0;j < tf.length;j++)
			{
				tf.matrix[i][j] = v++;
			}
		}
	}
	tf.GF2 = 0; //0 -arithmetic, <0 - multiplication of terms, >0 - logical
	tf.pipelinedelay = 0;
	tf.timescale = 1;
	tf.numberofchannels = 1;
	tf.nummux = 1;
	tf.dropthreshold = 0;
	tf.offset = 1000000;
	printf("opening log file %s\n", log_name);
	logfile = fopen(log_name,"w");
	if (logfile)
	{
		//printf("opening matrix file %s\n", orig_name);
		//if (0 > (rc = getmatrix(&tf, orig_name)))
		//{
		//	printf("error %d\n", rc);
		//	fclose(logfile);
		//	return rc;
		//}

		if (make_idempotent(&tf))
		{
			factorizator(&tf);
#if (0)
			print_matrix_i("patterns", tf.patterns, tf.number_of_patterns, 5);
#endif
#if (0)
			print_vector_d("kernel", tf.kernel, tf.kernel_size);
#endif
			rc = extract((char *)file_name, tf);
			printf("totaloffset=%d\n", tf.totaloffset);
			printf("length=%d\n", tf.length);
			printf("numberofchannels=%d\n", tf.numberofchannels);
			printf("pipelinedelay=%d\n", tf.pipelinedelay);
			printf("maxtap=%d\n", tf.maxtap);
			printf("extention=%d\n", tf.extention);

			save_tfr(&tf, terms_file_name);
			//tfr tfn;
			//load_tfr(&tfn, terms_file_name);
			//cleanup_tfr(&tfn);
		}
		else
			puts("wrong time scale for idempotent operation\n");

		//cleanup_tfr(&tf);
		fclose(logfile);
	}
	else
		puts("couldn't open logfile\n");

#if (0)
	puts("===complete===, hit any key\n");
	while(!_kbhit());
#else
	puts("===complete===\n");
#endif
	return rc;
}


void fastmatrixproductC(double *spectrum, double *testvector1,unsigned testvectorlength1, int numberofiterations)
{
        matrix_multiplicationConveyorC(spectrum,tf.patterns,tf.number_of_patterns,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.length,tf.kernel,tf.kernel_size,tf.extention,testvector1,testvectorlength1,tf.pipelinedelay,tf.GF2,numberofiterations);
        print_vector_d("Spectrum of the fast method for long test vector",spectrum, tf.width);
        matrix_multiplicationOneVector(spectrum,tf.patterns,tf.number_of_patterns,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.length,tf.kernel,tf.kernel_size,tf.extention,testvector1,testvectorlength1,tf.pipelinedelay,tf.numberofchannels,tf.GF2,numberofiterations);
        print_vector_d("Spectrum of the fast 2 method for long test vector",spectrum, tf.width);
}

void fastmatrixproduct(double *spectrum, double *testvector1,unsigned testvectorlength1, int numberofiterations)
{
	int offset = numberofiterations - tf.length;
	int numiter = numberofiterations + tf.extention;
	print_matrix_d("Matrix",tf.matrix,tf.width,tf.length);
	print_vector_d("Vector",testvector1,testvectorlength1);
	matrix_multiplication_conventional(tf.matrix,tf.width,tf.length,spectrum,testvector1,testvectorlength1,offset,tf.GF2);
print_vector_d("Spectrum of the conventional method",spectrum, tf.width);
        matrix_multiplicationConveyorCI(spectrum,tf.patterns,tf.number_of_patterns,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.length,tf.kernel,tf.kernel_size,tf.extention,testvector1,testvectorlength1,tf.pipelinedelay,tf.GF2,numiter);
        print_vector_d("Spectrum of the fast method with CI",spectrum, tf.width);
        matrix_multiplicationConveyor(spectrum,tf.patterns,tf.number_of_patterns,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.length,tf.kernel,tf.kernel_size,tf.extention,testvector1,testvectorlength1,numiter);
        print_vector_d("Spectrum of the fast method vector",spectrum, tf.width);
        matrix_multiplicationOneVector(spectrum,tf.patterns,tf.number_of_patterns,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.length,tf.kernel,tf.kernel_size,tf.extention,testvector1,testvectorlength1,tf.pipelinedelay,tf.numberofchannels,tf.GF2,numiter);
        print_vector_d("Spectrum of the fast 2 method",spectrum, tf.width);
}


void cleanall(void)
{
	cleanup_tfr(&tf);
}

int main(int argc, char *argv[])
{
	rc = synthesise(argc, argv);
        ///allocate memory for iconventional multiplier
        tf.spectrum=malloc1d(tf.width);
        tf.vector=malloc1d(tf.length);
	
	tf.truevectorindex = tf.length-1;
	tf.conveyerlength=tf.length+tf.extention;
        tf.conveyerwidth=tf.kernel_size+tf.number_of_patterns;
        tf.terms=malloc2d(tf.conveyerwidth,tf.conveyerlength);
	tf.conveyerindex = 0;
	
	int numberofiterations = 50; //tf.length; //tf.extention;
	//fastmatrixproduct(tf.spectrum, testvector1, testvectorlength1, numberofiterations);
	for(int i=0; i<numberofiterations; i++)
	{
		matrix_multiplication_conventional_step(tf.matrix,tf.width,tf.length,tf.spectrum,tf.vector,&tf.truevectorindex,tf.GF2,i+5);
		print_vector_d_offset("Vector of the conventional method",tf.vector,tf.length,tf.truevectorindex);
		print_vector_d("Spectrum of the conventional method",tf.spectrum,tf.width);
		//matrix_multiplication_iteration(tf.terms,tf.conveyerwidth,tf.patterns,tf.kernel,tf.kernel_size,tf.conveyerlength,i+5,tf.conveyerindex);
		matrix_multiplication_iteration(tf.terms,tf.conveyerwidth,tf.patterns,tf.number_of_patterns,tf.kernel,tf.kernel_size,tf.conveyerlength,i+5,tf.conveyerindex);
                extract_spectrumC(tf.spectrum,tf.terms,tf.spectrumterms,tf.spectrumtermtaps,tf.width,tf.conveyerlength,tf.conveyerindex);
                print_vector_d("Spectrum, fast method:",tf.spectrum,tf.width);
		printf("%d:",tf.conveyerindex);
		tf.conveyerindex = (tf.conveyerindex + 1) % tf.conveyerlength;
		print_matrix_d("Terms",tf.terms,tf.conveyerwidth,tf.conveyerlength);
		
	}

        free(tf.spectrum);
        free(tf.vector);
        free2d(tf.terms,tf.conveyerwidth);
	cleanall();
	return rc;
}


