#include "factorizator.h"

FILE *logfile = NULL;
//FILE *matrixfile = NULL;

bool make_idempotent(ptfr ptf)
{
	bool isok = true;
	if (ptf->timescale <= 0)
	{
		double **m = ptf->matrix;
		unsigned M = ptf->width;
		unsigned N = ptf->length;
		//check
		for (unsigned j = 0; j < N; j++)
		{
			bool notallzeros = false;
			bool allones = true;
			for (unsigned i = 0; i < M; i++)
			{
				if (m[i][j])
				{
					notallzeros = true;
					if (m[i][j] != 1)
					{
						allones = false;
						break;
					}
				}
			}
			isok = isok && allones && notallzeros;
			if (!isok)
				break;
		}
		//update
		if (isok)
		{
			for (unsigned j = 0; j < N; j++)
				for (unsigned i = 0; i < M; i++)
					if (m[i][j])
						m[i][j] = j+1;
		}
	}
	return isok;
}

void factorizator(ptfr tf)
{
	unsigned totaloffset,maxtap,pattern_size_max,pattern_size,kernel_size,kernel_size_max;
	int extention;
	unsigned ind0,ind1,initial_position,j;
	int **patterns,*spectrumterms,**commutator_image,*spectrumtermtaps,*patternoccurences;
	double *kernel,*spectrum,*originalspectrum;
	unsigned dropthreshold;
	kernel_size_max = (abs(tf->GF2)>0) ? (1<<(abs(tf->GF2)))-1 : tf->length*tf->width;
	dropthreshold = tf->dropthreshold;
	extention=0;
	pattern_size_max=tf->length*tf->width;
#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_matrix_d("input matrix", tf->matrix, tf->width, tf->length);
#endif
#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_matrix_d_for_matlab("matrix",tf->matrix, tf->width, tf->length);
#endif	
	spectrum=malloc1d(tf->width);
	originalspectrum=malloc1d(tf->width);
	/////matrix reduction
	kernel=malloc1d(kernel_size_max);
	commutator_image=malloc2i(tf->width, tf->length);
	spectrumterms=malloc1i(tf->width);
	spectrumtermtaps=malloc1i(tf->width);
	for (unsigned int i = 0;i < tf->width;i++)
	{
		spectrumterms[i] = -1;
		spectrumtermtaps[i] = -1;
	}
	patternoccurences=malloc1i(pattern_size_max);
	patterns=malloc2i(pattern_size_max,5);
	kernel_size=factorizeC(tf->matrix,kernel,kernel_size_max,commutator_image, tf->width, tf->length);
	if (kernel_size > kernel_size_max)
	{
		printf("Kernel too long, exceeding %d elements\n", kernel_size_max);
		//free memory
		free2d(tf->matrix, tf->width);
		free(spectrum);
		free(originalspectrum);
		free(kernel);
		free2i(commutator_image, tf->width);
		free(spectrumterms);
		free(spectrumtermtaps);
		free(patternoccurences);
		free2i(patterns, pattern_size_max);
		//exit
		return;
	}

#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_matrix_i("Commutator image",commutator_image, tf->width, tf->length);
	print_vector_d("Kernel",kernel,kernel_size);
#endif
	pattern_size=serialize_commutatorC(commutator_image, tf->width, tf->length, tf->pipelinedelay, tf->numberofchannels,patterns,pattern_size_max,patternoccurences,spectrumterms,spectrumtermtaps,&extention,dropthreshold);
	update_pattern_timescale(tf->timescale,patterns, pattern_size, spectrumtermtaps, tf->width);
#ifdef __SHOW_SYNTHESYS_PROCESS__
	printf("Extention = %d\n",extention);
	print_matrix_i("Commutator image after serialization",commutator_image, tf->width, tf->length);
#endif
#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_vector_i("extention", &extention, 1);
	print_matrix_i("Patterns",patterns,pattern_size,5);
	print_vector_i("Pattern occurences",patternoccurences,pattern_size);
	print_vector_i("Spectrum terms",spectrumterms, tf->width);
	print_vector_i("Spectrum term taps",spectrumtermtaps, tf->width);
#endif

/////\\\\==/\/\/\/\/\/\/\/\/\/\==
#ifdef __SHOW_SYNTHESYS_PROCESS__
	puts("Test: Will compare results of fast and conventional multiplication on the same arbitrary vector\n");
#endif
	double testvector[16]={1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    unsigned testvectorlength=16;

	//conventional mult.
	initial_position=0;
	matrix_multiplication_conventional_for_offset_search(tf->matrix, tf->width, tf->length,originalspectrum);
#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_vector_d("Spectrum of the conventional method",originalspectrum, tf->width);
#endif
	//fast multiplication.
	totaloffset=matrix_multiplication_find_offset(spectrum,patterns,pattern_size,spectrumterms,spectrumtermtaps, tf->width, tf->length,kernel,kernel_size,extention,testvector,testvectorlength,tf->pipelinedelay, tf->numberofchannels,originalspectrum);
#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_vector_d("Spectrum of the fast method",spectrum, tf->width);
#endif
	//compare
	double sumtotaldiff = 0.0;
	for (ind1 = 0;ind1 < tf->width;ind1++)
	{
		originalspectrum[ind1] -= spectrum[ind1];
		sumtotaldiff += fabs(originalspectrum[ind1]);
	}
#ifdef __SHOW_SYNTHESYS_PROCESS__
	print_vector_d("Difference between the fast and conventional method",originalspectrum, tf->width);
	printf("offset=%d\ntotaldifference=%g\n",totaloffset, sumtotaldiff);
#endif

/////\\\\==/\/\/\/\/\/\/\/\/\/\==	
	maxtap=0;
	for(j=0;j<tf->width;j++)
	{
		maxtap=maxtap<(unsigned)spectrumtermtaps[j]?(unsigned)spectrumtermtaps[j]:maxtap;
	}	
	//populate the result for synthesys.
	//original matrix dimensions
	//architecture
	//tf->pipelinedelay=pipelinedelay;
	//constraints
	tf->totaloffset=totaloffset;
	//internal parameters
	tf->maxtap=maxtap;
	tf->extention=extention;
	//kernel
	tf[0].kernel_size=kernel_size;
	tf[0].kernel=malloc1d(tf[0].kernel_size);
	if(tf[0].kernel)
		memcpy(tf[0].kernel,kernel,sizeof(double)*tf[0].kernel_size);
	free(kernel);
	//commutator
	tf[0].size= tf->width;
	tf[0].commutator_image=malloc2i(tf[0].size, tf->length);
	if(tf[0].commutator_image)
		for(ind0=0;ind0<tf[0].size;ind0++)
			if(tf[0].commutator_image[ind0])
				memcpy(tf[0].commutator_image[ind0],commutator_image[ind0], tf->length);
	free2i(commutator_image, tf->width);
	//patterns
	tf[0].number_of_patterns=pattern_size;
	tf[0].patterns=malloc2i(tf[0].number_of_patterns,5);
	if(tf[0].patterns)
		for(ind0=0;ind0<tf[0].number_of_patterns;ind0++)
			if(tf[0].patterns[ind0])
				memcpy(tf[0].patterns[ind0],patterns[ind0],sizeof(int)*5);
	free2i(patterns,pattern_size_max);
	tf[0].spectrumterms=spectrumterms;
	///free(spectrumterms);
	tf[0].spectrumtermtaps=spectrumtermtaps;
	////free(spectrumtermtaps);
	tf[0].patternoccurences=malloc1i(tf[0].number_of_patterns);
	if(tf[0].patternoccurences)
		memcpy(tf[0].patternoccurences,patternoccurences,sizeof(int)*tf[0].number_of_patterns);
	free(patternoccurences);

	////done
	free(spectrum);
	free(originalspectrum);
}
