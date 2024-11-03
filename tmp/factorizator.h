#ifndef __FACTORIZATOR__
#define __FACTORIZATOR__

#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <memory.h>
#include <math.h>

#define PI 3.141592653589793116
//#define __AUTONOMOUS_INDEXING__
#define __SHOW_RESULTS_FOR_EACH_SHIFT__


//#define __USE_OMP__

#ifdef __USE_OMP__
#include <omp.h>
#endif

typedef struct _destination
{
	char name;
	int index;
	int input;
	struct _destination *next;
} destination;

typedef struct _source
{
	char name;
	int index;
	int output;
} source;

typedef struct _node
{
	char name;
	int index[2];
	unsigned long numdel;
	source *src;
	destination *dst;
	struct _node *next;
} node;

typedef struct _trf
{
	// for synthesis
	int GF2;
	int timescale;
	//original matrix
	double **matrix;
	unsigned length; //Xdim
	unsigned width; //Ydim
	double precision;
	//architecture
	unsigned numberofchannels;
	unsigned pipelinedelay;
	unsigned nummux;
	//constraints
	unsigned totaloffset;
	//internal parameters
	unsigned maxtap;
	unsigned extention;
	//synthesis components
	double *kernel;
	unsigned kernel_size;
	int **commutator_image;
	unsigned size;
	int **patterns;
	unsigned number_of_patterns;
	int *spectrumterms;
	int *spectrumtermtaps;
	int *patternoccurences;
	unsigned dropthreshold;
	unsigned offset;
	//for fast matrix product
        unsigned conveyerlength;
	unsigned conveyerwidth;
	unsigned maxindex;
	unsigned conveyerindex;
        double **terms;
        double **rez1;
	double *spectrum;
	double *vector;
	int truevectorindex;
}	tfr, *ptfr;

extern FILE *logfile;

int getmatrix(tfr *ptf, char *matrixfilename);

int extract(char *model_name,tfr tf);
void factorizator(ptfr ptf);
bool make_idempotent(ptfr ptf);
void update_pattern_timescale(int timescale, int ** patterns, unsigned pattern_size, int *spectrumtermtaps, unsigned width);
void cleanup_tfr(ptfr ptf);
extern FILE *logfile;
extern FILE *matrixfile;
int *malloc1i(unsigned int size);
double *malloc1d(unsigned int size);
int **malloc2i(unsigned int rows,unsigned int cols);
void free2i(int **arrayPtr,unsigned rows);
double **malloc2d(unsigned rows,unsigned cols);
void free2d(double **arrayPtr,unsigned rows);
void print_matrix_i(char *name,int **matr,unsigned hights,unsigned length);
void print_matrix_d(char *name,double **matr,unsigned hights,unsigned length);
void print_matrix_d_for_matlab(char *name,double **matr,unsigned hights,unsigned length);
void print_vector_d(char *name,double *vect,unsigned length);
void print_vector_d_offset(char *name,double *vect,unsigned length,unsigned offset);
void print_vector_i(char *name,int *vect,unsigned length);
int powi(int arg,int pow);
unsigned powu(unsigned arg,unsigned pow);
//void form_sceleton_matrixC(int *alphabet,unsigned alphabet_size,int **sceleton,unsigned word_length,unsigned number_of_words);
//int upsample(int **sceleton,double **sceleton_with_context,unsigned word_length,unsigned number_of_words,unsigned nsamp);
//int calculate_response(double *rrcfilter,unsigned nsamp,unsigned filtorder,double rolloff);
//double apply_response(unsigned i,double *s,unsigned sn,double *r,unsigned rn);
//void accumulate_historyC(double **sceleton,unsigned s_width,unsigned s_length);
unsigned factorizeC(double **matr,double *kernel,unsigned kernel_size_max,int **commutator_image,unsigned m_width,unsigned m_length);
unsigned countpatternoccurence(int pattern1,int pattern2,int distancebetweenpatternelements,int **matrix,unsigned numberofvectors,unsigned vectorlength);
unsigned counttodesiredpatternoccurence(int pattern1,int pattern2,int distancebetweenpatternelements,int **matrix,unsigned numberofvectors,unsigned vectorlength,unsigned desired_count,int *tmp);
unsigned getmostcommonpatternC(int **matrix,unsigned numberofvectors,unsigned vectorlength,int *mcp);
unsigned getnextpatternC(int **matrix,unsigned numberofvectors,unsigned vectorlength,int *mcp,unsigned desired_repeat_num,int *tmp);
unsigned replacepatternC(int *pattern,int **matrix,unsigned numberofvectors,unsigned vectorlength);
unsigned reducetermsC(int **matrix,unsigned numberofvectors,unsigned vectorlength,int **patterns,unsigned pattern_size_max,int *occurences, unsigned dropthreshold);
int makemultichannelpatterns(int **patterns,unsigned pattern_size,unsigned number_of_channels);
int retardarguments(int **patterns,unsigned pattern_size,unsigned pipelinedelay,int *time_aligning_delays);
int calculatetotaldelay(int pattern,int **patterns,unsigned pattern_size,unsigned pipelinedelay,int *time_aligning_delays);
unsigned serialize_commutatorC(int **commutator_image,unsigned numberofvectors,unsigned vectorlength,unsigned pipelinedelay,unsigned number_of_channels,int **patterns,unsigned pattern_size_max,int *patternoccurences,int *spectrumterms,int *spectrumtermtaps,int *extention,unsigned dropthreshold);
void matrix_multiplication_iterationC(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement, int GF2);
void matrix_multiplication_iterationC(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement,int GF2,int conveyerindex);
void matrix_multiplication_iteration(double **terms,unsigned conveyerwidth,int **patterns,unsigned pattern_size,double *kernel,unsigned kernel_size,unsigned conveyerlength,double currentelement,int conveyerindex);
void extract_spectrumC(double *spectrum,double **terms,int *spectrumterms,int *spectrumtermtaps,unsigned numberofvectors,unsigned conveyerlength,int conveyerindex);
void matrix_multiplicationConveyorC(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,int GF2,int numberofiterations);
void matrix_multiplicationConveyorCI(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,int GF2,int numberofiterations);
void matrix_multiplicationConveyor(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,int numberofiterations);
void matrix_multiplicationOneVector(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,unsigned number_of_channels,int GF2,int numberofiterations);
double calculate_energy(double *vect0,double *vect1, unsigned len);
unsigned matrix_multiplication_find_offset(double *spectrum,int **patterns,unsigned pattern_size,int *spectrumterms,int *spectrumtermtaps,unsigned original_matrix_width,unsigned original_matrix_length,double *kernel,unsigned kernel_size,unsigned extention,double *vector,unsigned vectorlength,unsigned pipelinedelay,unsigned number_of_channels,double *testspectrum);
int matrix_multiplication_conventional_for_offset_search(double **matrix,unsigned matrix_width,unsigned matrix_length,double *originalspectrum);
int matrix_multiplication_conventional(double **matrix,unsigned matrix_width,unsigned matrix_length,double *originalspectrum,double *testvector,unsigned testvectorlength,int offset,int GF2);
int matrix_multiplication_conventional_step(double **matrix,unsigned matrix_width,unsigned matrix_length,double *originalspectrum,double *testvector,int *offset, int GF2, int currentvalue);

void print_source(source *s);
void print_destination(destination *d);
void print_destination_list(destination *l);
void print_node(node *n);
void print_node_list(node *l);
node *find_node(node *l,char name,int index0,int index1);
node *add_node(node *l,char name,int index0,int index1);
int add_source(node *l,char name,int index0,int index1,char src_name,int src_ind, int src_out);
int add_destination(node *l,char name,int index0,int index1,char dst_name,int dst_ind,int dst_inp);
void free_destination_list(destination *l);
void free_node_list(node *n);

FILE *create_subsystem(char* file_name, int input_ports_number,int output_ports_number);
int place_inport(FILE *file_handle, int port_number);
int place_reductor(FILE *file_handle, unsigned red_number, int GF2);
int place_gain(FILE *file_handle, unsigned number, double gain_r, double gain_i, int GF2);
int place_outport(FILE *file_handle, int port_number);
int close_subsystem(FILE *file_handle);
int place_delay(FILE *file_handle, unsigned delay_number,unsigned delay);
int place_unit_delay(FILE *file_handle, unsigned delay_number);
int place_sum(FILE *file_handle, unsigned sum_number, int GF2);
int place_sum_with_delay(FILE *file_handle,unsigned subsystem_number,unsigned pipeline_delay, int GF2);
int place_mux(FILE *file_handle, unsigned mux_number,unsigned num_ports);
int place_source(FILE *file_handle,char name,unsigned nam_num,int port);
int place_destination(FILE *file_handle,char name,unsigned nam_num,int port);
int end_destinations(FILE *file_handle);
void place_node(FILE *file_handle,node *n);
void place_node_list(FILE *file_handle,node *n);

int load_tfr(ptfr tf, const char *filename);
int save_tfr(ptfr tf, const char *filename);

#endif
