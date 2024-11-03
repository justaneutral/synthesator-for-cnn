#include <malloc.h>
typedef float flpt;
typedef unsigned uint;

int inmemory_matrix_product(flpt *spectrum, uint *patterns, uint *spectrumterms, uint *spectrumtermtaps, uint *commutator_image, flpt *kernel, uint extention, flpt *vector, uint vectorlength, uint matrixwidth, uint numpatterns)
{
	int rc = 0;
	int tindex = 0;
	int co;
	uint conveyorwidth, conveyorlength;
	flpt *terms = NULL;
	uint *staps = (uint*)malloc(sizeof(int)*matrixwidth);
	initialize_terms(terms, &conveyorwidth, &conveyorlength, matrixwidth, vectorlength, patterns, extention, numpatterns, &co);
	for (int vectorindex = 1;vectorindex <= vectorlength;vectorindex++)
	{
		flpt currentelement = vector[vectorindex - 1];
		matrix_multiplication_iteration(terms, &tindex, patterns, kernel, conveyorwidth, conveyorlength, currentelement, numpatterns, &co);
	}
	for (int spectrumindex = 0;spectrumindex < matrixwidth;spectrumindex++)
		staps[spectrumindex] = 1 + (tindex + spectrumtermtaps[spectrumindex]) % conveyorlength;
	extract_spectrum(spectrum, terms, spectrumterms, staps, conveyorwidth, conveyorlength, matrixwidth);
	if (staps)
		free(staps);
	if (terms)
		free(terms);
	return rc;
}

int initialize_terms(flpt *terms, uint *conveyerwidth, uint *conveyerlength, uint numberofvectors, uint vectorlength, uint *patterns, uint extention, uint psv, int *co)
{
	*conveyerlength = vectorlength + extention;
	*co = *conveyerlength;
	if (psv > 0)
		*conveyerwidth = patterns[psv * 5];
	else
		*conveyerwidth = numberofvectors;

	terms = (flpt*)malloc(sizeof(flpt)**conveyerwidth**conveyerlength);
	for (int i = 0;i < *conveyerwidth**conveyerlength;i++)
		terms[i] = 0;
}

int matrix_multiplication_iteration(flpt *terms, int *tindex, uint *patterns, flpt *kernel, uint conveyerwidth, uint conveyerlength, flpt currentelement, uint psv, int *co)
{
	//propagate main conveyer(by changing pointer to its load position)
	*co = 1 + *co % conveyerlength;
	//populate products
	for (int kernind = 0;kernind < patterns[0] - 1;kernind++)
		terms[kernind*conveyerwidth + *co - 1] = currentelement*kernel[kernind];
	//here we update terms one by one.
	if (psv > 0)
	{
		for (int termsrow = patterns[0]; termsrow <= patterns[psv * 5];termsrow += 5)
		{
			uint patternsrow = termsrow - patterns[0];
			uint p1 = patterns[patternsrow * 5 + 1];
			uint p1delay = patterns[patternsrow * 5 + 3];
			uint p1pos = 1 + (*co - p1delay - 1) % conveyerlength;
			flpt arg1 = terms[p1*conveyerwidth + p1pos - 1];
			uint p2 = patterns[patternsrow * 5 + 2];
			uint p2delay = patterns[patternsrow * 5 + 4];
			uint p2pos = 1 + (*co - p2delay - 1) % conveyerlength;
			flpt arg2 = terms[p2*conveyerwidth + p2pos - 1);
			flpt rez = arg1 + arg2;
			terms[termsrow*conveyerwidth + *co - 1] = rez;
		}
	}
	*tindex = *co%conveyerlength;
	return *tindex;
}

int extract_spectrum(flpt *spectrum, flpt *terms, uint *spectrumterms, uint *spectrumtermtaps, uint conveyorwidth, uint conveyorlength, uint numberofvectors)
{
	spectrum = (flpt*)malloc(numberofvectors * sizeof(flpt));
	for (int spectrumindex = 0; spectrumindex < numberofvectors;spectrumindex++)
		if (spectrumterms[spectrumindex] > 0)
			spectrum[spectrumindex] = terms[(spectrumterms[spectrumindex] - 1)*conveyorwidth + spectrumtermtaps[spectrumindex]];
}
