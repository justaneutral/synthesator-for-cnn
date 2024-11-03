#include "factorizator.h"

int getmatrix(tfr *ptf, char *matrixfilename)
{
	int rc = -7;
	unsigned long totcnt = 0, tot = 0;
	long numberofBlocks, numberofBits, blockLength;
	int reversed = 0;

	FILE *matrixfile = NULL;
	matrixfile = fopen(matrixfilename, "r");
	if (matrixfile)
	{
		rc++;
		rc += fscanf(matrixfile, "%lf",&(ptf->precision)); //precision or whatever it means here

		rc += fscanf(matrixfile, "%d", &(ptf->width)); //vertical size - output vector length

		rc += fscanf(matrixfile, "%d", &(numberofBits));  // horizontal size - input vector length 
		rc += fscanf(matrixfile, "%d", &(blockLength)); // horizontal block length
		if(blockLength == 0)
		{ 
			blockLength = 1;
		}
		else
			if (blockLength < 0)
			{
				blockLength = abs(blockLength);
				reversed = 1;
			}
		numberofBlocks = (numberofBits + blockLength - 1) / blockLength;
		ptf->length = numberofBlocks;
		printf("matrix structire: Y x X = %ld x %ld\n", ptf->width, numberofBits);
		printf("matrix line blocking structire: Length: %ld, Blocks: %ld, Block length: %ld, Last Block Length: %ld \n", numberofBits, ptf->length, blockLength, (numberofBits%blockLength==0)?(numberofBits%blockLength):blockLength);

		//ptf->matrix = (double *)malloc(sizeof(double)*ptf->length*ptf->width);
		ptf->matrix = (rc > -3) ? malloc2d(ptf->width, ptf->length) : NULL;
		if (ptf->matrix)
		{
			rc++;
			int v;
			for (unsigned i = 0;i < ptf->width;i++)
			{
				for (unsigned j = 0;j < ptf->length;j++)
				{
					tot++;
					totcnt += fscanf(matrixfile, "%d", &v);
					if(totcnt==tot)
						//ptf->matrix[i*ptf->length + j] = v;
						ptf->matrix[i][j] = v;
					else
					{
						fclose(matrixfile);
						//free(ptf->matrix);
						free2d(ptf->matrix, ptf->width);
						ptf->matrix = NULL;
						ptf->precision = 0;
						ptf->length = 0;
						ptf->width = 0;
						return rc;
					}
				}
			}
		}
		fclose(matrixfile);
	}
	return ++rc;
}
