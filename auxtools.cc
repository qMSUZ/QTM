#include "auxtools.h"


template<typename T, size_t SIZE>
void dump_table_to_file(simpleComplex<T> *tbl, const char *filename)
{
	size_t i;
	FILE *foutput;
	
	foutput = fopen(filename, "w");
	
	for( i = 0 ; i < SIZE ; i++)
	{
		fprintf(foutput, "alpha[%d]=make_simpleComplex(%-16.16lf, %-16.16lf);\n", i, tbl[i].re, tbl[i].im);
		
	}
	fprintf(foutput, "\n");

	fclose(foutput);
	
}


template<typename T, size_t SIZE>
void dump_table_to_file(simpleComplex<T> *tbl, const char *matname, const char *filename)
{
	size_t i;
	FILE *foutput;
	
	foutput = fopen(filename, "w");
	
	for( i = 0 ; i < SIZE ; i++)
	{
		fprintf(foutput, "%s[%d]=make_simpleComplex(%-16.16lf, %-16.16lf);\n", matname, i, tbl[i].re, tbl[i].im);
		
	}
	fprintf(foutput, "\n");

	fclose(foutput);
	
}


template<typename T>
void dump_uCSRMatrix_to_file(const uCSRMatrix< simpleComplex<T> > &tmat, const char *matname, const char *filename)
{
	unsigned int i;
	FILE *foutput;
	
	foutput = fopen(filename, "w");
	
	for( i = 0 ; i < tmat._values_size ; i++)
	{
		fprintf(foutput, "%s.values[%d]=make_simpleComplex(%-16.16lf,%-16.16lf);\n", matname, i, tmat.values[i].re, tmat.values[i].im);
	}
	fprintf(foutput, "\n");

	for( i = 0 ; i < tmat._col_ind ; i++)
	{
		fprintf(foutput, "%s.col_ind[%d]=%d;\n", matname, i, tmat.col_ind[i]);
	}
	fprintf(foutput, "\n");

	for( i = 0 ; i < tmat._row_ptr ; i++)
	{
		fprintf(foutput, "%s.row_ptr[%d]=%d;\n", matname, i, tmat.row_ptr[i]);
	}
	fprintf(foutput, "\n");
	
	
	fclose(foutput);
	
}