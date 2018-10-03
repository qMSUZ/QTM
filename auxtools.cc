/***************************************************************************
 *   Copyright (C) 2018 by Marek Sawerwain                         	       *
 *                                     <M.Sawerwain@gmail.com>             *
 *                                                                         *
 *   Part of the Quantum Trajectory Method:                                *
 *   https://github.com/qMSUZ/QTM                                          *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
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