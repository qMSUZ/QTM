#ifndef __QTM__
#define __QTM__

/***************************************************************************
 *   Copyright (C) 2017 -- 2018 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
 *							   and Joanna Wiśniewska                       *
 *                                         <jwisniewska@wat.edu.pl>        *
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
 
#include "complexnum.h"
 
#define OUTPUT_FILE 1000 
#define OUTPUT_FILE_PYTHON_STYLE 1001
#define OUTPUT_STATE_OF_TRJ_FILE 1500 
#define METADAMS 	2000
#define METBDF		3000

#define COLLAPSE_OPERATOR 4000
#define EXPECATION_OPERATOR 4001
#define COLLAPSE_FUNCTION 4002
#define EXPECATION_FUNCTION 4003

#define TRJ_PERFECT 80
#define TRJ_GOOD_WITH_COLLAPSE 71
#define TRJ_BAD 66
#define TRJ_UNKNOWN 85
 
typedef simpleComplex<double> dblcmplx;


typedef struct {
	int type_output;
	int state_of_trj_output;
	int only_final_trj;
	int rnd_test_retry;
	int ode_method;
	int verbose_mode;
	double tolerance;
	char *file_name;
	int (*fnc)(long int *NEQ, double *T, dblcmplx *Y, dblcmplx *YDOT, dblcmplx *RPAR, long int *IPAR);
	
} extra_options;


template<size_t N, size_t Ntrj, size_t _WV_LEAD_DIM, size_t _WV_LEAD_DIM_SQR, size_t _C_OPS_SIZE>
int mpi_main(int argc, char *argv[],
						 double _from_time, 
						 double _to_time,
						 int use_colappse_operator, int use_expecation_operator,
						 extra_options &opt);


 #endif // __QTM__