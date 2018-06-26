#ifndef __qtm__
#define __qtm__


/***************************************************************************
 *   Copyright (C) 2017 -- 2018 by Marek Sawerwain                         *
 *                                         <M.Sawerwain@gmail.com>         *
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
 
typedef simpleComplex<double> dblcmplx;

//extern uMatrix< simpleComplex<double>, 4> c_ops[ 1 ];
//extern template uVector< simpleComplex<double>, size_t > H;
//extern simpleComplex<double> alpha[2];
//extern uVector< simpleComplex<double>, 4> expect_operator;

template<size_t N, size_t Ntrj, size_t _WV_LEAD_DIM, size_t _WV_LEAD_DIM_SQR, size_t _C_OPS_SIZE>
int mpi_main(int argc, char *argv[], int verbose_mode,
						 double _from_time, 
						 double _to_time,
						 int use_colappse_operator, int use_expecation_operator);


 #endif // __qtm__