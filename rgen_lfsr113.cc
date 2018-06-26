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
 
#include "rgen_lfsr113.h"

typedef struct {
  unsigned int z1,
			   z2,
			   z3,
			   z4;
} lfsr113_internal_state;


/**
 * Internal state for LFST113 pseudo-random number generator,
 * used as default state for PR generator.
 * 
 */
 
static lfsr113_internal_state lfsr113_ds; 


typedef union {
    float s;
    unsigned int u;
} _float_hlp;

typedef union {
    double d;
    unsigned int u[2];
} _double_hlp; 

/**
 * helper variables for generation of
 * floating point numbers from zero to one
 */

_float_hlp _flt_mantissa;
_double_hlp _dbl_mantissa;
 
static void prepare_mantissa() 
{
	
	volatile double x = 1.0; 
	
	_double_hlp t;
	_float_hlp s;
	
	double y = 0.5;
	do
	{			    
	    t.d = x;
	    x += y;
	    y *= 0.5;
	}
	while (x != t.d && x < 2.0);

	volatile float xx = 1.0;
                            
	float yy = 0.5;
	do
	{
		s.s = xx;
		xx += yy;
	    yy *= 0.5;
	}
	while (xx != s.s && xx < 2.0);	

	_dbl_mantissa.d = 1.0;
	_dbl_mantissa.u[0] ^= t.u[0];
	_dbl_mantissa.u[1] ^= t.u[1];

	_flt_mantissa.s = 1.0;
	_flt_mantissa.u ^= s.u;
	
}
 
/** 
 * Very important remark for lfsr113 !!!!
 * The initial seeds z1, z2, z3, z4  must be larger than
 * 1, 7, 15, and 127 respectively. 
 */
 
unsigned int lfsr113_init()
{
	lfsr113_ds.z1 = 2;
	lfsr113_ds.z2 = 8;
	lfsr113_ds.z3 = 16;
	lfsr113_ds.z4 = 128;
	
	prepare_mantissa();
	
	return 0;
}

unsigned int  lfsr113_generator_init(unsigned int z1, unsigned int z2, unsigned int z3, unsigned int z4)
{
  lfsr113_ds.z1 = z1;
  lfsr113_ds.z2 = z2;
  lfsr113_ds.z3 = z3;
  lfsr113_ds.z4 = z4;
  
  prepare_mantissa();
  
  return 0;
}

unsigned int lfsr113_genRand()
{
   unsigned long int b;

   b = ((lfsr113_ds.z1 <<  6) ^ lfsr113_ds.z1) >> 13;
   lfsr113_ds.z1 = ((lfsr113_ds.z1 & 4294967294UL) << 18) ^ b;
   
   b = ((lfsr113_ds.z2 <<  2) ^ lfsr113_ds.z2) >> 27;
   lfsr113_ds.z2 = ((lfsr113_ds.z2 & 4294967288UL) <<  2) ^ b;
   
   b = ((lfsr113_ds.z3 << 13) ^ lfsr113_ds.z3) >> 21;
   lfsr113_ds.z3 = ((lfsr113_ds.z3 & 4294967280UL) <<  7) ^ b;
   
   b = ((lfsr113_ds.z4 <<  3) ^ lfsr113_ds.z4) >> 12;
   lfsr113_ds.z4 = ((lfsr113_ds.z4 & 4294967168UL) << 13) ^ b;
   
   return (lfsr113_ds.z1 ^ lfsr113_ds.z2 ^ lfsr113_ds.z3 ^ lfsr113_ds.z4);
} 


float lfsr113_genRand_asFloat()
{
	_float_hlp _flt_rslt;
	
    _flt_rslt.s = 1.0;
    _flt_rslt.u |= (lfsr113_genRand() & _flt_mantissa.u);
    _flt_rslt.s -= 1.0;
	
	return _flt_rslt.s;

}

double lfsr113_genRand_asDouble()
{
	_double_hlp _dbl_rslt;

    _dbl_rslt.d = 1.0;
    _dbl_rslt.u[0] |= (lfsr113_genRand() & _dbl_mantissa.u[0]);
    _dbl_rslt.u[1] |= (lfsr113_genRand() & _dbl_mantissa.u[1]);
    _dbl_rslt.d -= 1.0;

    return  _dbl_rslt.d;
}

