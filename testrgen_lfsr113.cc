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
 
#include <iostream>

#include "rgen_lfsr113.h"

using namespace std;


int main(int argc, char *argv[])
{
	int i;

	lfsr113_generator_init( 0xA234, 0xB123, 0xC123, 0xD123);

	cout << "Ten prng values as float:" << endl;
	for(i=0;i<10;i++)
	{
		cout << "rng: " << lfsr113_genRand_asFloat() << endl;
	}

	cout << "Ten prng values as double:" << endl;
	for(i=0;i<10;i++)
	{
		cout << "rng: " << lfsr113_genRand_asDouble() << endl;
	}

	return 0;
}