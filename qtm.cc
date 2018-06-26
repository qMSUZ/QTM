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
#include <cstdio>


#include <mpi.h>

#include "complexnum.h"

#ifdef _WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#endif

#ifndef _WIN32
#include <sys/socket.h>
#include <sys/ioctl.h>
#include <net/if.h>
#include <unistd.h>
#include <arpa/inet.h> 
#endif


#include "complexnum.h"
#include "rgen_lfsr113.h"
#include "qtm.h"


#define MASTER_NODE 0

#define cmd_NOP 0x0000

#define cmd_RNG_TEST 0x1000
#define cmd_TRJ_PROC 0x1001
#define cmd_COMPUTED_TRJ 0x2000
#define cmd_EXIT 0x1fff


#define CMD_DATA_TAG 0xA000
#define TRJ_DATA_TAG 0xA001

#define GLOBAL_CODE 0xAFFFFF


long int mf;
long int neq;
dblcmplx err;
long int liw, lrw, lzw, ipar;
double atol_val;
dblcmplx rpar;
long int itol, iopt;
double rtol;
long int iout;
double tout;
dblcmplx wtru;
double aberr, aemax;
long int itask;
long int saveitask;
long int iwork[30];
double dtout;
double rwork[22];
dblcmplx zwork[30];
long int istate;




void print_ip(const int host_rank, const char *processor_name)
{
    const char* google_dns_server = "8.8.8.8";
    int dns_port = 53;

    struct sockaddr_in serv;

    int sock = socket ( AF_INET, SOCK_DGRAM, 0);

    //Socket could not be created
    if(sock < 0)
    {
        perror("Socket error");
    }

    memset( &serv, 0, sizeof(serv) );
    serv.sin_family = AF_INET;
    serv.sin_addr.s_addr = inet_addr( google_dns_server );
    serv.sin_port = htons( dns_port );
 
    int err = connect( sock , (const struct sockaddr*) &serv , sizeof(serv) );

    struct sockaddr_in name;
    socklen_t namelen = sizeof(name);
    err = getsockname(sock, (struct sockaddr*) &name, &namelen);

    char buffer[128];
#ifdef _WIN32
    const char* p = inet_ntoa(name.sin_addr);
#else
    const char* p = inet_ntop(AF_INET, &name.sin_addr, buffer, 128);
#endif

    if(p != NULL)
    {
#ifdef _WIN32
        printf ("main: Node [%d] Proc.Name [%s]:  local ip=[%s] \n", host_rank, processor_name, buffer);
#else
        printf ("main: Node [%d] Proc.Name [%s]:  local ip=[%s] \n", host_rank, processor_name, buffer);
#endif
        //printf("Local IP is : %s \n" , buffer);
    }
    else
    {
        //Some error
        printf ("main: Error number : %d . Error message : %s \n" , err , strerror(err));
    }
#ifdef _WIN32
	closesocket(sock);
#else
    close(sock);
#endif
	
}


typedef int (*fncFP)( ... );

extern "C" int zvode_(fncFP, long int *, dblcmplx *,
	    double *, double *, long int *, double *, double *,
	     long int *, long int *, long int *, dblcmplx *, long int *,
	    double *, long int *, long int *, long int *, fncFP, long int *,
	    dblcmplx *, long int *);


int myfex_fnc_f1(	long int *NEQ,
			double *T,
			dblcmplx *Y,
			dblcmplx *YDOT,
			dblcmplx *RPAR,
			long int *IPAR)
{
    simpleComplex<double> o0, o1, out0, out1;

      o0.re=0.0;   o0.im=0.0; o1.re=0.0;   o1.im=0.0;
    out0.re=0.0; out0.im=0; out1.re=0.0; out1.im=0.0;

    o0 = Y[0] * H[0];
    o1 = Y[1] * H[1];

    out0 = o0 + o1;

    o0.re=0.0;   o0.im=0.0; o1.re=0.0;   o1.im=0.0;

    o0 = Y[0] * H[2];
    o1 = Y[1] * H[3];

    out1 = o0 + o1;

    YDOT[0] = out0;
    YDOT[1] = out1;

	return 0;
}


#if 0
int myfex(	long int *NEQ,
			double *T,
			dblcmplx *Y,
			dblcmplx *YDOT,
			dblcmplx *RPAR,
			long int *IPAR)
{
	YDOT[0] = mone * (*RPAR) * Y[0] * Y[0] * Y[1];
	YDOT[1] = (*RPAR) * Y[1];

	return 0;
}
#endif


int dummyjex(	long int *NEQ,
			double *T,
			dblcmplx *Y,
			long int *ML,
			long int *MU,
			dblcmplx *PD,
			long int *NRPD,
			dblcmplx *RPAR,
			long int *IPAR)
{
	return 0;
}




// cmd_code[0] -- code of command
// cmd_code [1] -- 0xAFFFFF all nodes
//              -- less than 0xAFFFFF number of nude
int cmd_code[2];

char processor_name[MPI_MAX_PROCESSOR_NAME];

void _small_delay(unsigned int mseconds)
{
    clock_t goal = mseconds + clock();
 
	while (goal > clock()) { };
}

int rng_test_1(int host_rank)
{
    int i;

	unsigned int seed0, seed1, seed2, seed3;
	
	seed0 = time(0) + host_rank; _small_delay( 567 );
	seed1 = time(0) + host_rank; _small_delay( 127 );
	seed2 = time(0) + host_rank; _small_delay( 723 );
	seed3 = time(0) + host_rank; _small_delay( 901 );
	
	lfsr113_generator_init(seed0, seed1, seed2, seed3);
	
	

    for(i = 0 ; i < 10 ; i++)
    {
        printf ("rng_test: Node [%d] Proc.Name [%s]: %d %f \n", host_rank, processor_name, i, lfsr113_genRand_asDouble());
    }

    return 0;
}


template <typename TYPE, size_t SIZE, size_t AlphaSize>
int zvode_method_for_mc(TYPE h, TYPE &_Tout_par, int steps,
               simpleComplex<TYPE> &_T_par,
               uVector< simpleComplex<TYPE>, SIZE > &_Y_par,
               int (*fnc)(long int *NEQ, TYPE *T, dblcmplx *Y, dblcmplx *YDOT, dblcmplx *RPAR, long int *IPAR) )
{
	int i = 0;

	int errLvl = 0;

		zvode_((fncFP)fnc,
			&neq, _Y_par.m, &_T_par.re, &_Tout_par, &itol, &rtol, &atol_val, &itask,
			&istate, &iopt, zwork, &lzw, rwork, &lrw, iwork, &liw,
			(fncFP)dummyjex,
			&mf, &rpar, &ipar);


		if(istate < 0)
		{
			cout << "Error halt.  ISTATE = " << istate << endl ;
			errLvl = -1;
			//break;
		}


	return errLvl;
}

template <typename TYPE, size_t SIZE>
uVector< simpleComplex<TYPE>, SIZE > fnc_t1(const simpleComplex<TYPE> T, const uVector< simpleComplex<TYPE>, SIZE > &Y)
{
    simpleComplex<TYPE> o0, o1, out0, out1;
    uVector< simpleComplex<TYPE>, SIZE > YDOT;

      o0.re=0.0;   o0.im=0.0; o1.re=0.0;   o1.im=0.0;
    out0.re=0.0; out0.im=0; out1.re=0.0; out1.im=0.0;

    o0 = Y[0] * H[0];
    o1 = Y[1] * H[1];

    out0 = o0 + o1;

    o0.re=0.0;   o0.im=0.0; o1.re=0.0;   o1.im=0.0;

    o0 = Y[0] * H[2];
    o1 = Y[1] * H[3];

    out1 = o0 + o1;

    YDOT[0] = out0;
    YDOT[1] = out1;

    return YDOT;

}


void send_to_all_nodes(int c1, int c2, int num_procs)
{
	int i=0;

	cmd_code[0] = c1;
	cmd_code[1] = c2;

	for( i = 1 ; i < num_procs ; ++i )
	{
		MPI_Send(&cmd_code[0], 2, MPI_INT, i, CMD_DATA_TAG, MPI_COMM_WORLD);
	}
}

template<size_t N, size_t Ntrj, size_t _WV_LEAD_DIM, size_t _WV_LEAD_DIM_SQR, size_t _C_OPS_SIZE>
int process_trajectories(double _from_time, double _to_time,
						 int host_rank,
						 int use_colappse_operator, int use_expecation_operator)
{
	int i = 0, j=0, k = 0, k_ons = 0, trj = 0, odesolverstate = 0, ode_norm_steps = 5, cnt = 0;
	int m = 5;

	double a =  _from_time;
	double b = _to_time;
	double h = (b - a) / ((double)N - 1);
	//double aa,bb;
	double hh;

	double mu = 0.0, nu = 0.0,
        norm2_prev = 0.0,
        norm2_psi = 0.0,
        norm2_guess= 0.0,
        ode_norm_tol = 0.0001, sump = 0.0 ;

		
    uVector< simpleComplex<double>, N > tlist;

    struct uMatrix<simpleComplex<double>, Ntrj> pt_trjs;
    pt_trjs.rows = Ntrj ;
    pt_trjs.cols = N ;
    uVector< simpleComplex<double>, N > pt_avg_trj;

	dblcmplx mone;
	dblcmplx mtwo;

	
	mone.re = -1;
	mone.im =  0;

	mtwo.re = -2;
	mtwo.im =  0;


	unsigned int seed0, seed1, seed2, seed3;
	
	seed0 = time(0) + host_rank; _small_delay( 567 );
	seed1 = time(0) + host_rank; _small_delay( 127 );
	seed2 = time(0) + host_rank; _small_delay( 723 );
	seed3 = time(0) + host_rank; _small_delay( 901 );
	
	lfsr113_generator_init(seed0, seed1, seed2, seed3);
	
	printf ("process_trajectories: Node [%d] Proc.Name [%s]: begin \n", host_rank, processor_name);


    uVector< simpleComplex<double>, N > P;

    simpleComplex<double> T ;
    simpleComplex<double> T_prev ;
    simpleComplex<double> T_final ;
    simpleComplex<double> T_guess ;


    uVector< simpleComplex<double>, _WV_LEAD_DIM > Y;
    uVector< simpleComplex<double>, _WV_LEAD_DIM > Y2;
    uVector< simpleComplex<double>, _WV_LEAD_DIM > Y_prev ;
    uVector< simpleComplex<double>, _WV_LEAD_DIM > Y_tmp;

    uVector< simpleComplex<double>, _WV_LEAD_DIM > out_psi;
	
	uVector< simpleComplex<double>, _WV_LEAD_DIM_SQR> id_operator;
	
	k=0;
	for(i=0;i<_WV_LEAD_DIM;i++)
	{
		for(j=0;j<_WV_LEAD_DIM;j++)
		{
			if( (i==j) )
			{
				id_operator[k] = make_simpleComplex( 1.0, 0.0 );
			}
			else
			{
				id_operator[k] = make_simpleComplex( 0.0, 0.0 );
			}
			k++;
		}
		k++;
	}
	k=0; i=0; j=0;

    struct simpleComplex<double>  ev;

	neq = _WV_LEAD_DIM;
	itol = 1;
	rtol = 1e-8;
	atol_val = 1e-7;
	itask = 5;
	istate = 1;
	iopt = 0;
	lzw = 30;
	lrw = 22;
	liw = 30;
	mf = 10;  // 10 for nonstiff (Adams) method, no Jacobian used.
	rpar.re = 0.0;
	rpar.im = 0.0;
	aemax = 0.0;

//          10 for nonstiff (Adams) method, no Jacobian used.
//          21 for stiff (BDF) method, user-supplied full Jacobian.
//          22 for stiff method, internally generated full Jacobian.
//          24 for stiff method, user-supplied banded Jacobian.
//          25 for stiff method, internally generated banded Jacobian.

		memset( iwork, 0, sizeof(iwork) );
		memset( rwork, 0, sizeof(rwork) );
		memset( zwork, 0, sizeof(zwork) );

		for( i = 0 ; i < N ; i++)
		{
			tlist[i].re = ( ( (i+1) - 1) * h ) ;
			tlist[i].im = 0.0 ;
		}


	printf ("process_trajectories: Node [%d] Proc.Name [%s]:  tlist done \n", host_rank, processor_name);

	for ( trj = 0 ; trj < Ntrj ; trj++)
    {
        mu = lfsr113_genRand_asDouble();
        nu = lfsr113_genRand_asDouble();

        for( i = 0 ; i < _WV_LEAD_DIM ; i++)
        {
            Y[i] = alpha[i];
        }

		istate = 1;

			if(use_expecation_operator == 1)
			{
				ev = expect_cnv_denmat< simpleComplex<double>, _WV_LEAD_DIM_SQR, _WV_LEAD_DIM>(_WV_LEAD_DIM, _WV_LEAD_DIM, expect_operator, Y);
				pt_trjs(trj, 0) = ev;
			} else {
				ev = expect_cnv_denmat< simpleComplex<double>, _WV_LEAD_DIM_SQR,_WV_LEAD_DIM>(_WV_LEAD_DIM, _WV_LEAD_DIM, id_operator, Y);
				pt_trjs(trj, 0) = ev;
			}


            int counterIter = 0;

			printf ("process_trajectories: Node [%d] Proc.Name [%s]:  trj %d \n", host_rank, processor_name, trj);

			T=tlist[0];

            for ( k = 0 ; k < N; k++)
            {
				if( k == 0)
					tout = T.re;
				else
					tout = tlist[k].re;

				rwork[0]=tout;

				norm2_psi = norm(Y);
                while( T.re < tout )
                {

                    T_prev = T;
                    Y_tmp = Y;

                    norm2_prev = norm2_psi;

					odesolverstate = zvode_method_for_mc<double,_WV_LEAD_DIM,_WV_LEAD_DIM>(h, tout, 1,  T, Y, &myfex_fnc_f1);

                    norm2_psi = norm( Y );

                    if(norm2_psi <= mu )
                    {
						printf ("process_trajectories: Node [%d] Proc.Name [%s]:  %s \n", host_rank, processor_name, "Collpase has occured.");

                        T_final = T ;
                        cnt = 0 ;

                        for( k_ons=0; k_ons < ode_norm_steps ; k_ons++)
                        {
                            T_guess = T_prev + log(norm2_prev/mu) / log(norm2_prev/norm2_psi) * (T_final-T_prev) ;

                            if ( (T_guess.re < T_prev.re) || (T_guess.re > T_final.re) )
                            {
                                T_guess.re = T_prev.re + 0.5 * (T_final.re - T_prev.re) ;
                            }

                            Y = Y_tmp ;
                            T = T_prev ;

							saveitask=itask;
							itask=1;
							odesolverstate = zvode_method_for_mc<double,_WV_LEAD_DIM,_WV_LEAD_DIM>(h, T_guess.re, 1, T, Y, &myfex_fnc_f1);
							itask=saveitask;
							printf ("process_trajectories: Node [%d] Proc.Name [%s]: istate=%d %s \n", host_rank, processor_name, istate, "zvode_method_for_mc collapse call 2");

                            norm2_guess = norm( Y ) ;
                            if ( abs(mu - norm2_guess) < (ode_norm_tol * mu) )
                            {
                                break ;
                            }
                            else if (norm2_guess < mu)
                            {
                                T_final = T_guess ;
                                norm2_psi = norm2_guess ;
                            }
                            else
                            {
                                T_prev = T_guess ;
                                Y_tmp = Y ;
                                norm2_prev = norm2_guess ;
                            }
                            cnt = cnt + 1 ;
                        } // for( k=1; k < ode_norm_steps ; k++)

                        if (cnt > ode_norm_steps)
                        {
							printf ("process_trajectories: Node [%d] Proc.Name [%s]:  %s \n", host_rank, processor_name, "Norm tolerance not reached. Increase accuracy of ODE solver or norm_steps.");
                            break;
                        }
                        else
                        {
							printf ("process_trajectories: Node [%d] Proc.Name [%s]:  %s \n", host_rank, processor_name, "Norm tolerance is reached.");
                        }

                        for ( j = 0 ; j < _C_OPS_SIZE; j++ )
                        {
                            Y_tmp = c_ops[j] * Y ;
                            P[j].re = norm ( Y_tmp );
                        }

                        P = P / sum(P).re;
                        sump = 0.0 ;

                        for ( j=0 ; j < _C_OPS_SIZE; j++ )
                        {
                            if ( (sump <= nu) && (nu < (sump + P[j].re) ) )
                            {
                                Y = c_ops[j] * Y ;  
                            }
                            sump = sump + P[j].re ;
                        } // for ( j=0 ; j < _C_OPS_SIZE; j++ )

                        // new random numbers
                        mu = lfsr113_genRand_asDouble();
                        nu = lfsr113_genRand_asDouble();

                        normalize(Y);

                     // reset, first call to zvode
                     istate = 1 ;

                    } // if(norm2_psi <= mu )


                } // while(T.re < tlist[k].re)

				Y_tmp = Y;
				normalize(Y_tmp);

                out_psi = Y_tmp / normsqrt(Y_tmp);

				if(use_expecation_operator == 1)
					ev = expect_cnv_denmat< simpleComplex<double>, _WV_LEAD_DIM_SQR, _WV_LEAD_DIM>(_WV_LEAD_DIM, _WV_LEAD_DIM, expect_operator, out_psi);
				else
					ev = expect_cnv_denmat< simpleComplex<double>, _WV_LEAD_DIM_SQR, _WV_LEAD_DIM>(_WV_LEAD_DIM, _WV_LEAD_DIM, id_operator, out_psi);

                pt_trjs(trj, k) = ev;

            } // for ( k = 1 ; k < N; k++)
    } // for ( trj = 0 ; trj < Ntrj ; trj++)


    /**
     *
     * Average of trajectories
     *
     */

    simpleComplex<double> tmp;

    for ( k = 0 ; k < N; k++)
    {
        tmp.re=0; tmp.im=0;
        for ( trj = 0 ; trj < Ntrj ; trj++)
        {
            tmp = tmp + pt_trjs(trj, k);
        }
        tmp = tmp / (double)Ntrj;
        pt_avg_trj[k] = tmp;
    }


	printf ("process_trajectories: Node [%d] Proc.Name [%s]:  send TRJ to master node \n", host_rank, processor_name);

	MPI_Send(&pt_avg_trj.m[0], 2*N, MPI_DOUBLE, MASTER_NODE, TRJ_DATA_TAG, MPI_COMM_WORLD);

	printf ("process_trajectories: Node [%d] Proc.Name [%s]:  end \n", host_rank, processor_name);

	return 0;
}


template<size_t N, size_t Ntrj, size_t _WV_LEAD_DIM, size_t _WV_LEAD_DIM_SQR, size_t _C_OPS_SIZE>
int process_codes(int c1, int c2, int host_rank, 
						 double _from_time, 
						 double _to_time,
						 int use_colappse_operator, int use_expecation_operator)
{
	int ret_code = 0;

	if(cmd_code[0] == cmd_RNG_TEST)
	{
		printf("Work Node Node [%d] Proc.Name [%s] -- cmd_RNG_TEST\n", host_rank, processor_name);
		rng_test_1( host_rank );
	} //  if(cmd_code[0] == cmd_RNG_TEST)


	if(cmd_code[0] == cmd_TRJ_PROC)
	{
		printf("Work Node Node [%d] Proc.Name [%s] -- cmd_RNG_TEST\n", host_rank, processor_name);
		process_trajectories<N, Ntrj, _WV_LEAD_DIM, _WV_LEAD_DIM_SQR, _C_OPS_SIZE>( _from_time, _to_time, 
			host_rank, 
			use_colappse_operator, use_expecation_operator );
	}

	if(cmd_code[0] == cmd_EXIT)
	{
		printf("Work Node Node [%d] Proc.Name [%s] -- cmd_EXIT\n", host_rank, processor_name);
		ret_code = cmd_EXIT;
	} // if(cmd_code[0] == cmd_EXIT)

	return ret_code;
} // process_codes


template<size_t N, size_t Ntrj, size_t _WV_LEAD_DIM, size_t _WV_LEAD_DIM_SQR, size_t _C_OPS_SIZE>
int mpi_main(int argc, char *argv[], int verbose_mode,
						 double _from_time, 
						 double _to_time,
						 int use_colappse_operator, int use_expecation_operator)
{
	int i, k, lvl = 0;
	int max_buf_len = 0, errlvl = 0;
	int host_rank, num_procs;
	
	struct uMatrix<simpleComplex<double>, Ntrj> glb_trjs;
	uVector< simpleComplex<double>, N > final_avg_trj;


	MPI_Init( &argc, &argv );

	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &host_rank);

	MPI_Get_processor_name(processor_name, &max_buf_len);

	if (num_procs < 2)
	{
		fprintf(stderr, "Process numbers must be greater than one for %s\n", argv[0]);
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	print_ip(host_rank, processor_name);

	if(host_rank == MASTER_NODE)
	{
		printf("Master Node Node [%d] Proc.Name [%s]\n", host_rank, processor_name);
		printf("Master Node Node [%d] Proc.Name [%s]: N=%d \n", host_rank, processor_name, N);

		send_to_all_nodes(cmd_TRJ_PROC, 0xAFFFFF, num_procs);

		// receive trajectories from computational nodes

		for( k=1; k < num_procs; k++ )
		{
			uVector< simpleComplex<double>, N > avg_trj;

			avg_trj.size=N;
			for( i = 0 ; i < N ; i++)
			{
				avg_trj[i].re = 0.0 ;
				avg_trj[i].im = 0.0 ;
			}

			MPI_Recv(&avg_trj.m[0], 2*N, MPI_DOUBLE, MPI_ANY_SOURCE, TRJ_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			for( i = 0 ; i < N ; i++)
			{
			 	 glb_trjs(k, i) = avg_trj[i];
			}


			printf ("main: Node [%d] Proc.Name [%s]:  recv TRJ from work node \n", host_rank, processor_name);
			for( i = 0 ; i < N ; i++)
			{
				printf("data from node %d in trj[%d] = %f \n", k, i, glb_trjs(k, i).re);
			}
			printf ("main: Node [%d] Proc.Name [%s]:  end data \n", host_rank, processor_name);
		}

		send_to_all_nodes(cmd_EXIT, GLOBAL_CODE, num_procs);

		simpleComplex<double> tmp;

		for ( i = 0 ; i < N ; i++)
		{
			tmp.re=0; tmp.im=0;
			for( k=1; k < num_procs; k++ )
			{
				tmp = tmp + glb_trjs(k, i).re;
			}
			tmp = tmp / (double)(num_procs - 1);
			glb_trjs(0, i) = tmp;
		}

		printf ("main: Node [%d] Proc.Name [%s]:  avg TRJ begin \n", host_rank, processor_name);
		for( i = 0 ; i < N ; i++)
		{
				printf("main: data from node in avg_trj[%d] = %f \n", i, glb_trjs(0, i).re);
		}
		printf ("main: Node [%d] Proc.Name [%s]:  avg TRJ end data \n", host_rank, processor_name);

		printf ("main: Node [%d] Proc.Name [%s]:  dump avg TRJ begin data \n", host_rank, processor_name);
		for( i = 0 ; i < N - 1 ; i++)
		{
			printf("%f, ", glb_trjs(0, i).re);
			if( (i % 10 ) == 0 ) printf("\n");
		}
		printf("%f \n", glb_trjs(0, i).re);

		printf ("main: Node [%d] Proc.Name [%s]:  dump avg TRJ end data \n", host_rank, processor_name);

	} // if(host_rank == MASTER_NODE)
	else
	{
		printf("Work Node Node [%d] Proc.Name [%s]\n", host_rank, processor_name);
		printf("Work Node Node [%d] Proc.Name [%s]: N=%d \n", host_rank, processor_name, N);

		while(1)
		{
			cmd_code[0] = 0;
			cmd_code[1] = 0;

			MPI_Recv(&cmd_code, 2, MPI_INT, MPI_ANY_SOURCE, CMD_DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			if(cmd_code[1] == GLOBAL_CODE) // global code
			{
				lvl = 0;
				lvl = process_codes<N,Ntrj,_WV_LEAD_DIM,_WV_LEAD_DIM_SQR,_C_OPS_SIZE>(cmd_code[0], cmd_code[1], 
						host_rank,
						_from_time, _to_time,
						use_colappse_operator, use_expecation_operator);
				if( lvl == cmd_EXIT)
					break;
			}
			else // single code for single node
			{
				if(host_rank == cmd_code[1])
				{
					lvl = 0;
					lvl = process_codes<N,Ntrj,_WV_LEAD_DIM,_WV_LEAD_DIM_SQR,_C_OPS_SIZE>(cmd_code[0], cmd_code[1], 
						host_rank,
						_from_time, _to_time,
						use_colappse_operator, use_expecation_operator);
					if( lvl == cmd_EXIT)
						break;
				}
			}
		}
	} // if(host_rank == MASTER_NODE) else

	MPI_Finalize( );

	return errlvl;
} // int mpi_main(int argc, char *argv[])
