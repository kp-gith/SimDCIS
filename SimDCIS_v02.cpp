// SimDCIS: Markov state machine for DCIS
// Copyright (C) 2024 Marcel Greuter / Keris Poelhekken

//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.

//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include <iostream> 
#include <cstdlib> 
#include <time.h>
#include <cmath>
#include <conio.h>

using namespace std; 

int state = 0 ;
int grade = 0 ;
int age = 0 ;
float p1[101] ;
float p2[101][3], p3[101][3], p4[101][3] ;
int Nmam = 0 ;
int Ndeath = 0 ;
int Nreg = 0 ;
int NSD[20][3] ;
int NCD[20][3] ;
int NIBC[20][3] ;
int TotSD ;
int TotCD ; 
int Ng[3] ;
int Nmami[13] ;
float TotDeathAge, TotDetected, TotClindet ;

double rnd1[100] ;
double rnd2[100] ;

float Compliance, Sensitivity, ClinicalDet ;
int Niterations, Npopulation ;

void ReadTransitionProbabilities() 
{
	FILE *fp;
	int i, j ;
	int dum ;
	
	if ((fp = fopen( "transition.txt", "r"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open input file transition.txt." << std::endl;
		exit(1);
	}
	printf( "Reading transitions parameters from transition.txt\n" ) ;
	for ( i = 0 ; i < 101 ; i++ ) 
	{							
		fscanf( fp, "%d\t%f\t", &dum, &p1[i] ) ;
		for ( j = 0 ; j < 3 ; j++ )
			fscanf( fp, "%f\t", &p2[i][j] ) ; 
		for ( j = 0 ; j < 3 ; j++ )
			fscanf( fp, "%f\t", &p3[i][j] ) ;
		for ( j = 0 ; j < 3 ; j++ )
			fscanf( fp, "%f\t", &p4[i][j] ) ;  	
	}
	fclose( fp ) ;
}

void ReadInputParams( void ) 
{
	FILE *fp;
	char dum[20] ;
	if ((fp = fopen( "inputparams.txt", "r"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open input file inputparams.txt." << std::endl;
		exit(1);
	}
	printf( "Reading input parameters from inputparams.txt\n" ) ;
	fscanf( fp, "%s %f\n", &dum, &Compliance ) ;
	fscanf( fp, "%s %f\n", &dum, &Sensitivity ) ;
	fscanf( fp, "%s %f\n", &dum, &ClinicalDet ) ;
	fscanf( fp, "%s %d\n", &dum, &Npopulation ) ;
	fscanf( fp, "%s %d\n", &dum, &Niterations ) ;
	fclose( fp ) ;
}

double Random( void )
{
	return (double)rand() / (double)((unsigned)RAND_MAX + 1 ) ;
}

void Healthy( int age )												// state = 0
{
	if ( rnd1[age] < p1[age] )
	{
		state = 1 ;
		Ndeath++ ;
	}
	else
	{
		if ( rnd1[age] < p1[age]+p2[age][0] )
		{
			state = 2 ;
			grade = 1 ;
			Ng[0]++ ;
		}
		else
		{
			if ( rnd1[age] < p1[age]+p2[age][0]+p2[age][1] )
			{	
				state = 2 ;
				grade = 2 ;
				Ng[1]++ ;
			}	
			else
			{
				if ( rnd1[age] < p1[age]+p2[age][0]+p2[age][1]+p2[age][2] )
				{	
					state = 2 ;
					grade = 3 ;
					Ng[2]++ ;
				}		
			}
		}
	} 
}

void DCIS( int age, int grade )										// state = 2
{
	int j = grade - 1 ;												
	if ( rnd1[age] < p1[age] )
	{
		state = 1 ;
		Ndeath++ ;
	}
	else
	{
		if ( rnd1[age] < p1[age]+p3[age][j] )
		{
			state = 0 ;
			grade = 0 ;
			Nreg++ ;
		}
		else
		{
			if ( rnd1[age] < p1[age]+p3[age][j]+p4[age][j] )
			{
				state = 3 ;
				NIBC[age/5][grade-1]++ ;
			}
		}
	}
}

int DoScreening( int age )
{
	double rndc = Random() ;
	int result = 0 ;
	int j ;
	if ( (age==50)||(age==52)||(age==54)||(age==56)||(age==58)||(age==60)||(age==62)||(age==64)||(age==66)||(age==68)||(age==70)||(age==72)||(age==74))
	{
		j = (age-50)/2 ;
		if ( j < 13 )
		{
			if ( rndc < Compliance )									// Compliance
			{
				Nmam++ ;
				Nmami[j]++ ;
				if ( ( rnd2[age] < Sensitivity ) && ( state == 2 ) )	// sensitivity for detecting DCIS
				{
					result = 4 ;
					NSD[age/5][grade-1]++;
					TotSD++ ;
				}
			}	
		}
	}
	return result ;
}

int ClinicalDetection( int age ) 
{
	int result = 0 ;
	if (  ( rnd2[age] < ClinicalDet ) && ( state == 2 ) )
	{
		result = 5 ;
		NCD[age/5][grade-1]++;
		TotCD++ ;
	}
	return result ;
}

void MakeRandomNumbers( int pop, int it )
{
	int i ;
	srand( pop+1+it*Npopulation ) ;
	for ( i = 0 ; i < 100 ; i++ )
	   rnd1[i] = Random() ;
	srand( 2*(pop+1)+it*Npopulation ) ;
	for ( i = 0 ; i < 100 ; i++ )
	   rnd2[i] = Random() ;
}

void InitParams( void )
{
	int i, j ;
	
	Nmam = 0 ;
	Ndeath = 0 ;
	Nreg = 0 ;
	
	for ( i = 0 ; i < 3 ; i++ ) 
		Ng[i] = 0 ;

	for ( i = 0 ; i < 20 ; i++ )
	for ( j = 0 ; j < 3 ; j++ )
	{
		NSD[i][j] = 0 ;
		NCD[i][j] = 0 ;
		NIBC[i][j] = 0 ;
	}
	TotSD = 0 ;
	TotCD = 0 ;
	
	TotDeathAge = 0 ;
	TotDetected = 0 ;
	TotClindet = 0 ;
	for ( i = 0 ; i < 13 ; i++ )
		Nmami[i] = 0 ;
}

int main()
{
	int pop ;
	FILE *fp1, *fp2 , *fp3;
	int death, IBC, dum, detected, clindet ;
	int it, i, j ;
	char buf[2], fn[80] ;

	printf( "DCIS Simulation SimDCIS\n" ) ;
	printf( "Version 14d - 20 apr 2024\n" ) ;

	ReadTransitionProbabilities() ;
	ReadInputParams() ;
	printf( "Compliance: %.2f\n", Compliance ) ;
	printf( "Sensitivity: %.2f\n", Sensitivity ) ;
	printf( "Clinical detection: %.2f\n", ClinicalDet ) ;
	printf( "Number of participants: %d\n", Npopulation ) ;
	printf( "Number of iterations: %d\n", Niterations ) ;
	printf( "simdcisx.out: for iteration x for every participant state per age\n" ) ;
	printf( "simdcis1.det: statistics per iteration\n" ) ;
	printf( "simdcis2.det: grades per age per iteration\n" ) ;
	if ((fp2 = fopen( "simdcis1.det", "w"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open output file simdcis1.det" << std::endl;
		exit(1);
	}
	if ((fp3 = fopen( "simdcis2.det", "w"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open output file simdcis2.det" << std::endl;
		exit(1);
	}
	fprintf( fp2, "Iter\tNmam\tNdeath\tNreg\t" ) ;
	for ( i = 0 ; i < 13; i++ )
		fprintf( fp2, "Nmam%d\t", 50+i*2 ) ;
	fprintf( fp2, "AvDeath\tAvSDAge\tAvCDAge\n")	;
	fprintf( fp3, "\tAge\tScreen\t\t\t\tClinical\t\t\tProgressed\n" ) ;
	fprintf( fp3, "Iter\tGrade\tAll\t1\t2\t3\tAll\t1\t2\t3\tAll\t1\t2\t3\n" ) ;
	for ( it = 0 ; it < Niterations ; it++ )
	{
		itoa( it, buf, 10 ) ;
  		strcpy( fn, "simdcis" ) ;
  		strcat( fn, buf ) ;
  		strcat( fn, ".out" ) ;
		if ((fp1 = fopen( fn, "w"))==NULL) {
			std::cerr << std::endl << fn << std::endl;
			exit(1);
		}
		InitParams() ;
		for ( pop = 0 ; pop < Npopulation ; pop++ )
		{
			if ((pop % 1000 == 0 ))
				printf( "." ) ;
			fprintf( fp1, "%d\t", pop ) ;
			MakeRandomNumbers( pop, it ) ;
			age = 0 ;
			state = 0 ;
			grade = 0 ;
			death = 0 ;
			IBC = 0 ;
			detected = 0 ;
			clindet = 0 ;
			while ( ( age < 101 ) && (!death) && (!IBC) && (state!=1) && (state!=3) && (state!=4) && (state!=5) ) 
			{
				switch ( state )
				{
					case 0: Healthy( age ) ; break ;
					case 1: death = 1 ; break ;
					case 2: DCIS( age, grade ) ; break ;
					case 3: IBC = 1 ; break ;
					case 4: detected = 1 ; break ;
					case 5: clindet = 1 ; break ; 
				} 
				if ( (state!=1) && (state!=3) && (state!=4) && (state!=5) )
				{
					dum = DoScreening( age ) ;
					if ( dum != 0 )
						state = dum ; 
					else
					{
						dum = ClinicalDetection( age ) ;
						if ( dum != 0 )
							state = dum ;
					}
				}
				if ( state == 1 )
				  TotDeathAge += age ;
				if ( state == 4 ) 
				  TotDetected += age ;
				if ( state == 5 ) 
				  TotClindet += age ; 
				fprintf( fp1, "%1d%1d\t", state, grade ) ;
				age++ ;
			}
			fprintf( fp1, "\n" ) ;
		}
		fclose( fp1 ) ;
		fprintf( fp2, "%d\t%d\t%d\t%d\t", it, Nmam, Ndeath, Nreg ) ;
		for ( i = 0 ; i < 13 ; i++ )
			fprintf( fp2, "%d\t", Nmami[i] ) ;
		fprintf( fp2, "%.2f\t%.2f\t%.2f\n", TotDeathAge / (float)Ndeath, TotDetected / (float)TotSD, TotClindet / (float)TotCD ) ;
		for ( i = 0 ; i < 20 ; i++ )
		{
			fprintf( fp3, "%d\t%d\t", it, i*5 ) ;
			fprintf( fp3, "%d\t", NSD[i][0] + NSD[i][1] + NSD[i][2] ) ;
			for ( j = 0 ; j < 3 ; j++ )
				fprintf( fp3, "%d\t", NSD[i][j] ) ;
			fprintf( fp3, "%d\t", NCD[i][0] + NCD[i][1] + NCD[i][2] ) ;
			for ( j = 0 ; j < 3 ; j++ )
				fprintf( fp3, "%d\t", NCD[i][j] ) ;	
			fprintf( fp3, "%d\t", NIBC[i][0] + NIBC[i][1] + NIBC[i][2] ) ;
			for ( j = 0 ; j < 3 ; j++ )
				fprintf( fp3, "%d\t", NIBC[i][j] ) ;	
			fprintf( fp3, "\n" ) ;
		}
	}
	fclose( fp2) ;
	fclose( fp3) ;
}