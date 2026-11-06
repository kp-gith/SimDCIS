// Markov state machine for DCIS
// Marcel Greuter / Keris Poelhekken
//
// Version 14e: time step of 0.5 year (instead of 1 year).
// - transition.txt contains 201 lines: age 0.0, 0.5, ..., 100.0 with
//   transition probabilities per 0.5 year
// - all arrays that were kept per age are now kept per half-year step (index t)
// - age in years = t * DT
// Version 14f: added
// - specificity
// - costs of screening, biopsy and treatment
// - discounting per ref_age
// - In the det1 file the total costs of screening, biopsy and treatment and total FP per iteration
// Version 15.01
// - FP calculation corrected
// Version 15.02
// - screening ages are read from inputparams.txt (line "ScreenAges")
// - discount exponent corrected: t*DT instead of t/2 (integer division rounded down for odd t)
// Version 15.03
// - 2nd screening modality added with its own sensitivity, specificity, costs
//   and screening ages (inputparams.txt: Sensitivity2, Specificity2, cost_screen2,
//   cost_biopsy2, cost_treatment2, ScreenAges2). These lines are optional: if they are
//   absent, the program runs as before with modality 1 only.
// - output of modality 2 in separate files: simdcis1_mod2.det and simdcis2_mod2.det
// - in the simdcisx.out files: state 4 = detected at screening with modality 1,
//   state 6 = detected at screening with modality 2
// - seeds corrected: every (iteration, participant, stream) gets a unique seed. In 15.02
//   the seeds overlapped (pop+1, pop+2, pop+3), so participant p+1 used the same
//   numbers for rnd1 as participant p for rnd2. Results therefore differ slightly
//   (statistically negligible) from 15.02 and are not bit-for-bit reproducible.
// - rand() (15 bit under Visual C++) replaced by std::mt19937; the compliance of
//   modality 1 now also has its own fixed stream (rnd0), so results for modality 1
//   do not depend on whether modality 2 is active.
// - costs are now summed in double (instead of float)


#include <iostream> 
#include <cstdlib> 
#include <cstdio>
#include <cstring>
#include <time.h>
#include <cmath>
#include <conio.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <random>

using namespace std; 

const int   NSTEP   = 201 ;		// number of time steps: 0.0 to 100.0 years in steps of 0.5
const float DT      = 0.5f ;	// time step in years
const int   NAGEGRP = 20 ;		// age groups of 5 years (for output)
const int   STEPS_PER_GRP = 10 ;// 5 years = 10 time steps
const int   NMOD    = 2 ;		// number of screening modalities
const int   NSTREAM = 7 ;		// number of random streams per participant

// ---------------------------------------------------------------------------
// SCREENING AGES (in years) are read from inputparams.txt, lines:
//   ScreenAges   50 52 54 56 58 60 62 64 66 68 70 72 74     (modality 1)
//   ScreenAges2  ...                                          (modality 2, optional)
// Half years are allowed (e.g. 52.5). Ages must be a multiple of
// 0.5 (DT) and sorted in ascending order. The number of rounds follows from the list.
// If both modalities fall on the same age, both are performed
// (modality 1 first; if modality 1 detects, modality 2 is skipped).
// ---------------------------------------------------------------------------
const char *INPUTPARAMS_FILE = "inputparams.txt" ;				// default path of the parameter file

// All data kept per screening modality
struct Modality
{
	bool   active ;												// modality in use?
	float  sens, spec ;											// sensitivity, specificity
	float  cost_screen, cost_biopsy, cost_treatment ;
	vector<float> age ;											// screening ages in years
	vector<int>   step ;										// corresponding time steps (t = age / DT)
	vector<int>   Nmami ;										// number of participants per screening round
	int    Nmam, NFP, TotSD ;
	int    NSD[NAGEGRP][3] ;									// screen-detected per age group and grade
	double TotCostScreen, TotCostBiopsy, TotCostTreatment ;
	double TotDetected ;										// sum of ages at screen detection
} ;
Modality M[NMOD] ;

int state = 0 ;
int grade = 0 ;
int t = 0 ;						// time step (age = t * DT)
float p1[NSTEP] ;
float p2[NSTEP][3], p3[NSTEP][3], p4[NSTEP][3] ;
int Ndeath = 0 ;
int Nreg = 0 ;

int NCD[NAGEGRP][3] ;
int NIBC[NAGEGRP][3] ;

int TotCD ; 

int Ng[3] ;

float TotDeathAge, TotClindet ;

double rnd1[NSTEP] ; // transition probability
double rnd2[NSTEP] ; // detection probability determined by sensitivity (modality 1) and clinical detection
double rnd3[NSTEP] ; // false-positive probability determined by specificity (modality 1)
double rnd0[NSTEP] ; // compliance modality 1
double rnd4[NSTEP] ; // compliance modality 2
double rnd5[NSTEP] ; // sensitivity modality 2
double rnd6[NSTEP] ; // specificity modality 2

float Compliance, ClinicalDet ;								// ClinicalDet: probability per 0.5 year
int Niterations, Npopulation ;

float cost_discount, cost_refage ;

// age group (5 years) at time step t; age 100 is added to the last group
int AgeGroup( int t )
{
	int g = t / STEPS_PER_GRP ;
	if ( g >= NAGEGRP ) g = NAGEGRP - 1 ;
	return g ;
}

// discount factor at time step t
double DiscountFactor( int t )
{
	double exponent = std::max(static_cast<double>(0.0), t*DT - static_cast<double>(cost_refage) );
	return std::pow(1.0 + cost_discount, exponent);
}

void ReadTransitionProbabilities() 
{
	FILE *fp;
	int i, j ;
	float dum ;														// age in years (0.0, 0.5, ...)
	
	if ((fp = fopen( "transition.txt", "r"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open input file transition.txt." << std::endl;
		exit(1);
	}
	printf( "Reading transitions parameters (per 0.5 year) from transition.txt\n" ) ;
	for ( i = 0 ; i < NSTEP ; i++ ) 
	{							
		if ( fscanf( fp, "%f\t%f\t", &dum, &p1[i] ) != 2 )
		{
			std::cerr << std::endl << "Error: transition.txt ended at line " << i+1 << ", expected " << NSTEP << " lines." << std::endl;
			exit(1);
		}
		for ( j = 0 ; j < 3 ; j++ )
			fscanf( fp, "%f\t", &p2[i][j] ) ; 
		for ( j = 0 ; j < 3 ; j++ )
			fscanf( fp, "%f\t", &p3[i][j] ) ;
		for ( j = 0 ; j < 3 ; j++ )
			fscanf( fp, "%f\t", &p4[i][j] ) ;  	
	}
	fclose( fp ) ;
}

// convert the numbers in 'line' to a list of ages
void ParseAges( const char *line, vector<float> &ages )
{
	ages.clear() ;
	char *pos = (char *)line, *end ;
	while ( true )
	{
		double age = strtod( pos, &end ) ;
		if ( end == pos ) break ;									// no more numbers
		ages.push_back( (float)age ) ;
		pos = end ;
	}
}

void ReadInputParams( void ) 
{
	FILE *fp;
	char dum[20], line[1024] ;
	if ((fp = fopen( INPUTPARAMS_FILE, "r"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open input file " << INPUTPARAMS_FILE << "." << std::endl;
		exit(1);
	}
	printf( "Reading input parameters from %s\n", INPUTPARAMS_FILE ) ;
	fscanf( fp, "%s %f\n", dum, &Compliance ) ;
	fscanf( fp, "%s %f\n", dum, &M[0].sens ) ;
	fscanf( fp, "%s %f\n", dum, &M[0].spec ) ;
	fscanf( fp, "%s %f\n", dum, &ClinicalDet ) ;
	fscanf( fp, "%s %f\n", dum, &M[0].cost_screen ) ;
	fscanf( fp, "%s %f\n", dum, &M[0].cost_biopsy ) ;
	fscanf( fp, "%s %f\n", dum, &M[0].cost_treatment ) ;
	fscanf( fp, "%s %f\n", dum, &cost_discount ) ;
	fscanf( fp, "%s %f\n", dum, &cost_refage ) ;
	fscanf( fp, "%s %d\n", dum, &Npopulation ) ;
	fscanf( fp, "%s %d\n", dum, &Niterations ) ;

	// screening ages modality 1: keyword ScreenAges followed by the ages on the same line
	if ( fscanf( fp, "%19s", dum ) != 1 || strcmp( dum, "ScreenAges" ) != 0 )
	{
		std::cerr << std::endl << "Error: expected line 'ScreenAges <age1> <age2> ...' after Niterations in " << INPUTPARAMS_FILE << "." << std::endl;
		exit(1);
	}
	if ( fgets( line, sizeof(line), fp ) == NULL )
	{
		std::cerr << std::endl << "Error: no screening ages found in " << INPUTPARAMS_FILE << "." << std::endl;
		exit(1);
	}
	ParseAges( line, M[0].age ) ;
	if ( M[0].age.empty() )
	{
		std::cerr << std::endl << "Error: no screening ages found in " << INPUTPARAMS_FILE << "." << std::endl;
		exit(1);
	}
	M[0].active = true ;

	// optional: modality 2
	M[1].active = false ;
	if ( fscanf( fp, "%19s %f", dum, &M[1].sens ) == 2 )
	{
		if ( fscanf( fp, "%19s %f", dum, &M[1].spec ) != 2 ||
		     fscanf( fp, "%19s %f", dum, &M[1].cost_screen ) != 2 ||
		     fscanf( fp, "%19s %f", dum, &M[1].cost_biopsy ) != 2 ||
		     fscanf( fp, "%19s %f", dum, &M[1].cost_treatment ) != 2 ||
		     fscanf( fp, "%19s", dum ) != 1 || strcmp( dum, "ScreenAges2" ) != 0 )
		{
			std::cerr << std::endl << "Error: incomplete modality 2 block in " << INPUTPARAMS_FILE << " (expected Sensitivity2, Specificity2, cost_screen2, cost_biopsy2, cost_treatment2, ScreenAges2)." << std::endl;
			exit(1);
		}
		if ( fgets( line, sizeof(line), fp ) != NULL )
			ParseAges( line, M[1].age ) ;
		M[1].active = !M[1].age.empty() ;
	}
	fclose( fp ) ;
}

// check screening ages and convert them to time steps
void SetupScreenSteps( Modality &m, int nr )
{
	int i, n = (int)m.age.size() ;
	m.step.assign( n, 0 ) ;
	for ( i = 0 ; i < n ; i++ )
	{
		m.step[i] = (int)floor( m.age[i] / DT + 0.5 ) ;
		if ( fabs( m.step[i] * DT - m.age[i] ) > 1e-4 || m.step[i] < 0 || m.step[i] >= NSTEP )
		{
			std::cerr << std::endl << "Error: screening age " << m.age[i] << " (modality " << nr << ") is not a multiple of " << DT << " or out of range." << std::endl;
			exit(1);
		}
		if ( i > 0 && m.step[i] <= m.step[i-1] )
		{
			std::cerr << std::endl << "Error: screening ages must be strictly ascending (modality " << nr << ", " << m.age[i] << ")." << std::endl;
			exit(1);
		}
	}
}

std::mt19937 rng ;													// Mersenne Twister (32 bit) instead of rand()

double Random( void )
{
	return (double)rng() / 4294967296.0 ;							// uniform on [0,1)
}

void Healthy( int t )												// state = 0
{
	if ( rnd1[t] < p1[t] )
	{
		state = 1 ;
		Ndeath++ ;
	}
	else
	{
		if ( rnd1[t] < p1[t]+p2[t][0] )
		{
			state = 2 ;
			grade = 1 ;
			Ng[0]++ ;
		}
		else
		{
			if ( rnd1[t] < p1[t]+p2[t][0]+p2[t][1] )
			{	
				state = 2 ;
				grade = 2 ;
				Ng[1]++ ;
			}	
			else
			{
				if ( rnd1[t] < p1[t]+p2[t][0]+p2[t][1]+p2[t][2] )
				{	
					state = 2 ;
					grade = 3 ;
					Ng[2]++ ;
				}		
			}
		}
	} 
}

void DCIS( int t, int grade )										// state = 2
{
	int j = grade - 1 ;												
	if ( rnd1[t] < p1[t] )
	{
		state = 1 ;
		Ndeath++ ;
	}
	else
	{
		if ( rnd1[t] < p1[t]+p3[t][j] )
		{
			state = 0 ;
			grade = 0 ;
			Nreg++ ;
		}
		else
		{
			if ( rnd1[t] < p1[t]+p3[t][j]+p4[t][j] )
			{
				state = 3 ;
				NIBC[AgeGroup(t)][grade-1]++ ;
			}
		}
	}
}

// Screening with modality m (0 = modality 1, 1 = modality 2).
// Returns 4 (detection by modality 1) or 6 (detection by modality 2) on screen detection, otherwise 0.
int DoScreening( int t, int m )
{
	Modality &mo = M[m] ;
	double rndc = ( m == 0 ) ? rnd0[t] : rnd4[t] ;
	const double *rsens = ( m == 0 ) ? rnd2 : rnd5 ;
	const double *rspec = ( m == 0 ) ? rnd3 : rnd6 ;
	int result = 0 ;
	int j, k ;
	double disc_factor, discounted ;
	
	if ( !mo.active ) return 0 ;
	j = -1 ;
	for ( k = 0 ; k < (int)mo.step.size() ; k++ )
		if ( mo.step[k] == t )
			j = k ;
	if ( j >= 0 )													// t is a screening moment
	{
		if ( rndc < Compliance )										// Compliance
		{
			mo.Nmam++ ;
			disc_factor = DiscountFactor( t ) ;
			discounted = mo.cost_screen / disc_factor;
			mo.TotCostScreen += discounted ;
			mo.Nmami[j]++ ;
			if ( ( rspec[t] > mo.spec ) && ( state == 0 ) )          // check for FP in healthy state
			{
				mo.NFP++ ;
				discounted = mo.cost_biopsy / disc_factor;
				mo.TotCostBiopsy += discounted ;
			}
			else
			{
				if ( ( rsens[t] < mo.sens ) && ( state == 2 ) )			// sensitivity for detecting DCIS
				{
					discounted = mo.cost_biopsy / disc_factor;
					mo.TotCostBiopsy += discounted ;
					discounted = mo.cost_treatment / disc_factor;
					mo.TotCostTreatment += discounted ;
					result = ( m == 0 ) ? 4 : 6 ;
					mo.NSD[AgeGroup(t)][grade-1]++;
					mo.TotSD++ ;
				}
			}
		}	
	}
	return result ;
}

// Clinical detection: costs (biopsy + treatment of modality 1) and counts belong to file 1
int ClinicalDetection( int t ) 
{
	double disc_factor, discounted ;
	
	int result = 0 ;
	if (  ( rnd2[t] < ClinicalDet ) && ( state == 2 ) )				// ClinicalDet per 0.5 year
	{
		disc_factor = DiscountFactor( t ) ;
		discounted = M[0].cost_biopsy / disc_factor;
		M[0].TotCostBiopsy += discounted ;
		discounted = M[0].cost_treatment / disc_factor;
		M[0].TotCostTreatment += discounted ;
		result = 5 ;
		NCD[AgeGroup(t)][grade-1]++;
		TotCD++ ;
	}
	return result ;
}

// Unique seed per (iteration, participant, random stream): no overlap between participants or streams
unsigned int SeedFor( int pop, int it, int stream )
{
	return (unsigned int)( ( (unsigned long long)it * Npopulation + pop ) * NSTREAM + stream + 1 ) ;
}

void FillRandom( double *r, int pop, int it, int stream )
{
	int i ;
	rng.seed( SeedFor( pop, it, stream ) ) ;
	for ( i = 0 ; i < NSTEP ; i++ )
	   r[i] = Random() ;
}

void MakeRandomNumbers( int pop, int it )
{
	FillRandom( rnd1, pop, it, 0 ) ;
	FillRandom( rnd2, pop, it, 1 ) ;
	FillRandom( rnd3, pop, it, 2 ) ;
	FillRandom( rnd0, pop, it, 6 ) ;
	if ( M[1].active )												// extra streams for modality 2
	{
		FillRandom( rnd4, pop, it, 3 ) ;
		FillRandom( rnd5, pop, it, 4 ) ;
		FillRandom( rnd6, pop, it, 5 ) ;
	}
}

void InitParams( void )
{
	int i, j, m ;
	
	Ndeath = 0 ;
	Nreg = 0 ;
	
	for ( i = 0 ; i < 3 ; i++ ) 
		Ng[i] = 0 ;

	for ( i = 0 ; i < NAGEGRP ; i++ )
	for ( j = 0 ; j < 3 ; j++ )
	{
		NCD[i][j] = 0 ;
		NIBC[i][j] = 0 ;
		for ( m = 0 ; m < NMOD ; m++ )
			M[m].NSD[i][j] = 0 ;
	}
	TotCD = 0 ;
	
	TotDeathAge = 0 ;
	TotClindet = 0 ;
	for ( m = 0 ; m < NMOD ; m++ )
	{
		M[m].Nmam = 0 ;
		M[m].NFP = 0 ;
		M[m].TotSD = 0 ;
		M[m].TotDetected = 0 ;
		M[m].Nmami.assign( M[m].age.size(), 0 ) ;
		M[m].TotCostScreen = 0 ;
		M[m].TotCostBiopsy = 0 ;
		M[m].TotCostTreatment = 0 ;
	}
}

int main()
{
	int pop ;
	FILE *fp1, *fp2 , *fp3, *fp4, *fp5 ;
	int death, IBC, dum, detected, clindet ;
	int it, i, j ;
	char buf[12], fn[80] ;

	printf( "DCIS Simulation SimDCIS\n" ) ;
	printf( "Version 15.03\n" ) ;

	ReadTransitionProbabilities() ;
	ReadInputParams() ;
	SetupScreenSteps( M[0], 1 ) ;
	if ( M[1].active )
		SetupScreenSteps( M[1], 2 ) ;

	printf( "Compliance: %.2f\n", Compliance ) ;
	printf( "Clinical detection (per 0.5 year): %.4f\n", ClinicalDet ) ;
	printf( "Discount percentage, ref age: %.2f %.2f\n", cost_discount, cost_refage ) ;
	for ( j = 0 ; j < NMOD ; j++ )
	{
		if ( !M[j].active ) 
		{
			printf( "Modality %d: not used\n", j+1 ) ;
			continue ;
		}
		printf( "Modality %d\n", j+1 ) ;
		printf( "  Screening ages:" ) ;
		for ( i = 0 ; i < (int)M[j].age.size() ; i++ )
			printf( " %g", M[j].age[i] ) ;
		printf( "\n" ) ;
		printf( "  Sensitivity, Specificity: %.2f %.2f\n", M[j].sens, M[j].spec ) ;
		printf( "  Cost of screen, biopsy, treatment: %.2f %.2f %.2f\n", M[j].cost_screen, M[j].cost_biopsy, M[j].cost_treatment ) ;
	}
	printf( "Number of participants: %d\n", Npopulation ) ;
	printf( "Number of iterations: %d\n", Niterations ) ;
	printf( "simdcisx.out: for iteration x for every participant state per 0.5 year (4 = screen-detected mod 1, 6 = screen-detected mod 2)\n" ) ;
	printf( "simdcis1.det: statistics per iteration (modality 1 + clinical detection)\n" ) ;
	printf( "simdcis2.det: grades per age per iteration (modality 1 + clinical detection)\n" ) ;
	if ( M[1].active )
	{
		printf( "simdcis1_mod2.det: statistics per iteration, modality 2\n" ) ;
		printf( "simdcis2_mod2.det: grades per age per iteration, modality 2\n" ) ;
	}
	if ((fp2 = fopen( "simdcis1.det", "w"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open output file simdcis1.det" << std::endl;
		exit(1);
	}
	if ((fp3 = fopen( "simdcis2.det", "w"))==NULL) {
		std::cerr << std::endl << "Error: Cannot open output file simdcis2.det" << std::endl;
		exit(1);
	}
	fp4 = NULL ;
	fp5 = NULL ;
	if ( M[1].active )
	{
		if ((fp4 = fopen( "simdcis1_mod2.det", "w"))==NULL) {
			std::cerr << std::endl << "Error: Cannot open output file simdcis1_mod2.det" << std::endl;
			exit(1);
		}
		if ((fp5 = fopen( "simdcis2_mod2.det", "w"))==NULL) {
			std::cerr << std::endl << "Error: Cannot open output file simdcis2_mod2.det" << std::endl;
			exit(1);
		}
		fprintf( fp4, "Iter\tNmam\t" ) ;
		for ( i = 0 ; i < (int)M[1].age.size(); i++ )
			fprintf( fp4, "Nmam%g\t", M[1].age[i] ) ;
		fprintf( fp4, "NSD\tAvSDAge\tCostScreen\tCostBiopsy\tCostTreatment\tNFP\n" ) ;
		fprintf( fp5, "\tAge\tScreen\n" ) ;
		fprintf( fp5, "Iter\tGrade\tAll\t1\t2\t3\n" ) ;
	}
	fprintf( fp2, "Iter\tNmam\tNdeath\tNreg\t" ) ;
	for ( i = 0 ; i < (int)M[0].age.size(); i++ )
		fprintf( fp2, "Nmam%g\t", M[0].age[i] ) ;
	fprintf( fp2, "AvDeath\tAvSDAge\tAvCDAge\tCostScreen\tCostBiopsy\tCostTreatment\tNFP\n")	;
	fprintf( fp3, "\tAge\tScreen\t\t\t\tClinical\t\t\tProgressed\n" ) ;
	fprintf( fp3, "Iter\tGrade\tAll\t1\t2\t3\tAll\t1\t2\t3\tAll\t1\t2\t3\n" ) ;
	for ( it = 0 ; it < Niterations ; it++ )
	{
		sprintf( buf, "%d", it ) ;
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
			t = 0 ;
			state = 0 ;
			grade = 0 ;
			death = 0 ;
			IBC = 0 ;
			detected = 0 ;
			clindet = 0 ;
			while ( ( t < NSTEP ) && (!death) && (!IBC) && (state!=1) && (state!=3) && (state!=4) && (state!=5) && (state!=6) ) 
			{
				switch ( state )
				{
					case 0: Healthy( t ) ; break ;
					case 1: death = 1 ; break ;
					case 2: DCIS( t, grade ) ; break ;
					case 3: IBC = 1 ; break ;
					case 4: detected = 1 ; break ;
					case 5: clindet = 1 ; break ; 
				} 
				if ( (state!=1) && (state!=3) && (state!=4) && (state!=5) && (state!=6) )
				{
					dum = DoScreening( t, 0 ) ;						// modality 1
					if ( dum == 0 && M[1].active )
						dum = DoScreening( t, 1 ) ;					// modality 2
					if ( dum != 0 )
						state = dum ; 
					else
					{
						dum = ClinicalDetection( t ) ;
						if ( dum != 0 )
							state = dum ;
					}
				}
				if ( state == 1 )
				  TotDeathAge += t * DT ;								// age in years
				if ( state == 4 ) 
				  M[0].TotDetected += t * DT ;
				if ( state == 6 ) 
				  M[1].TotDetected += t * DT ;
				if ( state == 5 ) 
				  TotClindet += t * DT ; 
				fprintf( fp1, "%1d%1d\t", state, grade ) ;
				t++ ;
			}
			fprintf( fp1, "\n" ) ;
		}
		fclose( fp1 ) ;
		// --- file 1: modality 1 + clinical detection ---
		fprintf( fp2, "%d\t%d\t%d\t%d\t", it, M[0].Nmam, Ndeath, Nreg ) ;
		for ( i = 0 ; i < (int)M[0].age.size() ; i++ )
			fprintf( fp2, "%d\t", M[0].Nmami[i] ) ;
		fprintf( fp2, "%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", TotDeathAge / (float)Ndeath, (float)(M[0].TotDetected / M[0].TotSD), TotClindet / (float)TotCD, (float)M[0].TotCostScreen, (float)M[0].TotCostBiopsy, (float)M[0].TotCostTreatment, M[0].NFP ) ;
		for ( i = 0 ; i < NAGEGRP ; i++ )
		{
			fprintf( fp3, "%d\t%d\t", it, i*5 ) ;
			fprintf( fp3, "%d\t", M[0].NSD[i][0] + M[0].NSD[i][1] + M[0].NSD[i][2] ) ;
			for ( j = 0 ; j < 3 ; j++ )
				fprintf( fp3, "%d\t", M[0].NSD[i][j] ) ;
			fprintf( fp3, "%d\t", NCD[i][0] + NCD[i][1] + NCD[i][2] ) ;
			for ( j = 0 ; j < 3 ; j++ )
				fprintf( fp3, "%d\t", NCD[i][j] ) ;	
			fprintf( fp3, "%d\t", NIBC[i][0] + NIBC[i][1] + NIBC[i][2] ) ;
			for ( j = 0 ; j < 3 ; j++ )
				fprintf( fp3, "%d\t", NIBC[i][j] ) ;	
			fprintf( fp3, "\n" ) ;
		}
		// --- files modality 2 ---
		if ( M[1].active )
		{
			fprintf( fp4, "%d\t%d\t", it, M[1].Nmam ) ;
			for ( i = 0 ; i < (int)M[1].age.size() ; i++ )
				fprintf( fp4, "%d\t", M[1].Nmami[i] ) ;
			fprintf( fp4, "%d\t%.2f\t%.2f\t%.2f\t%.2f\t%d\n", M[1].TotSD, (float)(M[1].TotDetected / M[1].TotSD), (float)M[1].TotCostScreen, (float)M[1].TotCostBiopsy, (float)M[1].TotCostTreatment, M[1].NFP ) ;
			for ( i = 0 ; i < NAGEGRP ; i++ )
			{
				fprintf( fp5, "%d\t%d\t", it, i*5 ) ;
				fprintf( fp5, "%d\t", M[1].NSD[i][0] + M[1].NSD[i][1] + M[1].NSD[i][2] ) ;
				for ( j = 0 ; j < 3 ; j++ )
					fprintf( fp5, "%d\t", M[1].NSD[i][j] ) ;
				fprintf( fp5, "\n" ) ;
			}
		}
	}
	fclose( fp2 ) ;
	fclose( fp3 ) ;
	if ( fp4 != NULL ) fclose( fp4 ) ;
	if ( fp5 != NULL ) fclose( fp5 ) ;
}
