#include <stdio.h>
#include <math.h>

#include "semiD.h"

/*

A wrapper function for prospect which just returns the leaf
reflectance and transmitance at a single wavelength

Originally hacked up from Mat's code.

*/



int prospect_monochromatic( wavelength, vai, cab, cw, cp, cc, refl, tran )
float wavelength;
double vai, cab, cw, cp, cc ;
double *refl, *tran;
/*
vai = leaf structure
cab = leaf chlorophyll
cw = leaf water equivelent thickness
cp = protein concentration
cc = cellulose and lignin
*/
{

	double  main_r[NWBANDS], main_t[NWBANDS];

	int upper, lower ;
	float fraction ;
	
	if( wavelength < 400 || wavelength > 2500 ) return( 1 ) ;

	/* do the do (as Mat would say) */

	leaftwo_( main_r, main_t, &vai, &cab, &cw, &cp, &cc );

	/*now get the leaf optics for the desired wavelength*/
	
	upper = 1 + (wavelength-400) / 5.0 ;
	lower = (wavelength-400) / 5.0 ;
		
	fraction = (float) ( wavelength -400) / 5.0 - lower ;
	
	
	*tran = main_t[ lower ] * ( 1 - fraction ) + main_t[ upper ] * fraction ;
	*refl = main_r[ lower ] * ( 1 - fraction ) + main_r[ upper ] * fraction ;
	
	
	
	return( 0 );
}



int prospect_fullSpectrum_interp1nm( vai, cab, cw, cp, cc, refl, tran )
double vai, cab, cw, cp, cc ;
double *refl, *tran;
/*
vai = leaf structure
cab = leaf chlorophyll
cw = leaf water equivelent thickness
cp = protein concentration
cc = cellulose and lignin
*/
{

	double  main_r[NWBANDS], main_t[NWBANDS];
        long wavelength,w ;
	int upper, lower ;
	float fraction ;
	
	/* call prospect*/

	leaftwo_( main_r, main_t, &vai, &cab, &cw, &cp, &cc );

	/*now interpolate the output to 1nm*/
	
	for (wavelength=400;wavelength<=2500;wavelength++){

	  w=wavelength-400;
	
	  upper = 1 + w / 5.0 ;
	  lower = w / 5.0 ;
		
	  fraction = (float) w / 5.0 - lower ;
	
	
	  *(tran+w) = main_t[ lower ] * ( 1 - fraction ) + main_t[ upper ] * fraction ;
	  *(refl+w) = main_r[ lower ] * ( 1 - fraction ) + main_r[ upper ] * fraction ;
	
	}
	
	return( 0 );
}



/*
main( )
{

	int w;
	
	double r, t;
	
	for( w = 400; w<=2500; w++ ){
		prospect_monochromatic( w, 1.0, 75.0, 0.01, 0.001, 0.001, &r, &t );
		printf( "%d %lf %lf\n", w, r, t );
	}
	
	return( 1 ) ;

};
*/
