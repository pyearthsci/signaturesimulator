#include <stdio.h>
#include <stdlib.h>
#include "semiD.h"

/*
Function to retreive interpolated soil spectra from a
set of principle components (which must be supplied)
Interpolates to 1nm.
*/
int price_soil_fullSpectrum_interp1nm( w1, w2, w3, w4, v1, v2, v3, v4, refl )
float w1, w2, w3, w4, *v1, *v2, *v3, *v4 ; 
double *refl ;
{


	int upper, lower ;
	float fraction ;
	long wavelength, w ;
	
	float rs_lower, rs_upper ;
	
	for (wavelength=400;wavelength<=2500;wavelength++){

	  w=wavelength-400;
	
	  upper = 1 + w / 5.0 ;
	  lower = w / 5.0 ;
		
	  fraction = (float) w / 5.0 - lower ;
	  rs_lower = w1 * *( v1 + lower ) + w2 * *( v2 + lower ) +  w3 * *( v3 + lower ) + w4 * *( v4 + lower ) ;
	  rs_upper = w1 * *( v1 + upper ) + w2 * *( v2 + upper ) +  w3 * *( v3 + upper ) + w4 * *( v4 + upper ) ;
	
	  *(refl+w) = rs_lower * ( 1 - fraction ) + rs_upper * fraction ;
	
	}
	
	return( 0 );
}

/*
Function to retreive interpolated soil reflectance from a
set of principle components (which must be supplied)
*/
int price_soil_monochromatic( wavelength, w1, w2, w3, w4, v1, v2, v3, v4, refl )
float wavelength;
float w1, w2, w3, w4, *v1, *v2, *v3, *v4, *refl ;
{


        int upper, lower ;
        float fraction ;
        
        float rs_lower, rs_upper ;
        
        if( wavelength < 400 || wavelength > 2500 ) return( 1 ) ;

        
        upper = 1 + (wavelength-400) / 5.0 ;
        lower = (wavelength-400) / 5.0 ;
                
        fraction = (float) ( wavelength - 400 ) / 5.0 - lower ;
        
        
        rs_lower = w1 * *( v1 + lower ) + w2 * *( v2 + lower ) +  w3 * *( v3 + lower ) + w4 * *( v4 + lower ) ;
        rs_upper = w1 * *( v1 + upper ) + w2 * *( v2 + upper ) +  w3 * *( v3 + upper ) + w4 * *( v4 + upper ) ;
        
        *refl = rs_lower * ( 1 - fraction ) + rs_upper * fraction ;
        
        
        return( 0 );
}



