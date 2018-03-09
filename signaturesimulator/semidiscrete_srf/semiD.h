#ifndef NADIM_H
#define NADIM_H

/*do not change either of these in this version*/
#define NWBANDS 421 /*internal for prospect and price: 400-2500 nm at 5nm res*/
#define NWAVELS 2101 /*400-2500 nm at 1nm res*/

#define MAX_LINE_LEN 1000
#define NANGLES 10

void leaftwo_( double *, double *, double *, double *, double *, double *, double * );
void nadimbrf_( float *, float *, int *, float *, float *, int *, float *, float *, float *, float *, float *, float *, float * );
void nadimbrfe_( float *, float *, int *, float *, float *, int *, float *, float *, float *, float *, float *, 
												float *, float * , float *, float *,  float * );
void energie_( float *, float *,  float *, float *,  float *  );
int prospect_monochromatic( float, double, double, double, double, double, double *, double * );
int prospect_fullSpectrum_interp1nm( double, double, double, double, double, double *, double * );
int price_soil_monochromatic( float, float, float, float, float, float *, float *, float *, float *, float * );
int price_soil_fullSpectrum_interp1nm( float, float, float, float, float *, float *, float *, float *, double * );

#endif /*NADIM_H*/
