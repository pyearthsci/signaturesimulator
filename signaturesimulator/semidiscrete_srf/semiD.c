/*
A front end to Nadine Gobron's CRM (semi discrete) that works on spectral respons functions.

It is coupled with PROSPECT to provide leaf spectra and it uses Price's EOFs to describe soil colour.

Nadine's original code has not been altered in anyway in
this model.
*/





/*
============================================
             Variable names
============================================

semiD:
-----

xrs = soil albedo
xhc = canopy height
xlai = LAI
rpl = leaf radius
xrl = leaf reflectance
xtl = leaf transmitance

lad = leaf angle distribution:
      1 - Planophile
      2 - Erectophile
      3 - Plagiophile
      4 - Extremophile
      5 - Uniform


PROSPECT:
--------

vai = leaf structure
cab = leaf chlorophyll
cw = leaf water equivelent thickness
cp = protein concentration
cc = cellulose and lignin


PRICE:
-----

soil_vector_N = the Nth spectral vector of the soil reflectance
rslN = the weight of the Nth vector


NB - the default values of the parameters for PROSPECT and PRICE
     are taken from Mat's kuusk code.
*/


#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <string.h>
#include <ctype.h>

#include "semiD.h"
#include "soil_rho.h"

#include<matrixio.h>


double *load_ffMatrix( long *xsize, long *ysize, char *filename )
/*Wrapper function for fread_ascii_ffMatrix that handles
the file opening and exits with failure if open fails*/
{
  FILE * fp ;
  double *ptr ;
        
  if( ( fp = fopen( filename, "r" ) ) == NULL ){
    fprintf( stderr, "%s: unable to open file %s\n", __FILE__, filename );
    exit( EXIT_FAILURE );
  }
        
  ptr = fread_ascii_ffMatrix( xsize, ysize, fp ) ;
  fclose( fp );
  return( ptr );
  
}




void usage( char *bin_name )
{


  fprintf( stderr, "usage: %s [options] < angles.dat\n\n", bin_name );
  fprintf( stderr, "Input data on the stdin is four columns of ascii:\nview_zenith view_azimuth solar_zenith solar azimuth\n\n"  );
  fprintf( stderr, "The options are:\n"  );
  fprintf( stderr, "-LAI arg\tset the leaf area index to arg\n"  );
  fprintf( stderr, "-lrad arg\tset the leaf radius to arg\n"  );
  fprintf( stderr, "-hc arg \tset the canopy height to arg\n"  );
  fprintf( stderr, "-lad arg\twhere arg is an int specifying the leaf angle distribution:\n"  );
  fprintf( stderr, "\t\t\t1 - Planophile\n"  );
  fprintf( stderr, "\t\t\t2 - Erectophile\n"  );
  fprintf( stderr, "\t\t\t3 - Plagiophile\n"  );
  fprintf( stderr, "\t\t\t4 - Extremeophile\n"  );
  fprintf( stderr, "\t\t\t5 - Uniform\n"  );  
  fprintf( stderr, "-N arg\tset the leaf structure variable to arg\n"  );
  fprintf( stderr, "-cab arg\tset the leaf chlorophyl concentration to arg\n"  );
  fprintf( stderr, "-cw arg \tset the equivelant leaf water thickness to arg\n"  );
  fprintf( stderr, "-cc arg \tset the celulose concentration to arg\n"  );
  fprintf( stderr, "-cp arg \tset the protein to arg\n"  );
  fprintf( stderr, "-rsl1 arg \tset the weight of the first soil vector to arg\n"  );
  fprintf( stderr, "-rsl2 arg \tset the weight of the first soil vector to arg\n"  );
  fprintf( stderr, "-rsl3 arg \tset the weight of the first soil vector to arg\n"  );
  fprintf( stderr, "-rsl4 arg \tset the weight of the first soil vector to arg\n"  );
  fprintf( stderr, "-u      \tprint this message\n"  );


}


double convolve(s1,s2,nw)
double *s1, *s2 ;
int    nw ;
{

  double out=0.0, sumS2=0.0;
  int i;
  
  for(i=0;i<nw;i++){
    out+=*(s1+i) * *(s2+i);
    sumS2+=*(s2+i);
  }
  return(out/sumS2);
}


int main( int argc, char **argv )
{

  void  cl_parser(  );

  char  line[ MAX_LINE_LEN ];
  int   nv=0,i, nw=NWAVELS, a, w,ww ;
  float theta_v[NANGLES], phi_v[NANGLES]; 
  float theta_i[NANGLES], phi_i[NANGLES]; 
  float brf[NANGLES][NWAVELS];
  float xrs=0.1, xhc=5.0, xlai=4.0, rpl=0.1, xrl=0.45, xtl=0.5;
  int   lad=2;

  /*leaf spectra file*/
  char    leaf_spectra_file_name[ MAX_LINE_LEN ];
  long    nColsLS, nRowsLS;
  double *leaf_spectra ;

  /*srf variables*/
  char  srf_file_name[ MAX_LINE_LEN ];
  long    nColsSrf, nRowsSrf, nSrf ;
  double *srf_data ;
  double  s1[NWAVELS], s2[NWAVELS];
  int     use_srf=0 ;

  /*for prospect*/
  double p_vai=1.0, p_cab=75.0, p_cw=0.01, p_cp=0.001, p_cc=0.001, p_r[2101], p_t[2101] ;
  int use_prospect=1;

  /*for price*/
  float *sv_1, *sv_2, *sv_3, *sv_4; 
  double p_s[2101] ;
  float rsl1=0.2, rsl2=0.1, rsl3=0.03726, rsl4=-0.002426 ;
  int use_price=1 ;
  
  
  int fast_mode=0;
  
  /*get command line options*/ 
  cl_parser( argc, argv, &xhc, &xlai, &xtl, &lad, &p_vai, &p_cab, &p_cw, &p_cp, &p_cc, &rsl1, &rsl2, &rsl3, &rsl4, srf_file_name, &use_srf, &fast_mode, leaf_spectra_file_name, &use_prospect );
  
        
  /*read the srf data*/
  if( use_srf ){
    srf_data = load_ffMatrix( &nColsSrf, &nRowsSrf, srf_file_name );
    if(nRowsSrf!=2101){
      fprintf(stderr,"Wrong number of samples in SRF file\n");
      exit(EXIT_FAILURE);
    }
  }
    
  /*read the angles in on the stdin*/
  while( fgets( line, MAX_LINE_LEN, stdin ) != NULL ){
  
    if( sscanf( line, "%f %f %f %f", &theta_v[nv], &phi_v[nv], &theta_i[nv], &phi_i[nv] ) != 4 ){
      fprintf( stderr, "%s: error on input, line %d\n", *argv, nv+1 );
      exit( EXIT_FAILURE );
    }
        
    /* zenith angles +ve */
    if ( theta_i[nv] < 0.0 ) {
      phi_i[nv] += 180 ;
      theta_i[nv] *= -1.0 ;
    }
    if ( theta_v[nv] < 0.0 ) {
      phi_v[nv] += 180 ;
      theta_v[nv] *= -1.0 ;
    }
    
    /* azimuth angles between 0 & 360 degrees */
    while ( phi_i[nv] > 360 ) phi_i[nv] -= 360;
    while ( phi_v[nv] > 360 ) phi_v[nv] -= 360;
    while ( phi_i[nv] < 0 ) phi_i[nv] += 360;
    while ( phi_v[nv] < 0 ) phi_v[nv] += 360;
    nv++;
  }


  if( use_prospect ){
    prospect_fullSpectrum_interp1nm( p_vai, p_cab, p_cw, p_cp, p_cc, p_r, p_t );  
  }  /*read the srf data*/
  else{
    leaf_spectra = load_ffMatrix( &nColsLS, &nRowsLS, leaf_spectra_file_name );
    if(nRowsLS!=2101||nColsLS!=2){
      fprintf(stderr,"Wrong number of samples in leaf spectra file\n");
      exit(EXIT_FAILURE);
    }
    for(w=0;w<2101;w++){
      p_r[w]=*(leaf_spectra+w*2);
      p_t[w]=*(leaf_spectra+w*2+1);
    }
  }


  if( use_price ){
    sv_1 = default_soil_vector_1 ;
    sv_2 = default_soil_vector_2 ;
    sv_3 = default_soil_vector_3 ;
    sv_4 = default_soil_vector_4 ;
    price_soil_fullSpectrum_interp1nm( rsl1, rsl2, rsl3, rsl4, sv_1, sv_2, sv_3, sv_4, p_s );    
  }    



  i=1;
  /*===================
  Loop over wavelengths
  =====================*/
  /*
  n.b. "fast" mode does the convolution on the leaf and soil spectra
  as opossed to doing monochromatic calculations of the canopy spectra 
  and then convolving against the SRF, which is very slow. The more
  narrow the SRF the more accurate the approximation.
  
  In fast mode nw is set to the number of sensor bands. Therefore "w"
  below becomes band number, rather than wavelength, which makes the code
  a little confusing to read.
  */
  
  if( fast_mode && use_srf ) nw=nColsSrf;
  for( w=0; w<nw; w++ ){    
    if( fast_mode && use_srf ){
      for( ww=0;ww<NWAVELS;ww++)
        s2[ww]=*(srf_data+w+nColsSrf*ww);
      xrs=convolve(p_s,s2,NWAVELS);
      xrl=convolve(p_r,s2,NWAVELS);
      xtl=convolve(p_t,s2,NWAVELS);
      //printf("** %f %f\n",xrl,xtl);
    }else{  
      xrs=p_s[w];
      xrl=p_r[w];
      xtl=p_t[w];
    }
    /*==============
    Loop over angles
    ================*/
    for( a=0; a<nv; a++ ){
      /*call the FORTRAN routine*/  
      //printf("**** %f %f %f %f %d %f %f %f %f %f %f\n", theta_i[a], phi_i[a], theta_v[a], phi_v[a], lad, xrs, xhc, xlai, rpl, xrl, xtl );
      nadimbrf_(&theta_i[a], &phi_i[a], &i , &theta_v[a], &phi_v[a], &lad, &xrs, &xhc, &xlai, &rpl, &xrl, &xtl, &brf[a][w]);
          
    }
  }


  /*SRF convolution *after* calculating full canopy spectra*/
  if( use_srf && fast_mode==0 ){
   for( a=0; a<nv; a++ ){ /*angles*/
     for(nSrf=0;nSrf<nColsSrf;nSrf++){ /*SRF bands*/
       /*
       Copy SRF data into arrays:
       */
       for( w=0;w<NWAVELS;w++){
         s1[w]=brf[a][w];
         s2[w]=*(srf_data+nSrf+nColsSrf*w);
       }
       xrs=convolve(s1,s2,NWAVELS);
       printf("%lf ",xrs);
     }
     printf("\n");
   }    
  
  /*this option is for when the SRF convolution is calculated
  on the leaf and soil spectra prior to */
  }else if( use_srf && fast_mode==1 ){
   for( a=0; a<nv; a++ ){ /*angles*/
     for(w=0;w<nColsSrf;w++)
       printf( "%f ", brf[a][w]);
     printf("\n");
   }    
  
  /*If not using SRF convolution print the whole spectra*/
  }else{
    /*print results*/
    for( w=0; w<NWAVELS; w++ ){ /*wavelengths*/
      printf("%d ",w+400);
      for( a=0; a<nv; a++ ){ /*angles*/
        printf( "%f ", brf[a][w]);
      }
      printf( "\n" );
    }
  }
  
  return( EXIT_SUCCESS ) ;

}




void cl_parser(  argc, argv, xhc, xlai, rpl, lad, p_vai, p_cab, p_cw, p_cp, p_cc, rsl1, rsl2, rsl3, rsl4, srf_file_name, use_srf, fast_mode, leaf_spectra_file_name, use_prospect )
int argc ;
char **argv ;
float *xhc, *xlai, *rpl ;
int *lad;
double *p_vai, *p_cab, *p_cw, *p_cp, *p_cc ;
float *rsl1, *rsl2, *rsl3, *rsl4;
char *srf_file_name ;
int  *use_srf, *fast_mode ;
char *leaf_spectra_file_name ;
int *use_prospect ;
{


  int i;
  void usage(  );
  
  for( i=1; i<argc; i++ ){


    if( *argv[ i ] == '-' ){
    
      /**/ if( !strncasecmp( argv[ i ], "-lad", 4 ) )   *lad     = atoi( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-LAI", 4 ) )   *xlai    = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-hc", 3 ) )    *xhc     = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-lrad", 3 ) )  *rpl     = atof( argv[ ++i ] ); 
            
      /*prospect*/
      else if( !strncasecmp( argv[ i ], "-N", 2 ) )   *p_vai     = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-cab", 4 ) ) *p_cab     = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-cw", 3 ) )  *p_cw      = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-cp", 3 ) )  *p_cp      = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-cc", 3 ) )  *p_cc      = atof( argv[ ++i ] ); 

      /*price*/
      else if( !strncasecmp( argv[ i ], "-rsl1", 5 ) )  *rsl1    = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-rsl2", 5 ) )  *rsl2    = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-rsl3", 5 ) )  *rsl3    = atof( argv[ ++i ] ); 
      else if( !strncasecmp( argv[ i ], "-rsl4", 5 ) )  *rsl4    = atof( argv[ ++i ] ); 
  
      
      /*srf file & fast srf mode*/
      else if( !strncasecmp( argv[ i ], "-fast", 5 ) )  *fast_mode = 1;
      else if( !strncasecmp( argv[ i ], "-srf", 4 ) ){
        strcpy( srf_file_name, argv[ ++i ] ) ; 
        *use_srf=1;
      }

      else if( !strncasecmp( argv[ i ], "-leaf_spectra", 8 ) ){
        strcpy( leaf_spectra_file_name, argv[ ++i ] ) ; 
        *use_prospect=0;
      }

      
      /*usage*/
      else if( !strncasecmp( argv[ i ], "-u", 2 ) ){  usage( argv[ 0 ] ); exit( EXIT_SUCCESS ) ; }
      
      else{
      
        fprintf( stderr, "%s: unknown option on command line: %s\n", argv[ 0 ], argv[ i ] );
        fprintf( stderr, "(use the option -u to see brief usage instructions)\n" );
        exit( EXIT_FAILURE );

      }

    }else{
    
      fprintf( stderr, "%s: unknown argument on command line: %s\n", argv[ 0 ], argv[ 1 ] );
      fprintf( stderr, "(use the option -u to see brief usage instructions)\n" );
      exit( EXIT_FAILURE );
    
    }

  }

  /*turn off fast mode if the user hasn't supplied an SRF file*/
  if(use_srf==0) fast_mode=0;

  return ;

}
