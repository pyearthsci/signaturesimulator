c**************************************************************
c  Name : nadimbrf.f
c***************************************************************
      SUBROUTINE nadimbrf(THETA_I,PHI_I,NV,THETA_V,PHI_V,
     * LAD,XRS,XHC,XLAI,RPL,XRL,XTL,BRF)
c
      INTEGER LAD,NV
      REAL BRF(1000)
      REAL THETA_V(1000),PHI_V(1000),THETA_I,PHI_I
      REAL XLAI,XHC,RPL,XRL,XTL,XRS
c***********************************************
c
c     Subroutine NADIMBRDF.f
c
c**********************************************
C
c  Output variable : BRF   = Bidirectional Reflectance Factor 
c
c 
c  Input variables :
c                       NV      = NUMBER OF VIEWING ANGLES
c                       THETA_I = SOLAR ZENITH ANGLE (IN DEGREES)
C                       PHI_I   = SOLAR AZIMUTH ANGLE (IN DEGREES)
c                       THETA_V = VIEWING ZENITH ANGLE (IN RADIANS)
c                       PHI_V   = VIEWING AZIMUTH ANGLE (IN RADIANS)
c
c                       LAD     = Leaf Angle Distribution
c                       = 1   <---> Planophile
c                       = 2   <---> Erectophile
c                       = 3   <---> Plagiophile
c                       = 4   <---> Extremophile
c                       = 5   <---> Uniforme
c
c                       XRS     = SOIL ALBEDO
c                       XHC     = HEIGHT OF THE CANOPY
C                       XLAI    = Leaf Area Index
c                       RPL     = radius of a single leaf
c                       XRL     = leaf reflectance
c                       XTL     = leaf transmittance
c
c************************************************************ 
c
      parameter (pi=3.141592653589793)
      real x_lambda_i,r1,r2,r3,xmeas
      real df,x_nf,a_f,n_c,h_c,x_ly,r 
      real c1
      real weights,points
      real Rs
      real tl,rl,lai
      real ag,bg,cg,dg
      integer number,ild
      real teta(1000),phi(1000)
      real teta_0,phi_0
c
c       COMMON DATAS
c
      common/canopee/df,x_nf,a_f,n_c,h_c,x_ly,r 
      common/hp/c1
      common/ga/weights(32),points(32)
      common/sol/Rs
      common/feuille/tl,rl,lai
      common/coef/ag,bg,cg,dg
      common/i/number,ild
      common/multi/I0(21,40),xIf(21,40),xI1u(21,40)
      common/limite/xI1,xImt
      common/angle_sol/teta_0,phi_0
c--------------------------------------------------
c
c       INPUT PARAMETERS
C
c       teta_0,phi_0            solar angles
c       ild                     leaf distribution
c       R_s                     soil albedo
c       h_c                     height of canopee (m)
c       df                      diameter of a single leaf (m)
c       lai                     leaf area index
c       rl                      leaf reflectance
c       tl                      leaf transmittance
c       na                      number of viewing angle
c       teta(i),phi(i)          viewing angles
c---------------------------------------------------------
c
c       INPUT PARAMETERS <------ ROUTINE
c
      lai = XLAI
      h_c = XHC
      df  = RPL*2.
      ild = lad
      rl  = XRL
      tl  = XTL
      Rs  = XRS
      na  = NV
c**************************************************
c
c       DEGREES ---> RADIANS
c
      Teta_0 = (180. - Theta_i )* pi / 180.
      Phi_0 =  Phi_i * pi / 180.
      DO I=1,na 
      teta(i) = THETA_V(i) * pi / 180.
      phi(i) = PHI_V(i) * pi / 180.
      ENDDO
c
c****************************************************
C
c       CALL SUBROUTINES WRITTING IN NADIMTOOLS.f
C
C***************************************************
c
c       16 GAUSS POINTS and WEIGHTS
c       to compute numerical integrals
c       and discrete ordinates methods
c
c
      number=16
      call gauleg(-1.,1.,points,weights,number)
c---------------------------------------------  
C
c      COEFFICIENTS FUNCTIONS OF DISTRIBUTION BUNNIK
C
c       (leaf angle distribution functions)
c
      call bunnik(ild)
c
C---------------------------------------------
C
c       GEOMETRIE CANOPY
C
      call archi(lai,teta_0)
c---------------------------------------------  
C
c       MULTIPLE INTENSITIES FOR ANGULAR POINTS GAUSS
c
c
      call multiple_dom(teta_0)
c
C-----------------------------------------------
      do i = 1 , na 
      x_lambda_i=0.01*abs(cos(teta_0))/G_ross(teta_0)*
     +         abs(cos(teta(i)))/G_ross(teta(i))
c
c   the three orders brdf   
c
      r1=rho_0_nad(teta(i),phi(i),x_lambda_i)
      r2=rho_1_nad(teta(i),phi(i),x_lambda_i)
      r3=rho_mult_nad(teta(i))
c
c so the sum ..to have the brdf total.
c
      xmeas= r3 + r1 +r2
      BRF(i)= xmeas
      enddo
      return
      end 



