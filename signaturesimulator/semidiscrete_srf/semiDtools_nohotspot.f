c**********************************************************
c nadimtools.f
c***********************************************************
c                                                          *
c This file (nadimtools.f) contains all routines and       * 
c functions used in the nadimbrf.f                         *    
c                                                          *    
c***********************************************************
c
c This function computes the "zero order" of scattering by the soil: 
c The downward radiation scattered once by the soil only.
C
      REAL FUNCTION rho_0_nad(teta_e,phi_e,x_lambda_i)
C
      parameter (pi=3.141592653589793)
      real x_lambda_i
      real teta_e,phi_e
      real Ki,Ke,xs1,xs2,rs_d
      real xh_p,xLi
      real c1
      real Rs
      real tl,rl,lai
      integer number,ild
      real teta_0,phi_0
c
      common/hp/c1
      common/sol/Rs
      common/feuille/tl,rl,lai
      common/i/number,ild
      common/angle_sol/teta_0,phi_0
c
      n_c=lai/x_lambda_i
      xg1=G_ross(teta_0)
      xg2=G_ross(teta_e)
      Ki=xg1/abs(cos(teta_0))
      Ke=xg2/cos(teta_e)
      xs1=(1.-x_lambda_i*Ki)**n_c
      xs2=1.
      xLi=c1/Geo(teta_e,teta_0,phi_e,phi_0)
c
c  actual optical path to account for the hot-spot effect
c
c      xh_p=hot_spot(Lai,xLi)
c
      xh_p=1.0
      do i=1,n_c
      xs2=xs2*(1.-x_lambda_i*Ke*xh_p)
      enddo
      rho_0_nad=Rs*xs2*xs1
      return
      end
c*****************************************************
c
c This function computes the "one order" of scattering by the leaves only
c
c 
c 
      REAL FUNCTION rho_1_nad(teta_e,phi_e,x_lambda_i)
c
      parameter (pi=3.141592653589793)
      real x_lambda_i
      real teta_e,phi_e
      real xLi
      real Ki,Ke,xga,xc1,rc_d
      real xg1,xg2,sum,x_hp,xL
      real c1
      real tl,rl,lai
      integer number,ild
      real teta_0,phi_0
c
      common/hp/c1
      common/feuille/tl,rl,lai
      common/i/number,ild
      common/angle_sol/teta_0,phi_0
c
      n_c=lai/x_lambda_i
      xg1=G_ross(teta_0)
      xg2=G_ross(teta_e)
      Ki=xg1/abs(cos(teta_0))
      Ke=xg2/cos(teta_e)
c
c  3D phase function
c
      xga=Gamma_leaf(teta_0,phi_0,teta_e,phi_e)
      xc1=(1.-x_lambda_i*Ki)
      xLi=c1/Geo(teta_e,teta_0,phi_e,phi_0)
      sum=0.
      do k=1,n_c
      xL=x_lambda_i*k
c      x_hp=hot_spot(xL,xLi)
      x_hp=1.0
      sum=sum+xc1**k*x_lambda_i*
     *(1.-x_lambda_i*Ke*x_hp)**k
      enddo
      rho_1_nad=sum/(cos(teta_e)*abs(cos(teta_0)))*xga
      return
      end
c***************************************************
c
c This function computes the multiple scattering term 
c averaged in azimuth. 
C       
      REAL FUNCTION rho_mult_nad(teta_u)
c
      parameter (pi=3.141592653589793)
      real xmu0,xm,xr,mu,Gu,xIm(21),m,dL
      real G0,teta_u,sum,sum1,sum2,x
      real xsu(21),xq0u(21),xq1(21),xgama_u(40)
      real fpar_tur,alb_tur,tr_tur,tr_dir,tran
      real weights,points
      real Rs
      real tl,rl,lai
      integer number,ild
      real I0,xIf,xI1u
      real xI1,xImt
      real teta_0,phi_0
c
      common/ga/weights(32),points(32)
      common/sol/Rs
      common/feuille/tl,rl,lai
      common/i/number,ild
      common/multi/I0(21,40),xIf(21,40),xI1u(21,40)
      common/limite/xI1,xImt
      common/angle_sol/teta_0,phi_0
c
c
      do i=1,21
      do j=1,20 
      enddo
      enddo
      xmu0=abs(cos(Teta_0))
      m=20
      dL=lai/m
      G0=G_ross(teta_0)
      xm=0.5*(1.-1.)
      xr=0.5*(1.+1.)
C
c       GAMMA (TETA(J) --> TETA_U)
c
      do j=1,number
      x=xm+xr*points(j)
      xgama_u(j)=fase_leaf(acos(x),teta_u)
      enddo
c
c
c       MULTIPLE SOURCE S(K) FOR VIEWING ANGLE = TETA_U
c           
c
      do k=1,m
      sum=0.
      do j=1,number
      x=xm+xr*points(j)
      sum=sum+2.*xgama_u(j)*xr*
     *weights(j)*(xIf(k+1,j)+xIf(k,j))/2.
      enddo
      xsu(k)=sum
      enddo
c
c
c       ZERO ORDER SOURCE
c
c
      do k=m,1,-1
      sum1=0.
      do j=number/2+1,number
      x=xm+xr*points(j)
      sum1=sum1+weights(j)*xr*2.*xgama_u(j)
     * *(I0(k+1,j)+I0(k,j))/2.
      enddo
      xq0u(k)=sum1
      enddo
c
c FIRST ORDER SOURCE
c
      do k=1,m
      sum2=0.
      do j=1,number
      x=xm+xr*points(j)
      sum2=sum2+2.*xgama_u(j)*xr*
     *weights(j)*(xI1u(k+1,j)+xI1u(k,j))/2.
      enddo
      xq1(k)=sum2
      enddo
c
      Gu= G_ross(teta_u)
      mu= cos(teta_u)
      xIm(m+1)=(xImt+xI1)
      do k=m,1,-1
      s1=xsu(k)+xq1(k)+xq0u(k)
      aa=Gu/2.-mu/dL
      bb=Gu/2.+mu/dL
      xIm(k)=(s1-xIm(k+1)*aa)/bb
      enddo
      rho_mult_nad=xIm(1)/(2.*abs(cos(teta_0)))
      return
      end
C***************************************************
c       
C FUNCTION G_ross = Ross function
C
      REAL FUNCTION G_ross(teta_p)
      parameter (pi=3.141592653589793)
      real weights,points
      integer number,ild
      common/ga/weights(32),points(32)
      common/i/number,ild
      real a,b,teta_p,ssc,xm,xr,x1,x2
      integer i
      ssc=0.
      a=0.
      b=pi/2.
      xm=0.5*(b+a)
      xr=0.5*(b-a)
      do 33 i=number/2+1,number
      dx=xr*points(i)
      x1=fg_ross(teta_p,xm+dx)
      x2=fg_ross(teta_p,xm-dx)
      ssc=ssc+weights(i)*(x1+x2)
 33   continue
      ssc=ssc*xr
      G_ross=ssc
      return
      end
c
C            ------------------------
c
      REAL FUNCTION fg_ross(teta_p,x)
      real teta_i,x
      teta_i=teta_p
      fg_ross=gl_bun(x)*psi_ross(teta_i,x)
      return
      end
c
c***************************************************
c
c       FONCTIONS BUNNIK (or gl_bun(teta_l)*sin(teta_l))
c    to compute the LAD 's functions.
C
      REAL FUNCTION gl_bun(x)
      parameter (pi=3.141592653589793)
      real ag,bg,cg,dg
      common/coef/ag,bg,cg,dg
      gl_bun=2./pi*(ag+bg*cos(2.*cg*x))+dg*sin(x)
      return
      end
c
c****************************************************
c
c Phase function.
c
      REAL FUNCTION Gamma_leaf(teta_p,phi_p,teta_i,phi_i)
      parameter (pi=3.141592653589793)
      real weights,points
      integer number,ild
      common/ga/weights(32),points(32)
      common/i/number,ild
      real teta_i,phi_i,teta_p,phi_p,gausg
      real tl,rl,lai
      common/feuille/tl,rl,lai
      real YM,YR,XM,XR,SY,dx,dy,sd
      xm=0.5*(pi/2.+0.)
      xr=0.5*(pi/2.-0.)
      ym=0.5*(2.*pi+0.)
      yr=0.5*(2.*pi-0.)
      gausg=0.
      do j=1,number
      dx=xm+xr*points(j)
      sy=0.
      do i=1,number
      dy=ym+yr*points(i)
      sd=fgamma_leaf(teta_p,phi_p,dx,dy,teta_i,phi_i)
      sy=sy+weights(i)*sd*xr
      enddo
      gausg=gausg+weights(j)*sy*yr
      enddo
      Gamma_leaf=gausg/2.
      return
      end
c
c             --------------------------
c
      REAL FUNCTION fgamma_leaf(teta_p,phi_p,x,y,teta_i,phi_i)
      parameter (pi=3.141592653589793)
      real weights,points
      real tl,rl,lai
      integer number,ild
      common/ga/weights(32),points(32)
      common/feuille/tl,rl,lai
      common/i/number,ild
      real f,g1,di,dp,dpp
      g1=gl_bun(x)
      dp=cos(x)*cos(teta_i)+sin(x)*sin(teta_i)*cos(abs(phi_i-y))
      dpp=cos(x)*cos(teta_p)+sin(x)*sin(teta_p)*cos(abs(phi_p-y))
      di=dp*dpp
      if (di.LT.0.) then
      f=rl
      else
      f=tl
      endif
      fgamma_leaf=g1*f/pi*abs(dp)*abs(dpp)
      return
      end
c
c               ------------------------------
c
      REAL FUNCTION psi_ross(teta_i,x)
      parameter (pi=3.141592653589793)
      real xmu,smu,xp,fit
      xmu=cos(x)
      smu=sin(x)
      if (xmu.EQ.1.) then
      psi_ross=cos(teta_i)
      else
      if (sin(teta_i).EQ.0.) then
      psi_ross=xmu
      else
      if (smu.eq.0.) then
      xp=0.
      else
      xp=1.*cos(teta_i)/sin(teta_i)*xmu/smu
      endif
      if (abs(xp).GT.1.) then
      psi_ross=cos(teta_i)*xmu
      else
      fit=acos(-xp)
      psi_ross=cos(teta_i)*xmu*(2.*fit/pi-1.)
     * +2./pi*sqrt(1.-cos(teta_i)**2)
     * *sqrt(1.-cos(x)**2)*sin(fit)
      endif
      endif
      endif
      psi_ross=abs(psi_ross)
      return
      end
c*****************************************************
c
c AZIMUTHALLY AVERAGED PHASE FUNCTION (after Shultis et al)
C
      REAL FUNCTION fase_leaf(teta_i,teta_e)
      parameter (pi=3.141592653589793)
      real weights,points
      real tl,rl,lai
      real teta_i,teta_e
      integer number,ild
      common/ga/weights(32),points(32)
      common/feuille/tl,rl,lai
      common/i/number,ild
      real xm,xr,sum
      xm=0.5*(pi/2.+0.)
      xr=0.5*(pi/2.-0.)
      sum=0.
      do j=1,number
      x=xm+xr*points(j)
      sum=sum+weights(j)*gl_bun(x)*
     * (tl*psi_leaf_p(teta_e,teta_i,x)+rl*psi_leaf_m(teta_e,teta_i,x))
      enddo
      fase_leaf=sum*xr
      return
      end
c
c----------------------------------------------------------
c
      REAL FUNCTION psi_leaf_p(teta_e,teta_i,teta_l)
      real ag,bg,cg,dg,xh1,xh2,xh3,xh4
      real teta_e,teta_i,teta_l,xmu_e,xmu_i,xmu_l
      common/coef/ag,bg,cg,dg
      xmu_e=cos(teta_e)
      xmu_i=cos(teta_i)
      xmu_l=cos(teta_l)
      xh1=xH_leaf(xmu_e,xmu_l)
      xh2=xH_leaf(xmu_i,xmu_l)
      xh3=xH_leaf(-xmu_e,xmu_l)
      xh4=xH_leaf(-xmu_i,xmu_l)
      psi_leaf_p=xh1*xh2+xh3*xh4
      return
      end
c
      REAL FUNCTION psi_leaf_m(teta_e,teta_i,teta_l)
      real ag,bg,cg,dg
      real teta_e,teta_i,teta_l,xmu_e,xmu_i,xmu_l
      common/coef/ag,bg,cg,dg
      xmu_e=cos(teta_e)
      xmu_i=cos(teta_i)
      xmu_l=cos(teta_l)
      xh1=xH_leaf(xmu_e,xmu_l)
      xh2=xH_leaf(-xmu_i,xmu_l)
      xh3=xH_leaf(-xmu_e,xmu_l)
      xh4=xH_leaf(xmu_i,xmu_l)
      psi_leaf_m=xh1*xh2+xh3*xh4
      return
      end
c
c       ------------------------------
c
      real function xH_leaf(xmu,xmu_l)
      parameter (pi=3.141592653589793)
      real ag,bg,cg,dg
      real x1,xmu_l,xmu,xfi,xh1,xh2,xh3
      common/coef/ag,bg,cg,dg
      if(xmu.eq.1.) then
      if(xmu_l.eq.1.) then
      xH_leaf=xmu*xmu_l
      return
      else
      if(xmu_l.eq.-1.) then
      xH_leaf=0.
      return
      else      
      if(xmu_l.gt.0.) then
      xH_leaf=xmu*xmu_l
      return            
      else
      xH_leaf=0.
      return
      endif
      endif
      endif
      else
      if(xmu_l.eq.1.) then
      if(xmu.eq.1.) then
      xH_leaf=xmu*xmu_l
      return
      else
      if(xmu.eq.-1.) then
      xH_leaf=0.
      return
      else
      if(xmu.gt.0.) then
      xH_leaf=xmu*xmu_l
      return
      else
      xH_leaf=0.
      endif
      endif
      endif
      else
      if(xmu.eq.-1.) then
      if (xmu_l.lt.0.) then
      xH_leaf=xmu*xmu_l
      return
      else
      xH_leaf=0.
      return
      endif
      else
      if(xmu_l.eq.-1.) then 
      if(xmu.gt.0.) then
      xH_leaf=0.
      return
      else
      xH_leaf=xmu*xmu_l
      return
      endif
      else
      x1=xmu*xmu_l/(sin(acos(xmu))*sin(acos(xmu_l)))
      if(x1.gt.1.) then
      xH_leaf=xmu*xmu_l
      return
      else
      if(x1.lt.-1.) then
      xH_leaf=0.
      return    
      else
      xfi=acos(-x1)
      xh1=xmu*xmu_l*xfi
      xh2=sqrt(1.-xmu**2.)*sqrt(1.-xmu_l**2.)
      xh3=sin(xfi)
      xH_leaf=abs(xh1+xh2*xh3)/pi
      return
      endif
      endif
      endif
      endif
      endif
      endif
      end
c
c********************************************************
c
c HOT_SPOT FUNCTION (after Verstraete et al)
C
      REAL FUNCTION hot_spot(x,xLi)
      parameter (pi=3.141592653589793)
      real x,xLi
      if (x.LT.xLi) then
      hot_spot=(1.-4./(3.*pi))*x/xLi
      else
      hot_spot=1.-4./(3.*pi)*xLi/x
      endif
      return
      end
c*******************************************************
c
      REAL FUNCTION Geo(teta_e,teta_0,phi_e,phi_0)
      real Li1,Li2
      real teta_e,teta_0,phi_0,phi_e
      real c1
      common/hp/c1
      Li1=tan(teta_0)**2+tan(teta_e)**2
      Li2=-2.*tan(teta_e)*tan(teta_0)*cos(phi_0-phi_e)
      Geo=sqrt(abs(Li1+Li2))
      if(Geo.lt.1.e-35) Geo=1.e-35
      return
      end
c
c*******************************************************
c
c HOT-SPOT DEEP
C
      REAL FUNCTION deep(teta_0)
      real teta_0,x_a
      x_a=0.005
      deep=-alog(x_a)*abs(cos(teta_0))/G_ross(teta_0)
      return
      end
c
c********************************************************
c
C RADIUS OF A SINGLE HOLE BETWEEN LEAVES 
C
      REAL FUNCTION sun_fleck(teta_0)
      parameter (pi=3.141592653589793)
      real teta_0,x_ly,x_nf,lai,x_l_opt
      integer n_c
      common/feuille/tl,rl,lai
      common/canopee/df,x_nf,a_f,n_c,h_c,x_ly,r
      x_l_opt=x_ly
      x_r2=abs(cos(teta_0))/(G_ross(teta_0)*x_l_opt)
      x_lambda_trou=-log(0.9)*abs(cos(teta_0))/G_ross(teta_0)
      x_trou=x_lambda_trou/a_f
      a_trou=(1.-x_trou*a_f*G_ross(teta_0)/abs(cos(teta_0)))/x_trou
      x_r1=sqrt(a_trou/pi)
      sun_fleck=sqrt(x_r2)*x_r1
      return
      end
c**********************************************************
c**********************************************************
c
c       SUBROUTINES
C
C*********************************************************
C*********************************************************
C
c       GEOMETRY OD THE VEGETATION CANOPY AND
C       COMPUTATION OF THE HOT-SPOT PARAMETER   
C
      SUBROUTINE ARCHI (x_lai,teta_s)
      parameter (pi=3.141592653589793)
      real x_lai,teta_s
      real weights,points
      real tl,rl,lai,Rs
      real ag,bg,cg,dg
      real teta_0,phi_0
      integer number,ild
      real c1
      real I0,xIf,xI1u,xI1,xImt
      real x_nf,a_f,h_c,x_ly,r
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
      a_f=(df/2.)**2*pi
      x_nf=x_lai/(a_f*h_c)
      x_ly=deep(teta_s)
      r=sun_fleck(teta_s)
      c1=2.*r*x_Lai/h_c
      return
      end
c
c*********************************************************
c
c       MULTIPLE SCATTERING INTENSITIES IN THE DIRECTION
C       CORRESPONDING TO THE GAUSS QUADRATURE.
C       
C 
      SUBROUTINE  multiple_dom(teta_0)
      parameter (pi=3.141592653589793)
      real weights,points
      real tl,rl,lai,Rs
      real ag,bg,cg,dg
      real teta_0
      integer number
      real I0,xIf,xI1u,xI1,xImt
      common/ga/weights(32),points(32)
      common/sol/Rs
      common/feuille/tl,rl,lai
      common/coef/ag,bg,cg,dg
      common/i/number,ild
      common/multi/I0(21,40),xIf(21,40),xI1u(21,40)
      common/limite/xI1,xImt
      real xmu0,G0,xm,xr,x
      real Gg(40),S(21,40)
      real Q0m(21,40),xI(21,40)
      real Q1(21,40),xgama_xy(40,40),xgama_0(40)
      real I00
      xmu0=abs(cos(Teta_0))
      m=20
      dL=lai/float(m)
      G0=G_ross(teta_0)
      xm=0.5*(1.+(-1.))
      xr=0.5*(1.-(-1.))
c
      I00=1.
c
c Computation of the G-function
C
      do j=1,number
      x=xm+xr*points(j)
      Gg(j)=G_ross(acos(x))
      enddo
c
C INITIALIZATION OF SOURCES AND INTENSITIES = 0
C
      do k=1,m
      do j=1,number
      S(k,j)=0.
      xIf(k,j)=0.
      enddo
      enddo
c
c COMPUTATION OF THE UNCOLLIDED INTENSITIES
c
      do k=1,m+1
      xL=float(k-1)*dL
      do i=1,number/2
      x=xm+xr*points(i)
      if((abs(x)).eq.xmu0) then
      I0(k,i)=I00*exp(-G0/xmu0*xL)
      else
      I0(k,i)=0.
      endif
      enddo     
      enddo
c
      do k=m+1,1,-1
      xL=float(k-1)*dL
      do i=number/2+1,number
      x=xm+xr*points(i)
      I0(k,i)=2.*Rs*I00*xmu0*exp(-G0/xmu0*lai)*
     *exp(-Gg(i)/x*(lai-xL))
      enddo
      enddo
c       
c TABLE OF 2-D PHASE FUNCTION FOR the GAUSS QUADRATURE ANGLES
c
      do i=1,number
      y=xm+xr*points(i)
      do j=1,number
      x=xm+xr*points(j)
      xgama_xy(i,j)=fase_leaf(acos(x),acos(y))
      enddo
      enddo
c
c COMPUTATION OF THE ZERO ORDER SOURCE
c
      do k=m,1,-1
      do i=1,number
      y=xm+xr*points(i)
      sum=0.
      do j=number/2+1,number
      x=xm+xr*points(j)
      sum=sum+weights(j)*xr*2.*xgama_xy(i,j)
     *  *(I0(k+1,j)+I0(k,j))/2.
      enddo
      Q0m(k,i)=sum  
      enddo
      enddo
c
c TABLE OF 2-D PHASE FUNCTION FOR ILLUMINATION --> QUADRATURE ANGLES
C       
      do i=1,number
      x=xm+xr*points(i)
      xgama_0(i)=fase_leaf(teta_0,acos(x))
      enddo
c
c COMPUTATIONS OF THE FIRST ORDER INTENSITIES
c
      do k=1,m+1
      xL=float(k-1)*dL
      do j=1,number/2
      x=xm+xr*points(j)
      if((abs(x)).ne.xmu0) then
      xI(k,j)=I00*2.*xgama_0(j)*xmu0*
     * (exp(-G0/xmu0*xL)
     * -exp(-Gg(j)/abs(x)*xL))/
     * (Gg(j)*xmu0-G0*abs(x))
      xI1u(k,j)=I00*2.*xgama_0(j)*xmu0*
     * (exp(-G0/xmu0*xL)
     * -exp(-Gg(j)/abs(x)*xL))/
     * (Gg(j)*xmu0-G0*abs(x))
      else
      xI(k,j)=I00*xL*2.*xgama_0(j)
     * *exp(-G0/xmu0*xL)/xmu0
      xI1u(k,j)=I00*xL*2.*xgama_0(j)
     * *exp(-G0/xmu0*xL)/xmu0
      endif
      enddo
      enddo     
c
      xI1=0.
      do j=1,number/2
      x=xm+xr*points(j)
      xI1=xI1+weights(j)*xr*abs(x)*2.*Rs*xI(m+1,j)
      enddo
c
      do k=m+1,1,-1
      xL=float(k-1)*dL
      do j=number/2+1,number
      x=xm+xr*points(j)
      xI(k,j)=I00*2.*xgama_0(j)*xmu0*
     * (exp(-G0/xmu0*xL)-
     * exp(-Gg(j)/x*(lai-xL))*exp(-G0/xmu0*lai))/
     * (G0*x+Gg(j)*xmu0)
c
      xI1u(k,j)=I00*2.*xgama_0(j)*xmu0*
     * (exp(-G0/xmu0*xL)-
     * exp(-Gg(j)/x*(lai-xL))*exp(-G0/xmu0*lai))/
     * (G0*x+Gg(j)*xmu0)
      enddo
      enddo
      do j=number/2+1,number
      xI(m+1,j)=0.
      xI1u(m+1,j)=0.
      enddo
c
      do k=1,m
      do i=1,number
      y=xm+xr*points(i)
      sum=0.
      do j=1,number
      x=xm+xr*points(j)
      sum=sum+weights(j)*xr*2.*xgama_xy(i,j)
     * *(xI(k+1,j)+xI(k,j))/2.
      enddo
      Q1(k,i)=sum
      enddo
      enddo
c
c MULTIPLE INTENSITIES
c
      do k=1,m+1
      do j=1,number/2
      xI(k,j)=0.
      enddo
      enddo
      l=0
  111 l=l+1
      do k=1,m
      do j=1,number/2
      x=xm+xr*points(j)
      xI(k+1,j)=(S(k,j)+Q0m(k,j)+Q1(k,j)-
     * xI(k,j)*(Gg(j)/2.+x/dL))/(Gg(j)/2.-x/dL)
      enddo
      enddo
      xImt=0.
      do j=1,number/2
      x=xm+xr*points(j)
      xImt=xImt+weights(j)*xr*2.*Rs*abs(x)*xI(m+1,j)
      enddo
      do j=number/2+1,number
      xI(m+1,j)=xImt+xI1
      enddo
      do k=m,1,-1
      do j=number/2+1,number
      x=xm+xr*points(j)
      xI(k,j)=(S(k,j)+Q0m(k,j)+Q1(k,j)-
     * xI(k+1,j)*(Gg(j)/2.-x/dL))/(Gg(j)/2.+x/dL)
      enddo
      enddo
c
      nt=0
      do k=1,m+1
      do j=1,number
      xnn=abs(xIf(k,j)-xI(k,j))
      if(xnn.lt.(1.e-4)) nt=nt+1
      xIf(k,j)=xI(k,j)
      enddo
      enddo
      if ((l.lt.1000).and.(nt.ne.(m+1)*number)) then
c
c MULTIPLE SOURCE 
c
      do k=1,m
      do i=1,number
      y=xm+xr*points(i)
      sum=0.
      do j=1,number
      x=xm+xr*points(j)
      sum=sum+weights(j)*xr*2.*xgama_xy(i,j)
     * *(xI(k+1,j)+xI(k,j))/2.
      enddo
      S(k,i)=sum
      enddo
      enddo
      goto 111
      endif
      return
      end
c
c******************************************************
c
c       ENERGY BALANCE
c
c
      SUBROUTINE energie(teta_i,phi_i,fpar,albedo_sys,trans_totale)
      parameter (pi=3.141592653589793)
c
      real weights,points
      real Rs
      real tl,rl,lai
      real ag,bg,cg,dg
      integer number,ild
      real I0,xIf,xI1u
      real xI1,xImt
      real teta_0,phi_0,teta_i,phi_i
      real fpar_tur,alb_tur,tran
c
      common/ga/weights(32),points(32)
      common/sol/Rs
      common/feuille/tl,rl,lai
      common/coef/ag,bg,cg,dg
      common/i/number,ild
      common/multi/I0(21,40),xIf(21,40),xI1u(21,40)
      common/limite/xI1,xImt
      common/angle_sol/teta_0,phi_0
      common/energie_tur/fpar_tur,alb_tur,tran
c
      real teta_s,phi_s,fpar,albedo_sys,trans_totale
      real Gs , xmus , xm , xr, ym ,yr
      real phi_w,teta_w,x,y,x_lambda(32,32),MM(32,32)
      real x_rho_0(32,32),x_rho_1(32,32),x_rho_m(32,32)
      real sum_y,sum_x,sum_k,gamma_0(32,32)
      real G_w(32),t_I1d(32,32)
      real tr_tur , tr_dir_tur
c
       teta_s = (180. - teta_i )* pi / 180.
       phi_s = (180.- phi_i)* pi / 180.
c
c BRF FOR THE GAUSS POINTS   
c
      Gs= G_ross (teta_s)
      xmus = cos(teta_s)
c    
      xm=0.5*(1.-1.)
      xr=0.5*(1.+1.)
      ym=0.5*(2.*pi-0.)
      yr=0.5*(2.*pi+0.)
c
c TABLE OF THE G_ross function
c Number of levels MM
c and BRF for all points of the Gauss quadrature
c
      do j=1,number
      x=xm+xr*points(j)
      teta_w=acos(x)
      G_w(j)= G_ross(teta_w)
      enddo 
c
      do i=1,number
      y=ym+yr*points(i)
      phi_w=y
      do j=1,number
      x=xm+xr*points(j)
      teta_w=acos(x)
      x_lambda(i,j)=0.01*abs(cos(teta_s))/Gs*
     *abs(cos(teta_w))/G_w(j)
      MM(i,j)=lai/x_lambda(i,j)
      x_rho_0(i,j)=rho_0_nad(teta_w,phi_w,x_lambda(i,j))  
      x_rho_1(i,j)=rho_1_nad(teta_w,phi_w,x_lambda(i,j))  
      x_rho_m(i,j)=rho_mult_nad(teta_w) 
      enddo
      enddo
c
c COMPUTATION OF THE ALBEDO
c
      sum_y=0.
      do i=1,number
      y=ym+yr*points(i)
      sum_x=0.
      do j=number/2+1,number 
      x=xm+xr*points(j)
      sum_x=sum_x+(x_rho_0(i,j)+x_rho_1(i,j)+x_rho_m(i,j))*
     *weights(j)*abs(x)*xr
      enddo
      sum_y = sum_y + sum_x*weights(i)*yr
      enddo     
      albedo_sys = sum_y  / (pi)
c
c DOWNWARD FIRST COLLIDED INTENSITIES WITH DISCRETE APPROACH   
c
      do i=1,number
      y = ym + yr * points (i)
      phi_w = y
      do j=1,number/2
      x = xm + xr * points (j)
      teta_w = acos(x)
      gamma_0(i,j) = Gamma_leaf (teta_s,phi_s,teta_w,phi_w)
      enddo
      enddo
c        
      do i=1,number
      y = ym + yr * points (i)
      phi_w = y
      do j=1,number/2
      x = xm + xr * points (j)
      teta_w = acos(x)
      sum_k = 0
      MMM=MM(i,j)
      do k=1,MMM
      sum_k = sum_k + (1.-x_lambda(i,j)*G_w(j)/abs(x))**(MMM-k)*
     *x_lambda(i,j)*
     *(1.- x_lambda(i,j)*Gs/abs(xmus)) **k 
      enddo 
      t_I1d (i,j) = sum_k * gamma_0 (i,j) * (1./ (pi * abs(x)))
      enddo 
      enddo 
c 
c TRANSMISSION OF THE FIRST COLLIDED INTENSITIES        
c
      sum_y=0.
      do i=1,number
      y=ym+yr*points(i)
      sum_x=0.
      do j=1,number/2
      x=xm+xr*points(j)
      sum_x=sum_x+ t_I1d(i,j)*
     +weights(j)*abs(x)*xr
      enddo
      sum_y = sum_y + sum_x * weights(i) * yr
      enddo     
      trans_1 = sum_y / abs(xmus) 
c
      sum_y=0.
      do i=1,number
      y=ym+yr*points(i)
      sum_x=0.
      do j=1,number/2
      x=xm+xr*points(j)
      sum_x=sum_x+ t_I1d(i,j)*rs*
     +weights(j)*abs(x)*xr
      enddo
      sum_y = sum_y + sum_x * weights(i) * yr
      enddo
      trans_1_in = sum_y / abs(xmus)
c
c TRANSMISSION OF MULTIPLY SCATTERED INTENSITIES
c
      sum_x=0.
      do j=1,number/2
      x=xm+xr*points(j)
      sum_x=sum_x+xIf(21,j)*
     *weights(j)*abs(x)*xr
      enddo
      trans_m = sum_x / ( abs (xmus) )
c
      sum_x=0.
      do j=number/2+1,number
      x=xm+xr*points(j)
      sum_x=sum_x+xImt*
     *weights(j)*abs(x)*xr
      enddo
      trans_in_m= sum_x / ( abs (xmus) )
      x_lambda_0 =  0.01 * *abs(xmus)/Gs*
     *abs(xmus)/Gs
      M0=lai/x_lambda_0
      xt_directe=(1.- x_lambda_0*Gs/abs(xmus))**M0 
      trans_totale= xt_directe + trans_m + trans_1
      fpar = 1.- albedo_sys - (1. - Rs) * trans_totale 
      return
      end
c
c**********************************************************
c
C       COEFFICIENTS FOR THE BUNNIK's FUNCTION
C
      subroutine bunnik(ild)
      real ag,bg,cg,dg
      common/coef/ag,bg,cg,dg
      ag=1.
      bg=1.
      cg=1.
      dg=1.
      if (ild.eq.1) then
      dg=0.
      endif
      if (ild.eq.2) then
      bg=-1.
      dg=0.
      endif
      if (ild.eq.3) then
      bg=-1.
      cg=2.
      dg=0.
      endif
      if (ild.eq.4) then
      cg=2.
      dg=0.
      endif
      if (ild.eq.5) then
      ag=0.
      bg=0.
      cg=0.
      endif
      return
      end
c
c****************************************************************
c 
c       CALCULATION OF THE QUADRATURE POINTS AND WEIGHTS (GAUSS)
c
      subroutine gauleg(x1,x2,x,w,n)
      parameter (pi=3.141592653589793)
      real weights,points
      real tl,rl,lai,Rs
      real ag,bg,cg,dg
      real teta_0,phi_0
      integer number,ild
      real c1
      real I0,xIf,xI1u,xI1,xImt
      real x_nf,a_f,h_c,x_ly,r
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
      integer n
      real*8 p1,p2,p3,pp,xm,z,z1
      real x1,x2,x(n),w(n),xl
      parameter (eps=3.d-14)
      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do 12 i=1,m
      z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1     continue
      p1=1.d0
      p2=0.d0
      do 11 j=1,n
      p3=p2
      p2=p1
      p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11    continue
      pp=n*(z*p1-p2)/(z*z-1.d0)
      z1=z
      z=z1-p1/pp
      if(abs(z-z1).gt.eps)go to 1
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
      w(n+1-i)=w(i)
12    continue
      return
      end
******************************************************************

