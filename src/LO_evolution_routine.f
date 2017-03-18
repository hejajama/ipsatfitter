* example of calling the gluon density evolved at LO in DGLAP
* with the evolution of the sinlget/gluon both coupled and de-coupled
C
C       implicit none
C       
C       double precision x,q2,cu2(15),gluon,gluondec(15,19),
C     1 gluoncou(15,51),xx(19)
C       integer i,j
C
C       data cu2/1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,20.,30.,40.,50.,60/
C       data xx/0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,
C     1 0.0009,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,
C     2 0.01/
C
C
C* TEST TO CHECK DIFFERENCES BETWEEN COUPLED AND DECOUPLED SINGLET/GLUON EVOLUTION
C
C       open(11,file='evolution_test.txt',status='unknown')
C              
C
C       call INIT		! initializes the Mellin momenta
C
C       do j=1,15
C       
C	       q2=cu2(j)
C
C    	   do i=1,19
C              
C		   x=xx(i)
C       
C    	   call LO_evol(x,q2,gluon,0)
C      
C	       gluondec(j,i)=gluon
C              
C    	   enddo
C
C       
C           do i=1,19
C              
C    	   x=xx(i)
C		   
C		   call LO_evol(x,q2,gluon,1)
C      
C           gluoncou(j,i)=gluon
C
C           enddo
C       
C* just write some output
C      
C           do i=1,19
C       
C           write(11,*) xx(i),gluondec(j,i),gluoncou(j,i)
C           write(6,*) xx(i),gluondec(j,i),gluoncou(j,i)
C
C           enddo
C       
C       enddo
C       
C       end

*************************  LO DGLAP evolution
** THE INVOKING ROUTINE MUST INCUDE
** 
** Call INIT: initializes the Mellin momenta
** 
** x : Bjorken momentum fraction
** Q2: virtuality
** coupling: 0 for decoupled evolution, 1 for coupled evolution
**
** This routine provides (x * gluon) density
**
** Any modification on the alpha_s parameters must be done here
**
** Modified by H.M.: This returns alphas(Q^2)*xg(x,Q^2)
*******************************************************************

       SUBROUTINE LO_evol(X, Q2, gluon, coupling, Ag, lambdag)
       
       implicit none
       
       DOUBLE PRECISION ALPHAS,C,EX,PA,WN(136),FR2,mur,ASMUR,MC,MB,MT,
     1 q2,fz,fun,x,ax,alpq,gluon,alps,alpc,alpb,alpt,Ag,lambdag
       double COMPLEX CC,XNM,CEX,N(136),FN(136)
       INTEGER NMAX, I1, M , coupling
       COMMON / CONTin  / C, CC                              
       COMMON / WEIGHTS / WN
       COMMON / MOMS    / N
       common / varalph / FR2, mur, ASMUR, MC, MB, MT  
       COMMON / FCOUPL  / ALPS, ALPC, ALPB,alpt

* THIS IS NEEDED FOR CALLING ALPHA_S, M_CHARM (MC) TO BE MODIFIED IF REQUIRED
* ALPHA_S VALUE AT THE INITIAL SCALE FROM LO MSTW2008

	   FR2 = 1.D0                ! ratio of mu_f^2 to mu_r^2c
	   MUR = 1.D0                ! input mu_r in GecV
 	   ASMUR = 0.68183d0         ! input value of alpha_s at mu_rc
	   MC = 1.27D0                ! charm quark mass
	   MB = 4.75D0               ! bottom quark mass
	   MT = 175d10               ! top quark mass

	   CALL INITALPHAS(0, FR2, mur, ASMUR, MC, MB, MT)

       ALPS = ALPHAS(mur)/(4.*3.141592654d0)
       ALPC = ALPHAS(mc)/(4.*3.141592654d0)
       ALPB = ALPHAS(mb)/(4.*3.141592654d0)
       ALPT = ALPHAS(mt)/(4.*3.141592654d0)
       ALPQ = ALPHAS(dSQRT(Q2))/(4.*3.141592654d0)

*...Q2- AND X-DEPENDENT QUANTITIES : 
        AX = dLOG (X)                                        
        EX = dEXP (-C * AX)                                 
*...INTEGRATION LENGTH PARAMETER FOR THE MELLIN INVERSION :     
             
        NMAX = 136                ! zmax = 36
                                             
*...CALCULATION OF THE GLUON DENSITY AND OUTPUT :

       CALL RENO (FN, N, ALPQ, NMAX, coupling, Ag, lambdag)

       FUN = 0.D0
       DO 1 I1 = 1, NMAX
       XNM = (C - N(I1)+1.) * AX
       CEX = zEXP (XNM) / 3.141592654d0 * CC
       FZ = dIMAG (FN(I1)*CEX)
       FUN = FUN + WN(I1) * FZ
  1    CONTINUE
       PA = FUN * EX

c NOTE: Here I have added *alphas!
       gluon=pa * ALPHAS(dSQRT(Q2))

       RETURN                           
       END

*...Mellin-n space Q**2 - evolution of the gluon at LO
*...The moments are calculated on an array of moments, N,  suitable 
*    for a (non-adaptive) Gauss quadrature. The necessary information
*    on amomalous dimensions, input and alpha(s)-values are taken from
*    the common-blocks ANOMS, INPUT and FCOUPL, respectively.
*
* Currently this takes the simplest possible fit form:
* xg = A_g x^(-lambdag) (1-x)^(5.6), following Amir&Raju
       SUBROUTINE RENO (FN, N, ALPQ, NMAX, coupling, Ag, lambdag)

       implicit none
       double complex ns8n,ns3n,gln,gl,sin,sg,svn,dsn,ssn,csn,bsn,tsn,
     1 ns35n,ns24n,ns15n,udbn,deln,usn,dvn,uvn,me2,ma2,me,ma,ep,epm,em,
     2 emp,oswi,ens,ANS(136,3:5),AM(136,3:5),AP(136,3:5),AL(136,3:5),
     3 BE(136,3:5),AB(136,3:5),AC(136,3:5),FN(136),N(136),XN,CBETA,one,
     4 two,zero
       Double precision XL,XL1,S,ALP,ALPS,ALPC,ALPB,ALPQ,alpt,
     1 B,b0,b1,alpq1,eqs3,q2,x,FR2,mur,ASMUR,MC,MB,MT,Ag,lambdag
       INTEGER F,K1,KK1,NMAX,coupling
       common/ varalph/ FR2, mur, ASMUR, MC, MB, MT         
       COMMON / ANOMS  / ANS, AM, AP, AL, BE, AB, AC
       COMMON / FCOUPL / ALPS, ALPC, ALPB,alpt

       DO 1 K1 = 1,NMAX  

*...INPUT MOMENTS OF THE PARTON DENSITIES :
*... AT THE LOW SCALE :
 
       XN=N(K1)
       
       zero=dcmplx(0.d0,0.d0)
       one=dcmplx(1.d0,0.d0)
       two=dcmplx(2.d0,0.d0)

c       GLN = 1.d0*(CBETA(XN+two,two+one))
c       SIN = 1.d0*(CBETA(XN+two,two+one))

        GLN = Ag * CBETA(XN - lambdag - one , 5.6+one)
c MSTW parameterization at LO
c       GLN = 0.0012216*(CBETA(XN-0.83657-one,2.3882+one)
c     1          -38.997*CBETA(XN-0.83657-0.5,2.3882+one)+
c     2           1445.5*CBETA(XN-0.83657,    2.3882+one))

       SIN = 0
c1.4335*(CBETA(XN+0.45232-one,3.0409+one)
c     1       -2.3737*CBETA(XN+0.45232-0.5,  3.0409+one)+
c     2        8.9924*CBETA(XN+0.45232,      3.0409+one))+
c     3       5.0903*(CBETA(XN+0.71978-one,5.1244+one)
c     4       -4.3654*CBETA(XN+0.71978-0.5,5.1244+one)+
c     5        7.4730*CBETA(XN+0.71978,    5.1244+one))+
c     6      0.59964*(CBETA(XN-0.16276-one,8.8801+one)
c     7       -2.9012*CBETA(XN-0.16276-0.5,8.8801+one)+
c     8        16.865*CBETA(XN-0.16276,    8.8801+one))
 

       if (ALPQ .GE. ALPC ) THEN ! EVOLUTION BELOW THE CHARM THRESHOLD :

       F = 3
       XL = ALPS / ALPQ
       B0 = 11.- 2./3.* F
       B1 = 102.- 38./3.* F
       B=B1/B0
       S   = dLOG (XL)
       XL1 = 1.- XL
       EM  = zEXP (-AM(K1,F)*S)
       EP  = zEXP (-AP(K1,F)*S)

       SG = SIN
       GL = GLN
      
       SIN = 0
cEM * ( AL(K1,F) * SG + BE(K1,F) * GL * coupling)
c     1     + EP * ( AC(K1,F) * SG - BE(K1,F) * GL * coupling)

       GLN = EM * ( AB(K1,F) * SG * coupling + AC(K1,F) * GL)
     1     + EP *( -AB(K1,F) * SG * coupling + AL(K1,F) * GL)

       else if ((ALPQ.LT.ALPC).and.(ALPQ.GE.ALPB)) then !between thresholds

       F = 3
       XL = ALPS / ALPC
       B0 = 11.- 2./3.* F
       B1 = 102.- 38./3.* F
       B=B1/B0
       S   = dLOG (XL)
       XL1 = 1.- XL
       EM  = zEXP (-AM(K1,F)*S)
       EP  = zEXP (-AP(K1,F)*S)

       SG = SIN
       GL = GLN

       SIN = 0
cEM * ( AL(K1,F) * SG + BE(K1,F) * GL * coupling)
c     1     + EP * ( AC(K1,F) * SG - BE(K1,F) * GL * coupling)

       GLN = EM * ( AB(K1,F) * SG * coupling + AC(K1,F) * GL)
     1     + EP *( -AB(K1,F) * SG * coupling + AL(K1,F) * GL)

       F = 4
       XL = ALPC / ALPQ

       B0 = 11.- 2./3.* F
       B1 = 102.- 38./3.* F
       B=B1/B0
       S   = dLOG (XL)
       XL1 = 1.- XL
       EM  = zEXP (-AM(K1,F)*S)
       EP  = zEXP (-AP(K1,F)*S)
       
       SG = SIN
       GL = GLN
      
       SIN = 0
cEM * ( AL(K1,F) * SG + BE(K1,F) * GL * coupling)
c     1     + EP * ( AC(K1,F) * SG - BE(K1,F) * GL * coupling)

       GLN = EM * ( AB(K1,F) * SG * coupling + AC(K1,F) * GL)
     1     + EP *( -AB(K1,F) * SG * coupling + AL(K1,F) * GL)

       elseif (ALPQ.LT.ALPB) then ! above bottom threshold

       F = 3
       XL = ALPS / ALPC
       B0 = 11.- 2./3.* F
       B1 = 102.- 38./3.* F
       B=B1/B0
       S   = dLOG (XL)
       XL1 = 1.- XL
       EM  = zEXP (-AM(K1,F)*S)
       EP  = zEXP (-AP(K1,F)*S)
       
       SG = SIN
       GL = GLN
      
       SIN = 0
cEM * ( AL(K1,F) * SG + BE(K1,F) * GL * coupling)
c     1     + EP * ( AC(K1,F) * SG - BE(K1,F) * GL * coupling)

       GLN = EM * ( AB(K1,F) * SG * coupling + AC(K1,F) * GL)
     1     + EP *( -AB(K1,F) * SG * coupling + AL(K1,F) * GL)

       F = 4
       ALP = ALPB
       XL = ALPC / ALPB
       B0 = 11.- 2./3.* F
       B1 = 102.- 38./3.* F
       B=B1/B0
       S   = dLOG (XL)
       XL1 = 1.- XL
       EM  = zEXP (-AM(K1,F)*S)
       EP  = zEXP (-AP(K1,F)*S)
       
       SG = SIN
       GL = GLN
      
       SIN = 0
cEM * ( AL(K1,F) * SG + BE(K1,F) * GL * coupling)
c     1     + EP * ( AC(K1,F) * SG - BE(K1,F) * GL * coupling)

       GLN = EM * ( AB(K1,F) * SG * coupling + AC(K1,F) * GL)
     1     + EP *( -AB(K1,F) * SG * coupling + AL(K1,F) * GL)

       F = 5
       XL = ALPB / ALPQ

       B0 = 11.- 2./3.* F
       B1 = 102.- 38./3.* F
       B=B1/B0
       S   = dLOG (XL)
       XL1 = 1.- XL
       EM  = zEXP (-AM(K1,F)*S)
       EP  = zEXP (-AP(K1,F)*S)
       
       SG = SIN
       GL = GLN
      
       SIN = 0
cEM * ( AL(K1,F) * SG + BE(K1,F) * GL * coupling)
c     1     + EP * ( AC(K1,F) * SG - BE(K1,F) * GL * coupling)

       GLN = EM * ( AB(K1,F) * SG * coupling + AC(K1,F) * GL)
     1     + EP *( -AB(K1,F) * SG * coupling + AL(K1,F) * GL)

       endif

       FN(K1) = GLN
             
   1   CONTINUE

       RETURN

       END
*
*
 
       double COMPLEX FUNCTION CBETA (Z1, Z2)

       implicit none
       double complex sub,zz1,z1,z2,zz2,lg1,lg2,lg12,lngam

       SUB = dCMPLX (0., 0.)
       ZZ1 = Z1
  1    CONTINUE
       IF ( dREAL (ZZ1) .LT. 15.) THEN
          SUB = SUB + zLOG ((ZZ1+Z2) / ZZ1)
          ZZ1 = ZZ1 + 1.
          GOTO 1
       END IF
       ZZ2 = Z2
  2    CONTINUE
       IF ( dREAL (ZZ2) .LT. 15.) THEN
          SUB = SUB + zLOG ((ZZ1+ZZ2) / ZZ2)
          ZZ2 = ZZ2 + 1.
          GOTO 2
       END IF
       LG1 = LNGAM (ZZ1)
       LG2 = LNGAM (ZZ2)
       LG12 = LNGAM (ZZ1 + ZZ2)
       CBETA = zEXP (LG1 + LG2 - LG12 + SUB) 
       RETURN
       END

***
       double complex function lngam (x)
 
       implicit none
       double complex x

       LNGAM = (X - dcmplx(0.5,0.0)) * zLOG (X) - X + 
     1  dcmplx(0.91893853,0.0) + dcmplx(1.0,0.0)/(dcmplx(12.0,0.0)* X)
     2  * (dcmplx(1.0,0.0)- dcmplx(1.0,0.0)/(dcmplx(30.,0.)* X*X) 
     3  * (dcmplx(1.,0.)- dcmplx(1.,0.0)/(dcmplx(3.5,0.0) * X*X)
     4  * (dcmplx(1.,0.0)- dcmplx(4.,0.0)/(dcmplx(3.,0.0)* X*X))))
       return
       end

 
*
*...Anomalous dimensions for leading order evolution of parton densities.
*   The moments are calculated on an externally given array of mellin 
*    moments, N, suitable for a (non-adaptive) quadrature.
*
*   Present version: the number of active flavours in the factorization
*    is fixed to ff=3, in the beta-function it varies between f=3 and 
*    f=5. The dimension of the moment array is 136.
*    
       SUBROUTINE ANOM (ANS, AM, AP, AL, BE, AB, AC, N)

       implicit none

       double complex QQI,QGF,GQI,GGI,GGF,XN,xac,xab,xbe,xal,xap,xam,
     1 xans,gm,gp,sq,gg,gq,qg,qq,ANS(136,3:5),AM(136,3:5),AP(136,3:5),
     2 AL(136,3:5),BE(136,3:5),AB(136,3:5),AC(136,3:5),N(136)
       Double precision B0, B1, B0F,b10,b02f,b02
       INTEGER F, FR, K1, K2, k3, k4, k5
       
       F = 3
       DO 1 K1 = 1, 136
       XN = N(K1)
       CALL ANCALC (QQI, QGF, GQI, GGI, GGF, XN)
       DO 2 K2 = 3, 5
       f=k2
*...ANOMALOUS DIMENSIONS AND RELATED QUANTITIES IN LEADING ORDER :
       B0 = 11.- 2./3.* F
       B0F = 11.- 2./3.* F
       B02 = 2.* B0
       B02F = 2.* B0F
       QQ = QQI
       QG = F * QGF
       GQ = GQI
       GG = GGI + F * GGF
       SQ = zSQRT ((GG - QQ) * (GG - QQ) + 4.* QG * GQ)
       GP = 0.5 * (QQ + GG + SQ)
       GM = 0.5 * (QQ + GG - SQ)
       XANS = QQ / B02
       XAM = GM / B02
       XAP = GP / B02
       XAL = (QQ - GP) / (GM - GP)
       XBE = QG / (GM - GP)
       XAB = GQ / (GM - GP)
       XAC = 1.- XAL
       ANS(K1,K2) = XANS
       AM(K1,K2) = XAM
       AP(K1,K2) = XAP
       AL(K1,K2) = XAL
       BE(K1,K2) = XBE
       AB(K1,K2) = XAB
       AC(K1,K2) = XAC
  2    CONTINUE

  1    continue
  
       RETURN
       END
*


*
       SUBROUTINE ANCALC (QQI, QGF, GQI, GGI, GGF, XN)

       implicit none 

       double complex ggf,ggi,gqi,qgf,qqi,cpsi,psifn,xnm,xn2,xn1,xns,xn
       Double precision ZETA2, ZETA3, CF, CA, TF
       XNS = XN * XN
       XN1 = XN + 1.
       XN2 = XN + 2.
       XNM = XN - 1.
       CF = 4./3.
       CA = 3.
       TF = 1./2.
*...LEADING ORDER :
       CPSI = PSIFN (XN1) + 0.577216
       QQI = (8./3.) * (-3.- 2./(XN * XN1) + 4.* CPSI)
       QGF = -4.* (XNS + XN +2.) / (XN * XN1 * XN2)
       GQI = -(16./3.) * (XNS + XN + 2.) / (XN * XN1 * XNM)
       GGI = -22.- 24./(XN * XNM) - 24./(XN1 * XN2) + 24.* CPSI
       GGF = 4./3.

       RETURN
       END

*
*
*...PSI - FUNCTION FOR COMPLEX ARGUMENT
*
       double COMPLEX FUNCTION PSIFN (Z)

       implicit none

       double COMPLEX Z, ZZ, RZ, DZ, SUB, zLOG
       SUB = dCMPLX (0.,0.)
       ZZ = Z
  1    CONTINUE
       IF (dREAL (ZZ) .LT. 10.) THEN
         SUB = SUB - 1./ ZZ
         ZZ = ZZ + 1.
         GOTO 1
       END IF
       RZ = 1./ ZZ
       DZ = RZ * RZ
       PSIFN = SUB + zLOG(ZZ) - RZ/2.- DZ/2520. * ( 210.+ DZ * (-21.+
     1         10.*DZ ))
       RETURN
       END
*

*
*...Initialization of support points in n-space and weights for the 
*    Gauss quadrature and of the anomalous dimensions for the RG 
*    evolution at these n-values. The outputs are written into the 
*    common-blocks CONT, MOMS, WEIGHTS and ANOMS, respectively. 
*    
       SUBROUTINE INIT

       implicit none
       
       double precision c,co,phi,si,sum,diff,z,ZS(8),WZ(8),DOWN(17),
     1 UP(17), WN(136)
       INTEGER I1, I2, I3, K
       double COMPLEX ANS(136,3:5),AM(136,3:5),AP(136,3:5),AL(136,3:5),
     1 BE(136,3:5),AB(136,3:5),AC(136,3:5),N(136),CC
       COMMON / WEIGHTS / WN
       COMMON / MOMS / N
       COMMON / CONTin / C, CC
       COMMON / ANOMS / ANS, AM, AP, AL, BE, AB, AC
*...WEIGHTS AND SUPPORT POINTS FOR NOMALIZED 8 POINT GAUSS QUADRATURE :
       DATA WZ /0.101228536290376,0.222381034453374,0.313706645877887,  
     1 0.362683783378362,0.362683783378362,0.313706645877887,
     2 0.222381034453374,0.101228536290376/
       DATA ZS/-0.960289856497536,-0.796666477413627,-0.525532409916329,
     1 -0.183434642495650,0.183434642495650,0.525532409916329,
     2 0.796666477413627,0.960289856497536/
*...INTEGRATION CONTOUR PARAMETERS :
       DATA DOWN / 0.D0, 0.5D0, 1.D0, 2.D0, 3.D0, 4.D0, 6.D0, 8.D0,
     1     1.D1, 1.2D1, 1.5D1, 1.8D1, 2.1D1, 2.4D1, 2.7D1, 3.D1, 3.3D1/
       C = 1
       PHI = 3.141592654 * 3./4.
       CO = DCOS (PHI)
       SI = DSIN (PHI)
       CC = dCMPLX (CO, SI)
       DO 1 I1 = 1, 16
         UP(I1) = DOWN(I1+1)
  1    CONTINUE
       UP(17) = 36.D0 
*...SUPPORT POINTS AND WEIGHTS FOR THE GAUSS INTEGRATION : 
*    (THE FACTOR (UP-DOWN)/2 IS INCLUDED IN THE WEIGHTS)
       K = 0
       DO 2 I2 = 1, 17
         SUM  = UP(I2) + DOWN(I2) 
         DIFF = UP(I2) - DOWN(I2) 
       DO 3 I3 = 1, 8
         K = K + 1
         Z = 0.5 * (SUM + DIFF * ZS(I3))
         WN(K) = DIFF / 2.* WZ(I3) 
         N(K)  = dCMPLX (C+CO*Z+1.,SI*Z)
  3    CONTINUE
  2    CONTINUE 
       CALL ANOM (ANS, AM, AP, AL, BE, AB, AC, N)
       RETURN
       END 

