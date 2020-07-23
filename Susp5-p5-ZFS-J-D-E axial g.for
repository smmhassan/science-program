C     PROGRAM TO CALCULATE SUSCEPTIBILITY FOR two identical spin 1/2 with ZFS of triplet
c	 Assume axial symmetry for g.  Take gx=gy=g_perp and gz=g_para
Hamiltonian -S1*J*S2 + muH*g*(S1+S2) in which J principle axis is not oriented with field

C   VARIABLES J, D and Eta, g_perp, g_para and TIP
C		MUTLIFIELD PROGRAM
      PROGRAM SUSCEP
      IMPLICIT REAL*8 (A-H,O-Z)
c		Dimmensions of 2D arrays: Max number of fields max number of temps at each field
	common SUSDAT(10,100),TEMP(10,100),NTemps(10),Field(10),Param(6),P_Lambda(6),
	1 NFields,NXs,NThetas,thfix,NPhis
	REAL*8 BM(10,100)!,PMatrix(7,6),PRow(7),Ycolumn(6) !, R_MuEFF(200)
	real*8 funk 
      EXTERNAL FUNK    
	  !,SX(144),SZ(144)
	 CALL ZERO(SUSDAT,1000)
       CALL ZERO(TEMP,1000)
	 CALL ZERO(BM,1000)
	 CALL INTZERO(NTemps,10)
	 CALL ZERO(Field,10)

	 OPEN(1,FILE='H_T_Sus_dat.txt',ACTION='READ')!Expect 3 column data H,T,Xm
	  NFields=0
	  NXs=0
	  NThetas=0
	  NPhis=0
	  NParams=6
	  NumT=0
	  PreviousH=0
	  ReadH=0
	  ReadTemp=0
	  ReadSus=0

	  DO WHILE (.NOT. EOF(1))
		NXs = NXs + 1
		READ (1, *) ReadH,ReadTemp,ReadSus
        	If (ABS(ReadH-PreviousH).gt.0.1) then			 !field different by more than .1G is a field change
			NTemps(NFields)=NumT !Record # of T at previous field, puts zero in zeroth element.  
			NFields=NFields+1
			NumT=0 !reset count of temps at a particular H
			Field(NFields)=ReadH
			PreviousH=ReadH
		endif
		NumT=NumT+1
		Temp(NFields,NumT)=ReadTemp
		SusDat(NFields,NumT)=ReadSus 	  
	  END DO
	  NTemps(NFields)=NumT !Record # of Ts at the final field
	 CLOSE(1)

	 OPEN(2,FILE='H_T_Sus_Init.TXT',ACTION='READ')
	    READ(2,*) RJ,D,Eta,g_perp,g_para,TIP,NThetas,NPhis	!Field is not read from this file 
!	    READ(2,*) RJ,D,E,g_perp,g_para,TIP,NThetas,NPhis
!	    Eta=3.*E/D
		READ(2,*) dRJ,dD,dEta,dg_perp,dg_para,dTIP
		if(NThetas.eq.1)  then
			read(2,*) thfix
			thfix=thfix*3.14159/180.
		endif
	 CLOSE(2)
	 

       WRITE(*,1001)'J seed =',RJ,'   dJ =',dRJ
	 WRITE(*,1001)'D seed =',D,'    dD = ',dD
	 WRITE(*,1001) 'Eta seed =',Eta,'    dEta =',dEta
	 WRITE(*,1001) 'E seed =',D*Eta/3.
	 WRITE(*,1001)'g_perp seed',g_perp,'    dg_perp =',dg_perp
	 WRITE(*,1001)'g_para seed',g_para,'    dg_para =',dg_para
	 WRITE(*,1001) 'Temp Independent Paramagnetism',TIP,'    dTIP = ',dTIP
	 WRITE(*,*)'Number of thetas =',NThetas
	 WRITE(*,*)'Number of phis =',NPhis
	 WRITE(*,*)''

	 !g_factor = 2.D0
	 Param(1)=RJ				 !Coupling factor for Hamiltonian S1*J*S2
	 Param(2)=D
	 Param(3)=Eta
!	 Param(3)=E
	 Param(4)=g_perp
	 Param(5)=g_para
	 Param(6)=TIP	
	 P_lambda(1) = dRJ
	 P_lambda(2) = dD
	 P_lambda(3) = dEta
	 P_lambda(4) = dg_perp			   
	 P_lambda(5) = dg_para
	 P_lambda(6) = dTIP
	 RKtoWavenum = 0.695035	  !k/hc in CGS units converts Kelvin to cm-1
	  
 	 ftol=1e-4
	 iter=0
C		FOLLOWING CODE IMPLEMENTS CUSTOM MINIMIZE (SLOW!)
	 call minimize(Param,yerr,P_lambda,NParams,iter)
C	  FOLLOWING CODE IMPLEMENTS CONJUGATE GRADIENT MINIMIZATION	 
!	 CALL frprmn(Param,NParams,ftol,iter,fret)
	 write(*,*) 'Funk(Param)=',Funk(Param)
C		FOLLOWING CODE IMPLEMENTS AMOEBA
!	 DO 12 i=1,NParams  
!12		PMatrix(1,i)=Param(i)
!	 DO 30 k=1,4
!		DO 10 i=2,NParams+1
!			DO 10 j=1,NParams
!10				PMatrix(i,j)=PMatrix(1,j)
!		DO 20 i=1,NParams
!20			PMatrix(i+1,i)=PMatrix(i+1,i)+P_Lambda(i)/k 
!	    DO 4 m=1,NParams+1
!			DO 3 n=1,NParams
!				PRow(n)=PMatrix(m,n)
!3				write(*,*) PRow(n)
!	        YColumn(m)=funk(PRow)
!4			write(*,*) YColumn(m)
!30		call amoeba(PMatrix,YColumn,NParams+1,NParams,NParams,ftol,funk,iter)
!	 DO 40 i=1,NParams
!40		Param(i)=PMatrix(1,i)
!	 fret=funk(Param)
1001	  FORMAT(3X,A,F9.3,A,F9.3,A,F9.3) 
1002	  FORMAT(3X,A,ES15.3E2)
1003	  FORMAT(F10.0,F10.2,ES15.6E2,ES15.6E2)

	  OPEN(3,FILE='SUSRPRT.TXT',ACTION='WRITE')
		WRITE(3,*) 'MODEL: Identical 1/2-1/2 spins, triplet ZFS and axially symmetric g'
		WRITE(3,*) 'FINAL PARAMETERS THROUGH ITERATIONS: ',ITER
		WRITE(3,1002) 'FINAL ERROR',yerr
		WRITE(3,1001) 'J=',Param(1),' K  =',Param(1)*RKtoWavenum,' cm^-1'
		WRITE(3,1001) 'D=',Param(2),' K  =',Param(2)*RKtoWavenum,' cm^-1' 
		WRITE(3,1001) 'E=',Param(2)*Param(3)/3.,' K  =',Param(3)*RKtoWavenum*Param(2)/3.,' cm^-1' 
	    WRITE(3,1001) 'Eta=',Param(3)
		WRITE(3,1001) 'g perp= ',Param(4)
		WRITE(3,1001) 'g para= ',Param(5)
	    WRITE(3,1002) 'TIP=',Param(6)

!		WRITE(3,*) 'Ho',H0
	  CLOSE(3)
C	 Duplicate output
	  WRITE(*,*) 'MODEL: Identical 1/2-1/2 spins, ZFS and axially symmetric g'
	  WRITE(*,*) 'FINAL PARAMETERS THROUGH ITERATIONS: ',ITER
	  WRITE(*,1002) 'FINAL ERROR',yerr
	  WRITE(*,1001) 'J=',Param(1),' K  =',Param(1)*RKtoWavenum,' cm^-1'
	  WRITE(*,1001) 'D=',Param(2),' K  =',Param(2)*RKtoWavenum,' cm^-1' 
	  WRITE(*,1001) 'E=',Param(2)*Param(3)/3.,' K  =',Param(3)*RKtoWavenum*Param(2)/3.,' cm^-1' 
	  WRITE(*,1001) 'Eta=',Param(3)
	  WRITE(*,1001) 'g perp= ',Param(4)
	  WRITE(*,1001) 'g para= ',Param(5)
	  WRITE(*,1002) 'TIP=',Param(6)
	  WRITE(*,*) ''

	  WRITE(*,*) ''
	  WRITE(*,*)  'ZERO FIELD ENERGIES'
	  WRITE(*,1001) 'E10 = -J/4-2D/3=',-Param(1)/4.-2.*Param(2)/3.
	  WRITE(*,1001) 'E00 = 3J/4=',3.*Param(1)/4.
	  WRITE(*,1001) 'E11 = -J/4+D/3-E=',-Param(1)/4.+(1-Param(3))*Param(2)/3.
	  WRITE(*,1001) 'E1-1 = -J/4+D/3+E=',-Param(1)/4.+(1+Param(3))*Param(2)/3.
!	  WRITE(*,*) 'E11 = J/4+D/3-E=',Param(1)/4.+Param(2)/3.-Param(3)
!	  WRITE(*,*) 'E1-1 = J/4+D/3+E=',Param(1)/4.+Param(2)/3.+Param(3)
	  WRITE(*,*) ''
	  WRITE(*,*) ''
!	  WRITE(*,*) 'Ho',H0
!	  WRITE(*,*)  'Temp','		Data','		Fit'
	  RJ=Param(1)
	  D=Param(2)
	  E=Param(3)*Param(2)/3.
! 	  E=Param(3)
	  g_perp=Param(4)
	  g_para=Param(5)
	  TIP=Param(6)

	  CALL SUSCAL(RJ,D,E,g_perp,g_para,TIP,BM,NFields,Field,NTemps,TEMP,NThetas,thfix,NPhis)
	  OPEN(4,FILE='BMOUT.TXT',ACTION='WRITE')
	  DO 50 JField = 1,NFields
     		DO 50 JTemp=1,NTemps(JField)
!		   WRITE(*,1003) Field(JField),TEMP(JField,JTemp),SUSDAT(JField,JTemp),BM(JField,JTemp)
  		   WRITE(4,1003) Field(JField),TEMP(JField,JTemp),SUSDAT(JField,JTemp),BM(JField,JTemp)
  50		CONTINUE
	  CLOSE(4)
	 END

	 FUNCTION FUNK(ParamTest)
        IMPLICIT REAL*8 (A-H,O-Z)
	  common SUSDAT(10,100),TEMP(10,100),NTemps(10),Field(10),Param(6),P_Lambda(6),
	1  NFields,NXs,NThetas,thfix,NPhis
     	  DIMENSION BMSIM(10,100),ParamTest(6)
	  REAL*8 FUNK

	  RJ=ParamTest(1)
	  D=ParamTest(2)
	  E=ParamTest(3)*ParamTest(2)/3.
!	  E=ParamTest(3)
	  g_perp=ParamTest(4)
	  g_para=ParamTest(5)
	  TIP=ParamTest(6)
	  CALL SUSCAL(RJ,D,E,g_perp,g_para,TIP,BMSIM,NFields,Field,NTemps,TEMP,NThetas,thfix,NPhis)
	  FUNK = CALCERR(BMSIM,SUSDAT,NTemps,NFields)
	  RETURN 
	 END

C	 Calculate partial derivatives of funk
	 SUBROUTINE DFUNK(P,DP)
       IMPLICIT REAL*8 (A-H,O-Z)
	 COMMON SUSDAT(10,100),TEMP(10,100),NTemps(10),Field(10),Param(6),P_Lambda(6),
	1  NFields,NXs,NThetas,thfix,NPhis
       DIMENSION P(6),DP(6)
	 NofParams=6
	 DO 10 I=1,NofParams
		IF (P_lambda(i).ne.0.0) THEN
			PSAVE=P(I)
			P(I)=P(I)+P_LAMBDA(I)
			X1=FUNK(P)
			P(I)=PSAVE-P_LAMBDA(I)
			X2=FUNK(P)
			P(I)=PSAVE
			DP(I)=(X1-X2)/(2.0*P_LAMBDA(I))
		ELSE
			dp(i)=0.0
		ENDIF
10	 CONTINUE
	 RETURN
	 END


	 FUNCTION CALCERR(BM,SUSDAT,NTemps,NFields)
        IMPLICIT REAL*8 (A-H,O-Z)
	  DIMENSION SUSDAT(10,100),BM(10,100),NTemps(10)
	  ERR=0
C	  !denominator of error ratio.  CalcErr should be sum(data-fit)^2/sum(data^2
	  SumSusDatSqrd=0
	  DO 10 JField=1,NFields
		DO 10 JTemp=1,NTemps(JField)
			SumSusDatSqrd=SumSusDatSqrd+SUSDAT(JField,JTemp)**2
10			ERR=ERR+(SUSDAT(JField,JTemp)-BM(JField,JTemp))**2
	  CALCERR=ERR/SumSusDatSqrd
	  RETURN 
	 END

      SUBROUTINE SUSCAL(RJ,D,E,g_perp,g_para,TIP,BM,NFields,Field,NTemps,TEMP,NThetas,thfix,NPhis)
       IMPLICIT REAL*8 (A-H,O-Z)
	 DIMENSION HR(4,4),HI(4,4),UR(4,4),UI(4,4),ENG(4),Field(10),NTemps(10)
       DIMENSION ENG2(4),SUS(10,100),BM(10,100),TEMP(10,100)  !,SX(144),SZ(144)
		
	 CALL ZERO(SUS,1000)
C		Identical Spin 1/2-1/2 system with 4 eigenstates
	 NStates=4
	 IFVEC=0  !flag that says don't want eigenvectors
	 DO 100 JField=1,NFields
!	 write(*,*) 'JField=',JField,'NFields=',NFields
		H0=Field(JField)
		dH=0.01*H0		!Try for approriate scale over range of H0
		H1=H0+dH
		Havg=(H0+H1)/2
		UBOHR=6.7171D-5	!Units K/Oe  Divided by Boltzman's constant k
		RAD=.017453
		Nav_Boltk=8.3115D+7  !Boltzman's constant times Avagadro's number
		Pi=3.141592654
		!dOmega=sin(theta)dThetaDPhi = d(Cos(theta))dPhi
		DO 110 JAngle=1,NThetas
			CosTH=(JAngle-0.5)/NThetas
			THETA=DACOS(CosTH)
			if(NThetas.eq.1)theta=thfix
			DO 120 KAngle=1,NPhis
				dPhi=Pi/(2*NPhis)  !theta->0-90  PHI->0-90 for 1st octant
	            PHI=(KAngle-0.5)*dPhi
				CALL HSET(RJ,D,E,g_perp,g_para,TIP,H0,THETA,PHI,HR,HI)
!				CALL DIAG(HR,ENG,UR,NStates)
				CALL DIAG(HR,HI,ENG,UR,UI,NStates,IFVEC) !What's IFVEC
				CALL HSET(RJ,D,E,g_perp,g_para,TIP,H1,THETA,PHI,HR,HI)
!				CALL DIAG(HR,ENG2,UR,NStates)
				CALL DIAG(HR,HI,ENG2,UR,UI,NStates,IFVEC) !What's IFVEC

				DO 130 JTemp=1,NTemps(JField)
					T=TEMP(JField,JTemp)
					Z=0.0
					Z1=0.0
					E1=ENG2(1)
					DO 10 K=1,NStates
						ARG=(ENG2(K)-E1)/T
						IF (ARG.GT.80.0) GOTO 10
						Z1=Z1+DEXP(-ARG)
						Z=Z+DEXP(-(ENG(K)-E1)/T)
   10					CONTINUE
C				Sus returns system molar susceptability 
				SUS(JField,JTemp)=SUS(JField,JTemp)+Nav_Boltk*T*(1/Havg)*(DLOG(Z1)-DLOG(Z))/dH
  130				CONTINUE !Temp loop
  120			CONTINUE	 !PHI loop
  110		CONTINUE !THETA loop
  100  CONTINUE  !Field loop
       
	 
	 DO 40 KField=1,NFields
		DO 40 KTemp=1,NTemps(KField)
			T=TEMP(KField,KTemp)
			SUS(KField,KTemp)=SUS(KField,KTemp)/(NThetas*NPhis)+TIP !scale molar susceptibility	model+TIP
C 40				BM(L)=2.828*DSQRT(T*SUS(L)/2) !Returns effective moment/ion of dimer"/2"
C 40				BM(L)=H0*SUS(L)/UBOHR !Returns average moment(BM) to fit
 40			BM(KField,KTemp)=SUS(KField,KTemp) !Returns molar susceptability to fit against data
	  RETURN
      END
C
C	Hamiltonian -S1*J*S2 + muH*g*(S1+S2) in which J principle axis is not oriented with field
C	Will assume axial symmetry for g: gx=gy=g_perp and gz=g_para for now.  
       SUBROUTINE HSET(RJ,D,E,g_perp,g_para,TIP,H,THETA,PHI,HR,HI)  !HR is real part of hamiltonian, HI is imaginary part
       IMPLICIT REAL*8 (A-H,O-Z)
       DIMENSION HR(4,4),HI(4,4)
       CALL ZERO(HR,16)
	 CALL ZERO(HI,16)
       UBOHR=6.7171D-5	  !Units K/Oe  Divided by Boltzman's constant k
	 gx=g_perp
	 gy=g_perp 
	 gz=g_para
!	 gy=g_perp
!	 gz=g_para
C		Hamiltonian in units of K
!	 WRITE(*,*) RJ,g_factor,H,Theta
	 Hx=gx*UBOHR*H*DSIN(THETA)*DCOS(PHI)
	 Hy=gy*UBOHR*H*DSIN(THETA)*DSIN(PHI)
       Hz=gz*UBOHR*H*DCOS(THETA)

C	SET DIAGONAL ELEMENTS OF HR
	 HR(1,1)=D/3.-RJ/4. + Hz
	 HR(2,2)=-(D/3.-RJ/4.)
	 HR(3,3)=-(D/3.-RJ/4.)
	 HR(4,4)=D/3.-RJ/4. - Hz
C	DIAGONAL ELEMENTS OF IR	are all zero
C	SET OFF-DIAGONAL ELEMENTS  OF HR
	 HR(1,2)=Hx/2.
	 HR(2,1)=Hx/2.
	 HR(1,3)=Hx/2.
	 HR(3,1)=Hx/2.
	 HR(1,4)=E
	 HR(4,1)=E
	 HR(2,3)=-RJ/2.-D/3.
	 HR(3,2)=-RJ/2.-D/3.
	 HR(2,4)=Hx/2.
	 HR(4,2)=Hx/2.
	 HR(3,4)=Hx/2.
	 HR(4,3)=Hx/2.
C	 SET OFF-DIAGONAL ELEMENTS OF IR
	 HI(1,2)=-Hy/2.
	 HI(2,1)=Hy/2.
	 HI(1,3)=-Hy/2.
	 HI(3,1)=Hy/2.
	 HI(1,4)=0.
	 HI(4,1)=0.
	 HI(2,3)=0.
	 HI(3,2)=0.
	 HI(2,4)=-Hy/2.
	 HI(4,2)=Hy/2.
	 HI(3,4)=-Hy/2.
	 HI(4,3)=Hy/2.
       RETURN
       END
C
	SUBROUTINE DIAG(HR,HI,ENG,UR,UI,N,IFVEC)
      IMPLICIT REAL*8 (A-H,O-Z)
	DIMENSION HR(N,N),HI(N,N),ENG(N),UR(N,N),UI(N,N),WK(N),TAU(2,N)
C     CONVERT TO REAL TRIDIAGONAL FORM
      CALL HTRIDI(N,N,HR,HI,ENG,WK,WK,TAU)
      N2=N*N
      CALL ZERO(UR,N2)
      DO 50 I=1,N
   50 UR(I,I)=1.
C     CALCULATE EIGENSYSTEM OF TRIDIAG MATRIX
      CALL TQL(N,N,ENG,WK,UR,IERR,IFVEC)
C     IF(IERR.EQ.0)GO TO 60
C      WRITE(*,1) IERR
C   1 FORMAT(' ERROR IN DIAG, IERR=',I2)
C     CONVERT VECTORS TO THOSE OF ORIGINAL MATRIX 
   60 IF(IFVEC.NE.0)CALL HTRIBK(N,N,HR,HI,TAU,N,UR,UI)
      RETURN
      END 
c      
      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
      TAU(1,N)=1.
      TAU(2,N)=0.
      DO 100 I=1,N
  100		D(I)=AR(I,I)
      DO 300 II=1,N
		I=N+1-II
		L=I-1
		H=0.
		SCALE=0.
		IF(L.LT.1)GO TO 130
C     SCALE ROW
		DO 120 K=1,L
  120			SCALE=SCALE+ ABS(AR(I,K))+ ABS(AI(I,K))
		IF(SCALE.NE.0.)GO TO 140
		TAU(1,L)=1.
		TAU(2,L)=0.
  130		CONTINUE
		E(I)=0.
		E2(I)=0.
		GO TO 290
  140		CONTINUE
		DO 150 K=1,L
			AR(I,K)=AR(I,K)/SCALE
			AI(I,K)=AI(I,K)/SCALE
			H=H+AR(I,K)*AR(I,K)+AI(I,K)*AI(I,K)
  150		CONTINUE
		E2(I)=SCALE*SCALE*H
		G= SQRT(H)
		E(I)=SCALE*G
		F=CABS(CMPLX(AR(I,L),AI(I,L)))
		IF(F.EQ.0.)GO TO 160
		TAU(1,L)=(AI(I,L)*TAU(2,I)-AR(I,L)*TAU(1,I))/F
		SI=(AR(I,L)*TAU(2,I)+AI(I,L)*TAU(1,I))/F
		H=H+F*G
		G=1.+G/F
		AR(I,L)=G*AR(I,L)
		AI(I,L)=G*AI(I,L)
		IF(L.EQ.1)GO TO 270
		GO TO 170
  160		CONTINUE
		TAU(1,L)=-TAU(1,I)
		SI=TAU(2,I)
		AR(I,L)=G
  170		CONTINUE
		F=0.
		DO 240 J=1,L
			G=0.
			GI=0.
C     FORM ELEMENT A*U
			DO 180 K=1,J
				G=G+AR(J,K)*AR(I,K)+AI(J,K)*AI(I,K)
				GI=GI-AR(J,K)*AI(I,K)+AI(J,K)*AR(I,K)
  180			CONTINUE
			JP1=J+1
			IF(L.LT.JP1)GO TO 220
			DO 200 K=JP1,L
				G=G+AR(K,J)*AR(I,K)-AI(K,J)*AI(I,K)
				GI=GI-AR(K,J)*AI(I,K)-AI(K,J)*AR(I,K)
  200		CONTINUE
C     FORM ELEMENT OF P
  220		E(J)=G/H
		TAU(2,J)=GI/H
		F=F+E(J)*AR(I,J)-TAU(2,J)*AI(I,J)
  240	CONTINUE
      HH=F/(H+H)
C     FORM REDUCED A
      DO 260 J=1,L
		F=AR(I,J)
		G=E(J)-HH*F
		E(J)=G
		FI=-AI(I,J)
		GI=TAU(2,J)-HH*FI
		TAU(2,J)=-GI
		DO 260 K=1,J
			AR(J,K)=AR(J,K)-F*E(K)-G*AR(I,K)+FI*TAU(2,K)+GI*AI(I,K)
			AI(J,K)=AI(J,K)-F*TAU(2,K)-G*AI(I,K)-FI*E(K)-GI*AR(I,K)
  260 CONTINUE
  270 CONTINUE
      DO 280 K=1,L
		AR(I,K)=SCALE*AR(I,K)
		AI(I,K)=SCALE*AI(I,K)
  280 CONTINUE
      TAU(2,L)=-SI
  290 CONTINUE
      HH=D(I)
      D(I)=AR(I,I)
      AR(I,I)=HH
      AI(I,I)=SCALE*SCALE*H
  300 CONTINUE
      RETURN
      END
      
	  SUBROUTINE TQL(NM,N,D,E,Z,IER,IFVEC)
      IMPLICIT REAL*8 (A-H,O-Z)      
      DIMENSION D(N),E(N),Z(NM,N)
      EPS=1.e-25
      IER=0
      IF(N.EQ.1)GO TO 1001
      DO 100 I=2,N
  100		E(I-1)=E(I)
      F=0.
      B=0.
      E(N)=0.
      DO 240 L=1,N
		J=0
		H=EPS*( ABS(D(L))+ ABS(E(L))) 
		IF(B.LT.H)B=H 
C     LOOK FOR SMALL SUBDIAG ELT
		DO 110 M=L,N
			IF( ABS(E(M)).LE.B)GO TO 120
  110		CONTINUE
  120		IF(M.EQ.L)GO TO 220 
  130		IF(J.EQ.30)GO TO 1000 
		J=J+1 
C     FORM SHIFT
		L1=L+1
		G=D(L)
		P=(D(L1)-G)/(2.*E(L)) 
		R= SQRT(P*P+1.) 
		D(L)=E(L)/(P+ SIGN(R,P))
		H=G-D(L)
		DO 140 I=L1,N 
  140			D(I)=D(I)-H 
		F=F+H 
C     QL TRANSFORMATION 
		P=D(M)
		C=1.
		S=0.
		MML=M-L 
		DO 200 II=1,MML 
			I=M-II
			G=C*E(I)
			H=C*P 
			IF( ABS(P).LT. ABS(E(I)))GO TO 150
			C=E(I)/P
			R= SQRT(C*C+1.) 
			E(I+1)=S*P*R
			S=C/R 
			C=1./R
			GO TO 160 
  150			C=P/E(I)
			R= SQRT(C*C+1.) 
			E(I+1)=S*E(I)*R 
			S=1./R
			C=C*S 
  160			P=C*D(I)-S*G
			D(I+1)=H+S*(C*G+S*D(I)) 
			IF(IFVEC.EQ.0)GO TO 200 
C     FORM VECTOR 
		DO 180 K=1,N
				H=Z(K,I+1)
				Z(K,I+1)=S*Z(K,I)+C*H 
				Z(K,I)=C*Z(K,I)-S*H 
  180			CONTINUE
  200 CONTINUE
      E(L)=S*P
      D(L)=C*P
      IF( ABS(E(L)).GT.B)GO TO 130
  220 D(L)=D(L)+F 
  240 CONTINUE
C     ORDER VALUES AND VECTORS
      DO 300 II=2,N 
		I=II-1
		K=I 
		P=D(I)
		DO 260 J=II,N 
			IF(D(J).GE.P)GO TO 260
			K=J 
			P=D(J)
  260		CONTINUE
		IF(K.EQ.I)GO TO 300 
		D(K)=D(I) 
		D(I)=P
		IF(IFVEC.EQ.0)GO TO 300
		DO 280 J=1,N
			P=Z(J,I)
			Z(J,I)=Z(J,K) 
			Z(J,K)=P
  280		CONTINUE
  300 CONTINUE
      GO TO 1001
C     SET ERROR, NO CONVERGENCE AFTER 30 ITERATIONS 
 1000 IER=L 
 1001 RETURN
      END 

      SUBROUTINE HTRIBK(NM,N,AR,AI,TAU,M,ZR,ZI) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AR(NM,N),AI(NM,N),TAU(2,N),ZR(NM,N),ZI(NM,N)
      DO 50 K=1,N 
		DO 50 J=1,M 
			ZI(K,J)=-ZR(K,J)*TAU(2,K) 
			ZR(K,J)=ZR(K,J)*TAU(1,K)
   50 CONTINUE
      IF(N.EQ.1)GO TO 200 
      DO 140 I=2,N
		L=I-1 
		H=AI(I,I) 
		IF(H.EQ.0.)GO TO 140
		DO 130 J=1,M
			S=0.
			SI=0.
			DO 110 K=1,L
				S=S+AR(I,K)*ZR(K,J)-AI(I,K)*ZI(K,J)
				SI=SI+AR(I,K)*ZI(K,J)+AI(I,K)*ZR(K,J)
  110			CONTINUE
			S=S/H
			SI=SI/H
			DO 120 K=1,L
				ZR(K,J)=ZR(K,J)-S*AR(I,K)-SI*AI(I,K)
				ZI(K,J)=ZI(K,J)-SI*AR(I,K)+S*AI(I,K)
  120			CONTINUE
  130		CONTINUE
  140 CONTINUE
  200 RETURN
      END
C
      DOUBLE PRECISION FUNCTION PYTHAG(A,B)
      DOUBLE PRECISION A,B
C
C     FINDS DSQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW

      DOUBLE PRECISION P,R,S,T,U
      P = DMAX1(DABS(A),DABS(B))
      IF (P .EQ. 0.0D0) GO TO 20
      R = (DMIN1(DABS(A),DABS(B))/P)**2
   10 CONTINUE
	 T = 4.0D0 + R
	  IF (T .EQ. 4.0D0) GO TO 20
       S = R/T
	 U = 1.0D0 + 2.0D0*S
	 P = U*P
	 R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
C
      subroutine zero(a,n)
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension a(n)
      do 1 i=1, n
    1 a(i)=0.0
      return
	end

	subroutine IntZero(K,n)
       IMPLICIT REAL*8 (A-H,O-Z)
       dimension K(n)
       do 1 i=1, n
   1   K(i)=0
       return
	   end
	  
	  subroutine minimize(Param,y,P_lambda,NParams,iter)
	  implicit real*8 (a-h,o-z)
	  DIMENSION Param(NParams),SeedP_Lambda(NParams),P_Lambda(NParams)
	  DIMENSION PSavePlus(NParams),PSaveMinus(NParams),PNew(NParams),FixParam(NParams)
	  real*8 ysav
	  integer iter, N_seeds !, SearchSteps
	  logical NewMin
	  !Fix Param(i) is a flag to indicate that the parameter is to be fixed at it's seed value it it =1, not fixed if 0
	  N_Seeds=0
	  iter=0
	  DO 1 i=1,NParams
		FixParam(i)=0
		SeedP_Lambda(i)=P_Lambda(i)
		Pnew(i)=Param(i)
	    PsavePlus(i)=Param(i)
1		PsaveMinus(i)=Param(i)


	  ysav=funk(Pnew)
20	  iter=iter+1
25	  NewMin=.False.

	  DO 2 i=1,NParams
		IF (P_Lambda(i).eq.0) goto 11 !if parameter held constant, don't waste time
		
		PSavePlus(i)=Param(i)+P_lambda(i)
!		IF (i.eq.3) THEN
!			IF (PSavePlus(i).gt.1) PSavePlus(i)=1	  !Check for eta upper range limit
!	    ENDIF
	    y_plus=funk(PSavePlus)
	    PSaveMinus(i)=Param(i)-P_lambda(i)
!		IF (i.eq.3) THEN
!			IF (PSaveMinus(i).lt.0) PSaveMinus(i)=0	!Check for eta lower range limt
!	    ENDIF
	    y_minus=funk(PSaveMinus)
C		If boundary is local min, then reset to endpoint else use linear approx to move toward min
          IF (y_minus.lt.ysav) THEN
			ysav=y_minus
	        PNew(i)=PSaveMinus(i)
			NewMin=.True.
		ENDIF
		IF (y_plus.lt.ysav) THEN
			ysav=y_plus
			PNew(i)=PSavePlus(i)
	        NewMin=.True.
		ENDIF
		PSavePlus(i)=PSavePlus(i)-P_Lambda(i)
	    PSaveMinus(i)=PSaveMinus(i)+P_Lambda(i)
11	    CONTINUE
2	  CONTINUE !look one step in all directions around current point  
	  DO 4 i=1,NParams
4		Param(i)=PNew(i)
	    IF (NewMin) THEN 
			write(*,*) N_Seeds,iter,ysav
			GOTO 25 
	    ENDIF
	    IF (iter.lt.4) THEN
			DO 5 i=1,NParams 
5				P_Lambda(i)=P_Lambda(i)*.5	
	 		GOTO 20
	    ENDIF	
	    IF (N_Seeds.lt.10) THEN
			DO 6 i=1,NParams
6				P_Lambda(i)=SeedP_Lambda(i)
			iter=0
			N_Seeds=N_Seeds+1
			GOTO 20
	    ENDIF

!	  do 40 j=1,SearchSteps
!		Psave(1)=Pnew(1)+1.0*(j-SearchSteps/2.0)*P_lambda(1)
!		do 35 k=1,SearchSteps
!			Psave(2)=Pnew(2)+1.0*(k-SearchSteps/2.0)*P_lambda(2)
!			do 30 l=1,SearchSteps
!				Psave(3)=Pnew(3)+1.0*(l-SearchSteps/2.0)*P_lambda(3)
!				y=funk(Psave)
!				if (y.lt.ysav) then
!					ysav=y
!					Param(1)=Psave(1)
!					Param(2)=Psave(2)
!					Param(3)=Psave(3)
!!					write(*,*) 'new min=',ysav, Param1, Param2
!					NoNewMin=.False.
!				endif
! 30			continue
! 35	    continue
! 40	  continue

!	  Psave(1)=Param(1)
 !       P_lambda(1)=P_lambda(1)*.4
!	  Psave(2)=Param(2)
 !       P_lambda(2)=P_lambda(2)*.4
!	  Psave(3)=Param(3)
!	  P_lambda(3)=P_lambda(3)*.4

!	  if (iter.eq.20) goto 60
! 	  if ((ymin-ysav).gt.1e-7) goto 20 !Measure of iteration's improvement
!	  goto 20

 60	  y=ysav

	  return
	  end
