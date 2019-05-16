	PROGRAM N2NARROW_PROCESS

	IMPLICIT NONE

	REAL XSIZE,YSIZE,H,HH !system size
        REAL HLAYER,HRDF,HDPROF  !profile step
	REAL*8 RDF(500,10),RDF_AV(500,10),RDFCUT
	REAL*8 RDFSF(500),RDFSF_AV(500)
	REAL*8 DENPROF(500),DENPROF_AV(100)
        REAL*8 SLAYER(10),SLAYER_AV(10)

	REAL X(50000),Y(50000),Z(50000)

	INTEGER LAYER(50000)

	LOGICAL IFSTRFAC
	INTEGER NEIB(50000,8),NNEIB
	INTEGER NEIBNUM(50000)
	REAL PSISUMR,PSISUMI,PSISUM1(50000),PSISUM2(50000),PHI
	REAL GNR(500,3),GNX(500,3),GNR_AV(500,3),GNX_AV(500,3)
	REAL RCUTNEIB 
        REAL DENS,AV_DENS

	INTEGER*4 I,J,K,L,N,M,INDEX_RDF,NSOLID,NLAYER

	REAL XX,YY,ZZ,XX1,YY1
	REAL PI
	CHARACTER*4 6, WORD  
        CHARACTER*20 FLNM, FLNM1  

	DATA PI/3.141592652/

        DO I=1,100 !MOMENTARY DENSITY PROFILE TO ZERO
	   DENPROF_AV(I)=0.
	 END DO
        DO I=1,500 !MOMENTARY DENSITY PROFILE TO ZERO
	   RDFSF_AV(I)=0.
	 END DO
        DO K=1,10
          SLAYER(K)=0.
	  DO J=1,500
	    RDF_AV(J,K)=0.
	   END DO
	 END DO

	OPEN(FILE='general.dat',UNIT=10) 
        READ(10,*)FLNM, FLNM1 
	READ(10,*)XSIZE,YSIZE,H   !CELL SIZE X AND Y DIMENSIONS,PORE WIDTH
	READ(10,*)NLAYER ! NUMBER OF LAYERS, MAXIMUM 4
	READ(10,*)RDFCUT
	READ(10,*)IFSTRFAC,NNEIB,RCUTNEIB
	HH=H/2.
	HLAYER=H/NLAYER
        HDPROF=H/100.
	CLOSE(10)
        HRDF=RDFCUT/500.

        PSISUMR=0.0
        PSISUMI=0.0

        !write(*,*)hrdf,hlayer,hdprof
	OPEN(FILE=FLNM,UNIT=10)
	DO I=1,3000  ! 00  !ALL CONFIGURATIONS
          READ(10,*)WORD
	  IF(WORD(1:3).EQ.'END')THEN  !TRAJECTORY IS OVER - EXIT
              CLOSE(10)
              GO TO 10001
           END IF 
          READ(10,*)
          READ(10,*)
          READ(10,*)N
          READ(10,*)
	  READ(10,*)XSIZE
          XSIZE=ABS(2*XSIZE)
c          YSIZE=XSIZE
          READ(10,*)YSIZE
          YSIZE=ABS(2*YSIZE)
          READ(10,*)
          READ(10,*)
	  M=N
	  DO J=1,N  !READ COORDINATES IN THE CURRENT CONFIGURATION
	    READ(10,*)L,K,X(J),Y(J),Z(J)
c	write(*,*)' conf ',i,' mol ',j,x(j),y(j),z(j)
	    IF (ABS(Z(J)).GE.HH) THEN
	      X(J)=X(M)
	      Y(J)=Y(M)
	      Z(J)=Z(M)
	      X(M)=0.
	      Y(M)=0.
	      Z(M)=0.
	      M=M-1
	     END IF
	   END DO 
           write(*,'(A,i5,A,2F8.3,A,i7)')
     * ' FRAME:',i,' size=',xsize,ysize,' N=',N 
          DENS=N/(XSIZE*YSIZE)        
	  IF (IFSTRFAC) THEN
	    DO J=1,N
	      NEIBNUM(J)=0
	      PSISUM1(J)=0.
	      PSISUM2(J)=0.
	     END DO
	    DO J=1,500
	      DO K=1,3
	        GNR(J,K)=0.
	        GNX(J,K)=0.
	       END DO
	     END DO
	   END IF
	  N=M
	  M=0
	  DO J=1,100 !MOMENTARY DENSITY PROFILE TO ZERO
	    DENPROF(J)=0.
	   END DO
c	write(*,*)'recordx configuration ',i
	  DO J=1,N
	    LAYER(J)=INT((Z(J)+HH)/HLAYER)+1  !LAYER NUMBER
            SLAYER(LAYER(J))=SLAYER(LAYER(J))+1
	    L=INT((Z(J)+HH)/HDPROF)+1      !LAYER IN THE DENSITY PROFILE
	    DENPROF(L)=DENPROF(L)+1           !+TO DENSITY PROFILE
c	if (layer(j).gt.1)write(*,*)'!!! ',j,z(j),layer(j)
c	write(*,*)'layer',layer(j)
	   END DO
c	write(*,*)'recordu configuration ',i
          AV_DENS=(AV_DENS*(I-1)+DENS)/I
          DO J=1,10
            SLAYER(J)=SLAYER(J)/(XSIZE*YSIZE)
            SLAYER_AV(J)=(SLAYER_AV(J)*(I-1)+SLAYER(J))/I
           END DO
          DO J=1,100
             DENPROF_AV(J)=(DENPROF_AV(J)*(I-1)+DENPROF(J))/I
           END DO
                     !SOLID-FLUID DRF
	  DO J=1,500
	    DO K=1,10
	      RDF(J,K)=0.
	     END DO
	   END DO

	  DO J=1,N-1    !LOOP BY PAIRS FOR FLUID-FLUID RDFS
	    DO K= J+1,N
	      L=LAYER(J)
	      M=LAYER(K)
	      IF(L.GT.M) THEN
	        M=LAYER(J)
	        L=LAYER(K)
	       END IF
               if (max(l,m).gt.1) write(*,*)'layers',m,l 
	      INDEX_RDF=NLAYER*(L-1)+M   !NUMBER OF THE RDF
	      XX=(X(J)-X(K))
	      IF(XX.GT.XSIZE*0.5)XX=XX-XSIZE
              IF(XX.LT.-XSIZE*0.5)XX=XX+XSIZE
              YY=(Y(J)-Y(K))
              IF(YY.GT.YSIZE*0.5)YY=YY-YSIZE
              IF(YY.LT.-YSIZE*0.5)YY=YY+YSIZE
              ZZ=SQRT(XX*XX+YY*YY)
	      IF (IFSTRFAC) THEN 
c              write(*,*)zz,rcutneib,xx,x(j),x(k),xsize 
	        IF ((ZZ.LT.RCUTNEIB).AND.(L.EQ.M)) THEN !THEY ARE NEARES NEIGHBOURS
C	write(*,*)'neib',zz,l,m
	          NEIBNUM(J)=NEIBNUM(J)+1
	          NEIBNUM(K)=NEIBNUM(K)+1
	          NEIB(J,NEIBNUM(J))=K
	          NEIB(K,NEIBNUM(K))=J
	          XX=XX/ZZ
	          YY=YY/ZZ
!        if (abs(zz).le.1e-3)write(*,*)'zz?',zz 
!	if(abs(xx*xx+yy*yy-1).gt.1e-4)write(*,*)' !! ',j,k
	          PHI=ACOS(XX)  !ANGLE BETWEEN KJ VECTOR AND X AXIS
	          IF (YY.LE.0)  PHI=-PHI   
	          PSISUM1(J)=PSISUM1(J)+COS(PHI*NNEIB) !REAL PART
	          PSISUM2(J)=PSISUM2(J)+SIN(PHI*NNEIB) !COMPLEX PART
	          IF (PHI.GE.0) THEN !ANGLE BETWEEN KJ VECTOR AND X AXIS
	            PHI=PHI-PI
	           ELSE
	            PHI=PHI+PI
	           END IF
	          PSISUM1(K)=PSISUM1(K)+COS(PHI*NNEIB) !REAL PART
	          PSISUM2(K)=PSISUM2(K)+SIN(PHI*NNEIB) !COMPLEX PART
	         END IF
	       END IF
	      IF (ZZ.LT.RDFCUT) THEN
	        L=INT(ZZ/HRDF)+1
	        IF((L.GT.500).OR.(INDEX_RDF.GT.9)) THEN
	          WRITE(*,*)I,J,K,L,INDEX_RDF
	         ELSE
	          RDF(L,INDEX_RDF)=RDF(L,INDEX_RDF)+1
	         END IF
	       END IF
	     END DO
	   END DO 
!        write(*,*)'strfac'
	  IF (IFSTRFAC) THEN
	    DO J=1,N  !FINAL PSI SUMMS 
              IF (NEIBNUM(J).GT.0)THEN
	        PSISUM1(J)=PSISUM1(J)/NEIBNUM(J) !REAL PART
	        PSISUM2(J)=PSISUM2(J)/NEIBNUM(J) !COMPLEX PART
               ELSE
                PSISUM1(J)=0. !REAL PART
                PSISUM2(J)=0. !COMPLEX PART
               ENDIF
!            if (neibnum(j).eq.0)write(*,*)neibnum(j)
!	write(*,*)' ** ', psisum1(j),psisum2(j),
!     *    sqrt(psisum1(j)**2+psisum2(j)**2),neibnum(j)
              PSISUMR=PSISUMR+sqrt(psisum1(j)**2+psisum2(j)**2)
	     END DO
            PSISUMI=PSISUMI+N
	    DO J=1,N-1    !LOOP BY PAIRS FOR STRUCTURE FUNCTION
	      DO K= J+1,N
	        IF (LAYER(J).EQ.LAYER(K)) THEN
	          M=LAYER(J)
	          XX=ABS(X(J)-X(K))
	          IF(XX.GT.XSIZE*0.5)XX=XSIZE-XX
	          YY=ABS(Y(J)-Y(K))
	          IF(YY.GT.YSIZE*0.5)YY=YSIZE-YY
	          ZZ=SQRT(XX*XX+YY*YY)
	          IF (ZZ.LE.RDFCUT) THEN
	            INDEX_RDF=INT(ZZ/HRDF)+1
	            ZZ=PSISUM1(J)*PSISUM1(K)+PSISUM2(J)*PSISUM2(K)
	            GNR(INDEX_RDF,M)=GNR(INDEX_RDF,M)+ZZ
	            GNX(INDEX_RDF,M)=GNX(INDEX_RDF,M)+1
	           END IF
	         END IF
	       END DO
	     END DO	          
	    DO J=1,500
!        write(*,*)GNR(J,1),GNX(J,1)
	      DO K=1,3
	        GNX_AV(J,K)=(GNX_AV(J,K)*(I-1)+GNX(J,K))/I
	        IF (GNX(J,K).GT.0) THEN
	           GNR_AV(J,K)=(GNR_AV(J,K)*(I-1)+GNR(J,K)/GNX(J,K))/I
	         END IF
	       END DO
	     END DO
	   END IF 
!            DO J=1,500
!        write(*,*)GNR_AV(J,1),GNX_AV(J,1)
!       end do 
	  DO J=1,10
	    DO K=1,500
	      RDF_AV(K,J)=(RDF_AV(K,J)*(I-1)+RDF(K,J))/I
	     END DO
	   END DO
c	write(*,*)'end configuration ',i
c	do k=1,n
c	write(*,*)k,neibnum(k)
c	end do
	 END DO

	CLOSE(10)

10001	OPEN(FILE=FLNM1,UNIT=10)	
	WRITE(10,50001)N,I,NLAYER,HLAYER,HRDF
	WRITE(10,*)'DENSITY PROFILE'
        XX=XSIZE*YSIZE*HDPROF
	DO I=1,100
	  WRITE(10,50002)HDPROF*(I-0.5)-HH,DENPROF_AV(I)/XX
	 END DO
50003	FORMAT(F8.4,10E13.6)
50001 FORMAT
     *(' N=',I4,' I=',I8,' NLAY=',I2,' hlayer=',f7.4,' HRFD=',F7.4)
        WRITE(10,*)'FLUID-FLUID RADIAL FUNCTIONS'
	DO I=1,500
          DO J=1,10
            L=((J-1)/NLAYER)+1
            M=J-(L-1)*NLAYER
            IF(L.GT.M)THEN
              K=M
              M=L
              L=K
             END IF
            IF(L.EQ.M)THEN
              XX=(SLAYER_AV(L)*XSIZE*YSIZE-1)*SLAYER_AV(L)*0.5
             ELSE
              XX=SLAYER_AV(M)*SLAYER_AV(L)*XSIZE*YSIZE
             END IF
C             WRITE(*,*)I,J,L,M,XX
             IF ((L.LE.NLAYER).AND.(M.LE.NLAYER))THEN
               RDF_AV(I,J)=RDF_AV(I,J)/(PI*(2*I-1)*HRDF*HRDF*XX)
              ELSE
               RDF_AV(I,J)=0.
              END IF
           END DO
	  WRITE(10,50003)HRDF*(I-0.5),(RDF_AV(I,J),J=1,10)
	 END DO
50002	FORMAT(F8.4,E13.6)
!        do j=1,500
!        write(*,*)'******',GNR_AV(j,1),gnx_av(j,1)
!        end do
	WRITE(10,*)' STRUCTURE FACTOR '
	DO I=1,500
          XX=(2*I-1)*HRDF*HRDF
	  WRITE(10,50004)HRDF*(I-0.5),GNR_AV(I,1),
     *      GNR_AV(I,2),GNR_AV(I,3),GNX_AV(I,1)
	 END DO
        WRITE(10,*)' HEX ORDER PARAMETER',PSISUMR/PSISUMI,
     &    ' DENSITY',AV_DENS
50004	FORMAT(F8.4,4E14.6)
	CLOSE(10)

	STOP
	END


   


	


