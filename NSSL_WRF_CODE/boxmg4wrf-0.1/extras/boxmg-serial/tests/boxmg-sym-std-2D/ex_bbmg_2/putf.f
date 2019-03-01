      SUBROUTINE PUTF (SO, QF, L, M, HX, HY, K)
C
C***PURPOSE  THIS PUTF ROUTINE DEFINES THE COEFFICIENT MATRIX AND THE
C            RIGHT HAND SIDE OF THE SYSTEM OF EQUATIONS FOR THE
C            STANDARD FIVE-POINT DIFFERENCE APPROXIMATION TO THE
C            DIVERGENCE EQUATION,
C               -(D/DX)[(X+Y)(DU/DX)] - (D/DY)[(X-Y)(DU/DY)]
C                           + X*Y U(X,Y) = X - Y + (X**2)*(Y**2)
C
C   THE EQUATION IS ASSUMED TO HOLD ON THE RECTANGLE
C
C         { (X,Y):  3 .LE. X .LE. 4  ,  1 .LE. Y .LE.  2 }.
C
C   THE BOUNDARY CONDITIONS ARE:
C
C       1.  U(3,Y) = 3*Y,
C
C       2.  (DU/DN)(4,Y) = Y,
C
C       3.  (DU/DN)(X,1) - X U(X,1) = X - X**2,
C
C       4.  U(X,2) = 2*X
C
C   THE APPROXIMATION IS ACCOMPLISHED BY MAKING SECOND-ORDER
C   APPROXIMATIONS TO ALL FIRST-DERIVATIVES ON A UNIFORM FINITE
C   DIFFERENCE GRID WITH GRID SPACINGS DX AND DY. ALL DIFFERENCE
C   EQUATIONS ARE SCALED BY DX*DY.
C
C   THE GENERAL EQUATION AT AN INTERIOR GRIDPOINT IS :
C
C        -(DY/DX){ P(X(I+1/2),Y(J))[U(X(I+1),Y(J))-U(X(I),Y(J))]
C                -P(X(I-1/2),Y(J))[U(X(I),Y(J))-U(X(I-1),Y(J))] }
C
C        -(DX/DY){ Q(X(I),Y(J+1/2))[U(X(I),Y(J+1))-U(X(I),Y(J))]
C                -Q(X(I),Y(J-1/2))[U(X(I),Y(J))-U(X(I),Y(J-1))] }
C
C       + DX*DY*R(X(I),Y(J))U(X(I),Y(J)) = DX*DY*F(X(I),Y(J)),
C
C   WHERE THE NOTATION X(I+1/2) MEANS THE POINT MIDWAY BETWEEN THE TWO
C   GRIDPOINTS X(I) AND X(I+1) WITH SIMILAR DEFINITIONS FOR LIKE TERMS.
C
C   AT THE BOUNDARIES OF THE RECTANGLE NORMAL DERIVATIVES ARE
C   APPROXIMATED BY A SECOND-ORDER CENTRAL DIFFERENCE OVER A SPACING OF
C   2*DX OR 2*DY.  THIS EQUATION IS USED TO ELIMINATE THE UNKNOWN AT A
C   VIRTUAL GRIDPOINT THAT LIES JUST OUTSIDE THE RECTANGLE.
C
C    THE GENERAL EQUATION AT GRIDPOINTS (2,J) IS :
C
C        -(DY/DX){ P(X(5/2),Y(J))[U(X(3),Y(J))-U(X(2),Y(J))]
C                -P(X(3/2),Y(J))[U(X(2),Y(J))]}
C
C        -(DX/DY){ Q(X(2),Y(J+1/2))[U(X(2),Y(J+1))-U(X(2),Y(J))]
C                -Q(X(2),Y(J-1/2))[U(X(2),Y(J))-U(X(2),Y(J-1))] }
C
C       + DX*DY*R(X(2),Y(J))U(X(2),Y(J)) = DX*DY*F(X(2),Y(J)) +
C                                        (DY/DX)*(P(X(3/2),Y(J))*U(XS,Y)
C
C    THE GENERAL EQUATION AT GRIDPOINTS ON SOUTH BOUNDARY IS :
C
C        -(DY/DX){ P(X(I+1/2),Y(1))[U(X(I+1),Y(1))-U(X(I),Y(1))]
C                -P(X(I-1/2),Y(1))[U(X(I),Y(1))-U(X(I-1),Y(1))] }
C
C        -(DX/DY){ Q(X(I),Y(3/2))[U(X(I),Y(2))-U(X(I),Y(1))]
C                -Q(X(I),Y(1/2))[U(X(I),Y(1)) - U(X(I),Y(2)) + 
C				2.*B*DY*U(X(I),Y(1)/ALPHA(X(I),Y(1))]}
C
C       + DX*DY*R(X(I),Y(1))U(X(I),Y(1)) = DX*DY*F(X(I),Y(1)) -
C                               2.*DY*U(X(I),Y(YS))/ALPHA(X(I),Y(1))
C                          
C
C    THE GENERAL EQUATION AT GRIDPOINTS ON EAST BOUNDARY IS :
C
C        -(DY/DX){ P(X(L+1/2),Y(J))[U(X(L+1),Y(J))-U(X(L),Y(J))]
C                -P(X(I-1/2),Y(J))[U(X(L),Y(J))-U(X(L-1),Y(J))] }
C
C        -(DX/DY){ Q(X(L),Y(J+1/2))[U(X(L),Y(J+1))-U(X(L),Y(J))]
C                -Q(X(L),Y(J-1/2))[U(X(L),Y(J))-U(X(L),Y(J-1))] }
C
C         -(DY/DX){P(X(L+1/2),Y(J))[U(X(L-1),Y(J)]
C
C       + DX*DY*R(X(L),Y(J))U(X(L),Y(J)) = DX*DY*F(X(L),Y(J)) +
C                     2.(DY/DX)*P(X(L+1/2),Y(J))*G(X(XF),Y(J))
C
C    THE GENERAL EQUATION AT GRIDPOINTS (I,M-1) IS :
C
C        -(DY/DX){ P(X(I+1/2),Y(M-1))[U(X(I+1),Y(M-1))-U(X(I),Y(M-1))]
C                -P(X(I-1/2),Y(M-1))[U(X(I),Y(M-1))-U(X(I-1),Y(M-1))] }
C
C        -(DX/DY){ Q(X(I),Y(M-1/2))[-U(X(I),Y(M-1))]
C                -Q(X(I),Y(M-3/2))[U(X(I),Y(M-1))-U(X(I),Y(M-2))] }
C
C       + DX*DY*R(X(I),Y(M-1))U(X(I),Y(M-1)) = DX*DY*F(X(I),Y(M-1)) +
C                              (DX/DY){Q(X(I),Y(M-1/2)[U(X(I),Y(YF))]}
C
C
C***INPUTS
C
C   L		THE NUMBER OF UNKNOWNS IN THE X-DIMENSION
C   M		THE NUMBER OF UNKNOWNS IN THE Y-DIMENSION
C   HX		THE X-DIRECTION GRID SPACING.
C   HY		THE Y-DIRECTION GRID SPACING.
C   K		THE GRID LEVEL FOR THE DIFFERENCE OPERATOR AND RIGHT HAND SIDE
C		TO BE COMPUTED ON. K = 0, IS THE FINEST GRID.
C
C***INPUT/OUTPUT
C
C   SO (I,J,KK)	REAL ARRAY CONTAINING THE COEFFICIENTS OF THE MATRIX, OBTAINED
C		FROM THE DIFFENCE OPERATOR. THE INDEXES ARE DEFINED AS THE
C		STENCIL CENTERED AT GRID POINT (I,J) 5 POINT STENCIL POINT KK.
C               SEE PARAMETERS FOR DEFINITION OF KK.
C   QF (I,J)	REAL ARRAY CONTAINING THE RIGHT HAND SIDE FOR THE EQUATION
C		ASSOCIATED WITH THE  5 POINT STENCIL CENTERED AT GRID 
C               POINT (I,J).
C
C
C---PARAMETERS
C
C   KO = 1	INDEX FOR THE CENTER POINT COEFFICIENT
C   KW = 2	INDEX FOR THE WEST POINT COEFFICIENT
C   KS = 3	INDEX FOR THE SOUTH POINT COEFFICIENT
C   KNW = 4	INDEX FOR THE NORTH-WEST POINT COEFFICIENT
C   KSW = 5	INDEX FOR THE SOUTH-WEST POINT COEFFIcient
C
C

      IMPLICIT NONE

      REAL*8   XS, XF, YS, YF
      COMMON   /BBMGDM/ XS, XF, YS, YF
C
C
      INTEGER  KO, KW, KS
      PARAMETER ( KO = 1, KW = 2, KS = 3)
      INTEGER  P, Q, R, F, GW, GE, BS, GS, GN
      PARAMETER( P  = 1,
     &           Q  = 2,
     &           R  = 3,
     &           F  = 4,
     &           GW = 6,
     &           GE = 8,
     &           BS = 9,
     &           GS = 10,
     &           GN = 12)
C
      INTEGER   K, L, M, I, J
      REAL*8    SO(0:L+1,0:M+1,3), PCOEF, QF(0:L+1,0:M+1),
     &          HX, HY, H2, P5, YXH, XYH,
     &          X, Y
C
C***FIRST EXECUTABLE STATEMENT
C
C-- INITIALIZE RIGHT HAND SIDE AND DIFFERENCE OPERATOR ON GRID K.
C
      P5=0.5D0
      H2 = HX*HY
      YXH = HY/HX
      XYH = HX/HY
C
C-- COMPUTE ALL SOUTH, WEST, AND CENTER STENCIL COEFFICIENTS, AND
C   COMPUTE THE RIGHT HAND SIDE FOR EACH EQUATION.
C   NOTE: THE SOUTH COEFFICIENTS FOR THE POINTS (I,1) I=1,l ARE SET To
C        ZERO, AND THE COEFFICIENTS FOR THE POINTS (1,J) J=1,M ARE SET TO
C        ZERO BECAUSE THEY REFER TO FICTICIOUS POINTS.
C

      DO 10 J = 2, M
        Y = YS + (J-P5)*HY - HY
        DO 10 I = 1, L
          X = XS + I*HX 
          SO(I,J,KS) = XYH*PCOEF(Q,X,Y)
   10   CONTINUE

      DO 20 J = 1, M
        Y = YS + J*HY - HY
        DO 20 I = 2, L
          X = XS + (I-P5)*HX 
          SO(I,J,KW) = YXH*PCOEF(P,X,Y)
   20   CONTINUE

      DO 30 J = 2, M
        Y = YS + J*HY - HY
        DO 30 I = 1, L-1
          X = XS + I*HX 
          SO(I,J,KO) = XYH*(PCOEF(Q,X,Y-P5*HY) +
     &                     PCOEF(Q,X,Y+P5*HY)) +
     &                 YXH*(PCOEF(P,X-P5*HX,Y) +
     &                     PCOEF(P,X+P5*HX,Y)) +
     &                 H2*PCOEF(R,X,Y)
          QF(I,J) = H2*PCOEF(F,X,Y)
   30 CONTINUE
C
C--- INITIALIZE EDGE (BOUNDARY) POINTS
C
C--- INITIALIZE SOUTH BOUNDARY ( YS EDGE OF GRID )
C
        Y = YS + P5*HY
        DO 32 I = 1, L-1
          X = XS + I*HX 
          SO(I,1,KO) = 2.*XYH*PCOEF(Q,X,Y) +
     &                 YXH*(PCOEF(P,X-.5*HX,YS) + 
     &                     PCOEF(P,X+.5*HX,YS)) + 
     &                 H2*PCOEF(R,X,YS)
          QF(I,1) = H2*PCOEF(F,X,YS)
   32   CONTINUE
C
C--- INITIALIZE EAST BOUNDARY ( XF EDGE OF GRID )
C
        X = XS + (L-P5)*HX 
        DO 38 J = 2, M
          Y = YS + J*HY - HY
          SO(L,J,KO) = XYH*(PCOEF(Q,XF,Y-P5*HY) +
     &                     PCOEF(Q,XF,Y+P5*HY)) +
     &                 2*YXH*PCOEF(P,X,Y) +
     &                 H2*PCOEF(R,XF,Y)
          QF(L,J) = H2*PCOEF(F,XF,Y)
  38    CONTINUE
C
C--- INITIALIZE CORNER POINTS
C
C--- SE CORNER - ( YS IS ROBIN'S BC )
C
          X = XS + (L-P5)*HX 
          Y = YS + P5*HY
          SO(L,1,KO) = 2*XYH*PCOEF(Q,XF,Y) +
     &                 2*YXH*PCOEF(P,X,YS) +
     &                 H2*PCOEF(R,XF,YS)
          QF(L,1) = H2*PCOEF(F,XF,YS)
C
C--- ADJUST THE RHS AND EQUATIONS FOR THE BOUNDARY CONDITIONS
C
C--- ADJUSTMENTS FOR THE SOUTH BOUNDARY
C
        Y = YS + P5*HY
        DO 44 I= 1, L
          X = XS + I*HX 
          SO(I,1,KO) = SO(I,1,KO) - 2*HX*PCOEF(Q,X,Y)*PCOEF(BS,X,YS)
          QF(I,1) = QF(I,1) - 2*HX*PCOEF(Q,X,Y)*PCOEF(GS,X,YS)
   44   CONTINUE
C
C--- ADJUSTMENTS FOR THE NORTH BOUNDARY
C
        Y = YS + (M+P5)*HY - HY
        DO 50 I = 1, L
          X = XS + I*HX 
          QF(I,M) = QF(I,M) + XYH*PCOEF(Q,X,Y)*PCOEF(GN,X,YF)
   50   CONTINUE
C
C--- ADJUSTMENTS FOR THE WEST BOUNDARY
C
        X = XS + P5*HX 
        DO 60 J = 1, M
          Y = YS + J*HY - HY
          QF(1,J) = QF(1,J) + YXH*PCOEF(P,X,Y)*PCOEF(GW,XS,Y)
   60   CONTINUE
C
C--- ADJUSTMENTS FOR THE EAST BOUNDARY
C
        X = XS + (L-P5)*HX 
        DO 72 J = 1, M
          Y = YS + J*HY - HY
          QF(L,J) = QF(L,J) + 2*HY*PCOEF(P,X,Y)*PCOEF(GE,XF,Y)
   72   CONTINUE
C
C
C--- SCALE EQUATIONS FOR NEUMANN AND ROBINS BOUNDARY CONDITION
C
        DO 80 I = 1, L
          SO(I,1,KW) = P5*SO(I,1,KW)
          SO(I,1,KO) = P5*SO(I,1,KO)
          QF(I,1) = P5*QF(I,1) 
   80 CONTINUE
        DO 87 J = 1, M
          SO(L,J,KS) = P5*SO(L,J,KS)
          SO(L,J,KO) = P5*SO(L,J,KO)
          QF(L,J) = P5*QF(L,J) 
   87 CONTINUE
C
      RETURN
      END



      REAL*8 FUNCTION PCOEF(IFLAG,X,Y)

C***PURPOSE   THIS PCOEF ROUTINE DETERMINES THE COEFFICIENTS OF THE  
C             PARTIAL DIFFERENTIAL EQUATION AND BOUNDARY CONDITIONS.
C
C     COEFFICIENTS OF THE DIVERGENCE EQUATION, P(X,Y), Q(X,Y),
C     R(X,Y), AND F(X,Y)
C
      IMPLICIT NONE

      INTEGER  IFLAG
      REAL*8   X, Y

      IF(IFLAG.EQ.1) THEN
	 PCOEF= X + Y
      ELSE IF(IFLAG.EQ.2) THEN
         PCOEF = X - Y
      ELSE IF(IFLAG.EQ.3) THEN
          PCOEF= X * Y
      ELSE IF(IFLAG.EQ.4) THEN
          PCOEF= X - Y + (X**2)*(y**2)

C     XS BONDARY COEFFICIENT G(XS,Y)

      ELSE IF(IFLAG.EQ.6) THEN 
          PCOEF = 3 * Y

C     XF BOUNDARY COEFFICIENT G(XF,Y)

      ELSE IF(IFLAG.EQ.8) THEN 
          PCOEF = Y

C     YS BOUNDARY COEFFICIENTS B(X,YS) AND G(X,YS)

      ELSE IF(IFLAG.EQ.9) THEN
           PCOEF = - X
      ELSE IF(IFLAG.EQ.10) THEN
         PCOEF = X - X**2

C     YF BOUNDARY COEFFICIENT G(X,YF)
 
      ELSE IF(IFLAG.EQ.12) THEN
	   PCOEF = 2 * X
      ENDIF 

      RETURN
      END










