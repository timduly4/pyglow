
      FUNCTION COSD(X)
          IMPLICIT NONE
          REAL :: COSD, UMR
          REAL, INTENT(IN) :: X
          UMR  = ATAN(1.0)*4./180. 
          COSD = COS(X*UMR)

      END FUNCTION COSD

      FUNCTION SIND(X)
          IMPLICIT NONE
          REAL :: SIND, UMR
          REAL, INTENT(IN) :: X
          UMR  = ATAN(1.0)*4./180. 
          SIND = SIN(X*UMR)

      END FUNCTION SIND

