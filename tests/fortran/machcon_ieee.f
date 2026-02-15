C     IEEE 754 double precision machine constants for TOMS 644.
C     Replaces machcon.f (which uses DATA statements for specific machines).
C     Also provides XERROR/FDUMP stubs needed by TOMS 644 error handling.
C
C     This file is part of the complex-bessel test infrastructure.
C     It is NOT a reference file — it is our own code for generating
C     reference values from the original Fortran TOMS 644 routines.
C
      DOUBLE PRECISION FUNCTION D1MACH(I)
      INTEGER I
C
C     IEEE 754 double precision machine constants:
C       D1MACH(1) = smallest normalized positive number
C       D1MACH(2) = largest finite number
C       D1MACH(3) = smallest relative spacing (0.5 * epsilon)
C       D1MACH(4) = largest relative spacing (epsilon)
C       D1MACH(5) = log10(radix) = log10(2)
C
      DOUBLE PRECISION DMACH(5)
      DATA DMACH(1) / 2.2250738585072014D-308 /
      DATA DMACH(2) / 1.7976931348623157D+308 /
      DATA DMACH(3) / 1.1102230246251565D-016 /
      DATA DMACH(4) / 2.2204460492503131D-016 /
      DATA DMACH(5) / 3.0102999566398120D-001 /
C
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'D1MACH: I OUT OF BOUNDS', I
         STOP
      END IF
      D1MACH = DMACH(I)
      RETURN
      END
C
      INTEGER FUNCTION I1MACH(I)
      INTEGER I
C
C     IEEE 754 / LP64 machine integer constants:
C       I1MACH( 1) = standard input unit   = 5
C       I1MACH( 2) = standard output unit  = 6
C       I1MACH( 3) = standard punch unit   = 6
C       I1MACH( 4) = standard error unit   = 0
C       I1MACH( 5) = bits per integer      = 32
C       I1MACH( 6) = chars per integer     = 4
C       I1MACH( 7) = base for integers     = 2
C       I1MACH( 8) = digits in integer base = 31
C       I1MACH( 9) = largest integer       = 2147483647
C       I1MACH(10) = base for floats       = 2
C       I1MACH(11) = digits in float significand = 24 (single)
C       I1MACH(12) = min float exponent    = -125 (single)
C       I1MACH(13) = max float exponent    = 128 (single)
C       I1MACH(14) = digits in double significand = 53
C       I1MACH(15) = min double exponent   = -1021
C       I1MACH(16) = max double exponent   = 1024
C
      INTEGER IMACH(16)
      DATA IMACH( 1) /    5 /
      DATA IMACH( 2) /    6 /
      DATA IMACH( 3) /    6 /
      DATA IMACH( 4) /    0 /
      DATA IMACH( 5) /   32 /
      DATA IMACH( 6) /    4 /
      DATA IMACH( 7) /    2 /
      DATA IMACH( 8) /   31 /
      DATA IMACH( 9) / 2147483647 /
      DATA IMACH(10) /    2 /
      DATA IMACH(11) /   24 /
      DATA IMACH(12) / -125 /
      DATA IMACH(13) /  128 /
      DATA IMACH(14) /   53 /
      DATA IMACH(15) / -1021 /
      DATA IMACH(16) /  1024 /
C
      IF (I .LT. 1 .OR. I .GT. 16) THEN
         WRITE(*,*) 'I1MACH: I OUT OF BOUNDS', I
         STOP
      END IF
      I1MACH = IMACH(I)
      RETURN
      END
C
C     Stub for XERROR — TOMS 644 calls this for error reporting.
C     We just print and return (no abort).
C
      SUBROUTINE XERMSG(LIBRAR, SUBROU, MESSG, NERR, LEVEL)
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      INTEGER NERR, LEVEL
C     Silent stub — errors are reported via IERR return code
      RETURN
      END
C
      SUBROUTINE FDUMP
C     Stub — no dump needed
      RETURN
      END
