C     Generate reference values for K Bessel function tests.
C     Calls ZBESK (upper interface) and ZBKNU (core subroutine) from
C     TOMS 644 (zbsubs.f) with IEEE 754 machine constants.
C
C     Output format: structured text parsed by txt_to_json.py
C     Each test case is delimited by BEGIN/END markers.
C
C     Usage: compile with zbsubs.f + machcon_ieee.f, then run.
C
      PROGRAM GEN_BESK_REFS
      IMPLICIT NONE
C
      DOUBLE PRECISION ZR, ZI, FNU, TOL, ELIM, ALIM
      DOUBLE PRECISION CYR(10), CYI(10)
      DOUBLE PRECISION D1MACH
      INTEGER I1MACH
      INTEGER KODE, N, NZ, IERR, I
C
C     Compute machine constants needed by ZBKNU
C     (same formulas as in ZBESK lines 57-68 of zbsubs.f)
C
      DOUBLE PRECISION R1M5, RL, DIG, FNUL, BB, AZ
      R1M5 = D1MACH(5)
      TOL = MAX(D1MACH(4), 1.0D-18)
      RL = 1.0D+3*D1MACH(4)*4.7D-1
C     ELIM = -(ln(smallest) - 3)
      ELIM = -2.302585093D0*(DBLE(I1MACH(15))*R1M5 + 3.0D0)
C     ALIM: threshold for exponential underflow range
      DIG = DMIN1(-DLOG10(TOL), DBLE(I1MACH(14))*R1M5)
      ALIM = ELIM + DLOG(DMAX1(1.0D-15,
     *       0.5D0*DSQRT(2.0D0)*10.0D0**(-DIG)))
C
      WRITE(*,'(A)') '# Reference values for complex-bessel K tests'
      WRITE(*,'(A)') '# Generator: TOMS 644 zbsubs.f (rev 930101)'
      WRITE(*,'(A)') '# Machine: IEEE 754 double precision'
      WRITE(*,'(A)') '#'
C
C     ============================================================
C     SECTION 1: ZBESK (upper interface) tests
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbesk'
C
C     --- Real axis tests ---
C
C     K_0(1), unscaled
      CALL BESK_TEST(1.0D0, 0.0D0, 0.0D0, 1, 1, 'K_0(1.0)')
C     K_1(1), unscaled
      CALL BESK_TEST(1.0D0, 0.0D0, 1.0D0, 1, 1, 'K_1(1.0)')
C     K_0(0.01), small argument (series path)
      CALL BESK_TEST(0.01D0, 0.0D0, 0.0D0, 1, 1, 'K_0(0.01)')
C     K_0(0.1), small argument
      CALL BESK_TEST(0.1D0, 0.0D0, 0.0D0, 1, 1, 'K_0(0.1)')
C     K_0(2), medium argument
      CALL BESK_TEST(2.0D0, 0.0D0, 0.0D0, 1, 1, 'K_0(2.0)')
C     K_1(2), medium argument
      CALL BESK_TEST(2.0D0, 0.0D0, 1.0D0, 1, 1, 'K_1(2.0)')
C     K_0(5), medium-large argument
      CALL BESK_TEST(5.0D0, 0.0D0, 0.0D0, 1, 1, 'K_0(5.0)')
C     K_0(10), large argument (asymptotic path)
      CALL BESK_TEST(10.0D0, 0.0D0, 0.0D0, 1, 1, 'K_0(10.0)')
C     K_0(50), very large argument
      CALL BESK_TEST(50.0D0, 0.0D0, 0.0D0, 1, 1, 'K_0(50.0)')
C     K_{0.5}(1), half-integer order
      CALL BESK_TEST(1.0D0, 0.0D0, 0.5D0, 1, 1, 'K_0.5(1.0)')
C     K_{0.25}(1), fractional order
      CALL BESK_TEST(1.0D0, 0.0D0, 0.25D0, 1, 1, 'K_0.25(1.0)')
C     K_{2.5}(3), medium fractional order
      CALL BESK_TEST(3.0D0, 0.0D0, 2.5D0, 1, 1, 'K_2.5(3.0)')
C     K_5(2), higher integer order
      CALL BESK_TEST(2.0D0, 0.0D0, 5.0D0, 1, 1, 'K_5(2.0)')
C
C     --- Complex argument tests ---
C
C     K_0(1+i)
      CALL BESK_TEST(1.0D0, 1.0D0, 0.0D0, 1, 1, 'K_0(1+i)')
C     K_1(1+i)
      CALL BESK_TEST(1.0D0, 1.0D0, 1.0D0, 1, 1, 'K_1(1+i)')
C     K_0(1.5+0.5i)
      CALL BESK_TEST(1.5D0, 0.5D0, 0.0D0, 1, 1, 'K_0(1.5+0.5i)')
C     K_{0.25}(1.5+0.5i), fractional order + complex
      CALL BESK_TEST(1.5D0, 0.5D0, 0.25D0, 1, 1, 'K_0.25(1.5+0.5i)')
C     K_0(5+3i), large complex
      CALL BESK_TEST(5.0D0, 3.0D0, 0.0D0, 1, 1, 'K_0(5+3i)')
C     K_0(0.5+2i), imaginary-dominant
      CALL BESK_TEST(0.5D0, 2.0D0, 0.0D0, 1, 1, 'K_0(0.5+2i)')
C     K_2(3+4i), medium order + large complex
      CALL BESK_TEST(3.0D0, 4.0D0, 2.0D0, 1, 1, 'K_2(3+4i)')
C
C     --- Scaled (KODE=2) tests ---
C
C     K_0(1), scaled
      CALL BESK_TEST(1.0D0, 0.0D0, 0.0D0, 2, 1, 'K_0(1.0) scaled')
C     K_0(10), scaled
      CALL BESK_TEST(10.0D0, 0.0D0, 0.0D0, 2, 1, 'K_0(10.0) scaled')
C     K_0(1+i), scaled
      CALL BESK_TEST(1.0D0, 1.0D0, 0.0D0, 2, 1, 'K_0(1+i) scaled')
C
C     --- Sequence tests (N>1) ---
C
C     K_0(2), K_1(2), K_2(2) (n=3 starting from fnu=0)
      CALL BESK_TEST(2.0D0, 0.0D0, 0.0D0, 1, 3, 'K_{0,1,2}(2.0)')
C     K_1(1+i), K_2(1+i), K_3(1+i)
      CALL BESK_TEST(1.0D0, 1.0D0, 1.0D0, 1, 3,
     *               'K_{1,2,3}(1+i)')
C     K_{0.5}(3), K_{1.5}(3), K_{2.5}(3)
      CALL BESK_TEST(3.0D0, 0.0D0, 0.5D0, 1, 3,
     *               'K_{0.5,1.5,2.5}(3.0)')
C
C     ============================================================
C     SECTION 2: ZBKNU (core subroutine) tests
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbknu'
C
C     --- Series path (|z| small) ---
C
C     zbknu K_0(0.5), series path
      CALL BKNU_TEST(0.5D0, 0.0D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(0.5) series')
C     zbknu K_0(1.0), near boundary |z|~2
      CALL BKNU_TEST(1.0D0, 0.0D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(1.0)')
C     zbknu K_1(0.5), series path, integer order
      CALL BKNU_TEST(0.5D0, 0.0D0, 1.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_1(0.5) series')
C     zbknu K_{0.25}(0.3), fractional order, small z
      CALL BKNU_TEST(0.3D0, 0.0D0, 0.25D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0.25(0.3) series')
C     zbknu K_0(0.1+0.1i), complex small z
      CALL BKNU_TEST(0.1D0, 0.1D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(0.1+0.1i) series')
C
C     --- Asymptotic/Miller path (|z| larger) ---
C
C     zbknu K_0(5), asymptotic
      CALL BKNU_TEST(5.0D0, 0.0D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(5.0) asymptotic')
C     zbknu K_0(10), large real
      CALL BKNU_TEST(10.0D0, 0.0D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(10.0) asymptotic')
C     zbknu K_{2.5}(3), fractional medium
      CALL BKNU_TEST(3.0D0, 0.0D0, 2.5D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_2.5(3.0)')
C     zbknu K_0(5+3i), complex asymptotic
      CALL BKNU_TEST(5.0D0, 3.0D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(5+3i) asymptotic')
C     zbknu K_0(1+i), complex medium
      CALL BKNU_TEST(1.0D0, 1.0D0, 0.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_0(1+i)')
C
C     --- Sequence (N>1) ---
C
C     zbknu K_{0,1}(2), integer sequence
      CALL BKNU_TEST(2.0D0, 0.0D0, 0.0D0, 1, 2,
     *               TOL, ELIM, ALIM, 'K_{0,1}(2.0)')
C     zbknu K_{0.5,1.5}(3), fractional sequence
      CALL BKNU_TEST(3.0D0, 0.0D0, 0.5D0, 1, 2,
     *               TOL, ELIM, ALIM, 'K_{0.5,1.5}(3.0)')
C
C     --- Scaled (KODE=2) ---
C
C     zbknu K_0(1) scaled
      CALL BKNU_TEST(1.0D0, 0.0D0, 0.0D0, 2, 1,
     *               TOL, ELIM, ALIM, 'K_0(1.0) scaled')
C     zbknu K_0(10) scaled
      CALL BKNU_TEST(10.0D0, 0.0D0, 0.0D0, 2, 1,
     *               TOL, ELIM, ALIM, 'K_0(10.0) scaled')
C     zbknu K_5(2), high order
      CALL BKNU_TEST(2.0D0, 0.0D0, 5.0D0, 1, 1,
     *               TOL, ELIM, ALIM, 'K_5(2.0)')
C
      STOP
      END
C
C     ============================================================
C     Helper: call ZBESK and print results
C     ============================================================
C
      SUBROUTINE BESK_TEST(ZR, ZI, FNU, KODE, N, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU
      INTEGER KODE, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION CYR(10), CYI(10)
      INTEGER NZ, IERR, J
C
      CALL ZBESK(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zbesk ', LABEL
      WRITE(*,'(A,E25.17)') 'z_re  ', ZR
      WRITE(*,'(A,E25.17)') 'z_im  ', ZI
      WRITE(*,'(A,E25.17)') 'fnu   ', FNU
      WRITE(*,'(A,I6)')     'kode  ', KODE
      WRITE(*,'(A,I6)')     'n     ', N
      WRITE(*,'(A,I6)')     'nz    ', NZ
      WRITE(*,'(A,I6)')     'ierr  ', IERR
      DO J = 1, N
         WRITE(*,'(A,I2,A,E25.17)') 'cy_re(', J, ') ', CYR(J)
         WRITE(*,'(A,I2,A,E25.17)') 'cy_im(', J, ') ', CYI(J)
      END DO
      WRITE(*,'(A)') 'END'
C
      RETURN
      END
C
C     ============================================================
C     Helper: call ZBKNU and print results
C     ============================================================
C
      SUBROUTINE BKNU_TEST(ZR, ZI, FNU, KODE, N,
     *                     TOL, ELIM, ALIM, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU, TOL, ELIM, ALIM
      INTEGER KODE, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION YR(10), YI(10)
      INTEGER NZ, J
C
      CALL ZBKNU(ZR, ZI, FNU, KODE, N, YR, YI, NZ, TOL, ELIM, ALIM)
C
      WRITE(*,'(A,A)') 'BEGIN zbknu ', LABEL
      WRITE(*,'(A,E25.17)') 'z_re  ', ZR
      WRITE(*,'(A,E25.17)') 'z_im  ', ZI
      WRITE(*,'(A,E25.17)') 'fnu   ', FNU
      WRITE(*,'(A,I6)')     'kode  ', KODE
      WRITE(*,'(A,I6)')     'n     ', N
      WRITE(*,'(A,I6)')     'nz    ', NZ
      DO J = 1, N
         WRITE(*,'(A,I2,A,E25.17)') 'y_re(', J, ') ', YR(J)
         WRITE(*,'(A,I2,A,E25.17)') 'y_im(', J, ') ', YI(J)
      END DO
      WRITE(*,'(A)') 'END'
C
      RETURN
      END
