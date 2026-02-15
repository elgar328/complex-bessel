C     Generate reference values for Hankel function tests.
C     Calls ZBESH (upper interface) from TOMS 644 (zbsubs.f)
C     with IEEE 754 machine constants.
C
C     Output format: structured text parsed by txt_to_json.py
C     Each test case is delimited by BEGIN/END markers.
C
C     Usage: compile with zbsubs.f + machcon_ieee.f, then run.
C
      PROGRAM GEN_BESH_REFS
      IMPLICIT NONE
C
      WRITE(*,'(A)') '# Reference values for complex-bessel H tests'
      WRITE(*,'(A)') '# Generator: TOMS 644 zbsubs.f (rev 930101)'
      WRITE(*,'(A)') '# Machine: IEEE 754 double precision'
      WRITE(*,'(A)') '#'
C
C     ============================================================
C     SECTION 1: ZBESH H^(1) tests (M=1, Im(z)>=0 or real axis)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbesh'
C
C     --- H^(1) Real axis ---
C
C     H^(1)_0(0.5)
      CALL BESH_TEST(0.5D0, 0.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(0.5)')
C     H^(1)_0(1.0)
      CALL BESH_TEST(1.0D0, 0.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(1.0)')
C     H^(1)_0(5.0)
      CALL BESH_TEST(5.0D0, 0.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(5.0)')
C     H^(1)_0(10.0)
      CALL BESH_TEST(10.0D0, 0.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(10.0)')
C     H^(1)_1(1.0)
      CALL BESH_TEST(1.0D0, 0.0D0, 1.0D0, 1, 1, 1,
     *               'H1_1(1.0)')
C     H^(1)_{0.5}(2.0), half-integer order
      CALL BESH_TEST(2.0D0, 0.0D0, 0.5D0, 1, 1, 1,
     *               'H1_0.5(2.0)')
C     H^(1)_{2.5}(3.0), fractional order
      CALL BESH_TEST(3.0D0, 0.0D0, 2.5D0, 1, 1, 1,
     *               'H1_2.5(3.0)')
C
C     --- H^(1) Upper half-plane ---
C
C     H^(1)_0(1+i)
      CALL BESH_TEST(1.0D0, 1.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(1+i)')
C     H^(1)_0(1+2i)
      CALL BESH_TEST(1.0D0, 2.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(1+2i)')
C     H^(1)_0(5+3i)
      CALL BESH_TEST(5.0D0, 3.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(5+3i)')
C     H^(1)_0(0.1+0.1i), small z
      CALL BESH_TEST(0.1D0, 0.1D0, 0.0D0, 1, 1, 1,
     *               'H1_0(0.1+0.1i)')
C     H^(1)_1(1+i)
      CALL BESH_TEST(1.0D0, 1.0D0, 1.0D0, 1, 1, 1,
     *               'H1_1(1+i)')
C     H^(1)_{0.25}(2+i)
      CALL BESH_TEST(2.0D0, 1.0D0, 0.25D0, 1, 1, 1,
     *               'H1_0.25(2+i)')
C
C     --- H^(1) Negative real axis (should work without ZACON) ---
C
C     H^(1)_0(-1.0)
      CALL BESH_TEST(-1.0D0, 0.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(-1.0)')
C     H^(1)_1(-2.0)
      CALL BESH_TEST(-2.0D0, 0.0D0, 1.0D0, 1, 1, 1,
     *               'H1_1(-2.0)')
C
C     --- H^(1) Pure imaginary ---
C
C     H^(1)_0(2i)
      CALL BESH_TEST(0.0D0, 2.0D0, 0.0D0, 1, 1, 1,
     *               'H1_0(2i)')
C
C     ============================================================
C     SECTION 2: ZBESH H^(2) tests (M=2, Im(z)<=0 or real axis)
C     ============================================================
C
C     --- H^(2) Real axis ---
C
C     H^(2)_0(1.0)
      CALL BESH_TEST(1.0D0, 0.0D0, 0.0D0, 1, 2, 1,
     *               'H2_0(1.0)')
C     H^(2)_0(5.0)
      CALL BESH_TEST(5.0D0, 0.0D0, 0.0D0, 1, 2, 1,
     *               'H2_0(5.0)')
C     H^(2)_{0.5}(2.0)
      CALL BESH_TEST(2.0D0, 0.0D0, 0.5D0, 1, 2, 1,
     *               'H2_0.5(2.0)')
C
C     --- H^(2) Lower half-plane ---
C
C     H^(2)_0(1-i)
      CALL BESH_TEST(1.0D0, -1.0D0, 0.0D0, 1, 2, 1,
     *               'H2_0(1-i)')
C     H^(2)_0(5-3i)
      CALL BESH_TEST(5.0D0, -3.0D0, 0.0D0, 1, 2, 1,
     *               'H2_0(5-3i)')
C     H^(2)_0(3-i)
      CALL BESH_TEST(3.0D0, -1.0D0, 0.0D0, 1, 2, 1,
     *               'H2_0(3-i)')
C     H^(2)_1(1-i)
      CALL BESH_TEST(1.0D0, -1.0D0, 1.0D0, 1, 2, 1,
     *               'H2_1(1-i)')
C     H^(2)_{0.25}(2-i)
      CALL BESH_TEST(2.0D0, -1.0D0, 0.25D0, 1, 2, 1,
     *               'H2_0.25(2-i)')
C
C     --- H^(2) Pure imaginary ---
C
C     H^(2)_0(-2i)
      CALL BESH_TEST(0.0D0, -2.0D0, 0.0D0, 1, 2, 1,
     *               'H2_0(-2i)')
C
C     ============================================================
C     SECTION 3: Scaled (KODE=2) tests
C     ============================================================
C
C     H^(1)_0(1+i) scaled
      CALL BESH_TEST(1.0D0, 1.0D0, 0.0D0, 2, 1, 1,
     *               'H1_0(1+i) scaled')
C     H^(2)_0(1-i) scaled
      CALL BESH_TEST(1.0D0, -1.0D0, 0.0D0, 2, 2, 1,
     *               'H2_0(1-i) scaled')
C     H^(1)_0(5.0) scaled
      CALL BESH_TEST(5.0D0, 0.0D0, 0.0D0, 2, 1, 1,
     *               'H1_0(5.0) scaled')
C
C     ============================================================
C     SECTION 4: Sequence (N>1) tests
C     ============================================================
C
C     H^(1)_{0,1,2}(2.0), sequence of 3
      CALL BESH_TEST(2.0D0, 0.0D0, 0.0D0, 1, 1, 3,
     *               'H1_{0,1,2}(2.0)')
C     H^(2)_{0,1,2}(2.0), sequence of 3
      CALL BESH_TEST(2.0D0, 0.0D0, 0.0D0, 1, 2, 3,
     *               'H2_{0,1,2}(2.0)')
C     H^(1)_{0.5,1.5,2.5}(3+i)
      CALL BESH_TEST(3.0D0, 1.0D0, 0.5D0, 1, 1, 3,
     *               'H1_{0.5,1.5,2.5}(3+i)')
C
      STOP
      END
C
C     ============================================================
C     Helper: call ZBESH and print results
C     ============================================================
C
      SUBROUTINE BESH_TEST(ZR, ZI, FNU, KODE, M, N, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU
      INTEGER KODE, M, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION CYR(10), CYI(10)
      INTEGER NZ, IERR, J
C
      CALL ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zbesh ', LABEL
      WRITE(*,'(A,E25.17)') 'z_re  ', ZR
      WRITE(*,'(A,E25.17)') 'z_im  ', ZI
      WRITE(*,'(A,E25.17)') 'fnu   ', FNU
      WRITE(*,'(A,I6)')     'kode  ', KODE
      WRITE(*,'(A,I6)')     'm     ', M
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
