C     Generate reference values for I Bessel function tests.
C     Calls ZBESI (upper interface) from TOMS 644 (zbsubs.f)
C     with IEEE 754 machine constants.
C
C     Output format: structured text parsed by txt_to_json.py
C     Each test case is delimited by BEGIN/END markers.
C
C     Usage: compile with zbsubs.f + machcon_ieee.f, then run.
C
      PROGRAM GEN_BESI_REFS
      IMPLICIT NONE
C
      WRITE(*,'(A)') '# Reference values for complex-bessel I tests'
      WRITE(*,'(A)') '# Generator: TOMS 644 zbsubs.f (rev 930101)'
      WRITE(*,'(A)') '# Machine: IEEE 754 double precision'
      WRITE(*,'(A)') '#'
C
C     ============================================================
C     SECTION 1: Real axis tests
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbesi'
C
C     I_0(1.0)
      CALL BESI_TEST(1.0D0, 0.0D0, 0.0D0, 1, 1, 'I_0(1.0)')
C     I_1(1.0)
      CALL BESI_TEST(1.0D0, 0.0D0, 1.0D0, 1, 1, 'I_1(1.0)')
C     I_0(0.5)
      CALL BESI_TEST(0.5D0, 0.0D0, 0.0D0, 1, 1, 'I_0(0.5)')
C     I_0(5.0)
      CALL BESI_TEST(5.0D0, 0.0D0, 0.0D0, 1, 1, 'I_0(5.0)')
C     I_0(10.0)
      CALL BESI_TEST(10.0D0, 0.0D0, 0.0D0, 1, 1, 'I_0(10.0)')
C     I_{0.5}(2.0), half-integer order
      CALL BESI_TEST(2.0D0, 0.0D0, 0.5D0, 1, 1, 'I_0.5(2.0)')
C     I_{2.5}(3.0), fractional order
      CALL BESI_TEST(3.0D0, 0.0D0, 2.5D0, 1, 1, 'I_2.5(3.0)')
C
C     ============================================================
C     SECTION 2: Complex argument tests
C     ============================================================
C
C     I_0(1+i)
      CALL BESI_TEST(1.0D0, 1.0D0, 0.0D0, 1, 1, 'I_0(1+i)')
C     I_0(1+2i)
      CALL BESI_TEST(1.0D0, 2.0D0, 0.0D0, 1, 1, 'I_0(1+2i)')
C     I_1(2+i)
      CALL BESI_TEST(2.0D0, 1.0D0, 1.0D0, 1, 1, 'I_1(2+i)')
C     I_0(0.1+0.1i), small z
      CALL BESI_TEST(0.1D0, 0.1D0, 0.0D0, 1, 1, 'I_0(0.1+0.1i)')
C     I_{0.25}(1.5+0.5i)
      CALL BESI_TEST(1.5D0, 0.5D0, 0.25D0, 1, 1, 'I_0.25(1.5+0.5i)')
C
C     ============================================================
C     SECTION 3: Negative real axis tests (analytic continuation)
C     ============================================================
C
C     I_0(-1.0)
      CALL BESI_TEST(-1.0D0, 0.0D0, 0.0D0, 1, 1, 'I_0(-1.0)')
C     I_1(-2.0)
      CALL BESI_TEST(-2.0D0, 0.0D0, 1.0D0, 1, 1, 'I_1(-2.0)')
C     I_0(-1+i)
      CALL BESI_TEST(-1.0D0, 1.0D0, 0.0D0, 1, 1, 'I_0(-1+i)')
C
C     ============================================================
C     SECTION 4: Scaled tests (KODE=2)
C     ============================================================
C
C     I_0(5.0) scaled
      CALL BESI_TEST(5.0D0, 0.0D0, 0.0D0, 2, 1, 'I_0(5.0) scaled')
C     I_0(1+i) scaled
      CALL BESI_TEST(1.0D0, 1.0D0, 0.0D0, 2, 1, 'I_0(1+i) scaled')
C
C     ============================================================
C     SECTION 5: Sequence (N>1) tests
C     ============================================================
C
C     I_{0,1,2}(2.0)
      CALL BESI_TEST(2.0D0, 0.0D0, 0.0D0, 1, 3, 'I_{0,1,2}(2.0)')
C     I_{0,1,2}(1+i)
      CALL BESI_TEST(1.0D0, 1.0D0, 0.0D0, 1, 3, 'I_{0,1,2}(1+i)')
C     I_{0.5,1.5,2.5}(3+i)
      CALL BESI_TEST(3.0D0, 1.0D0, 0.5D0, 1, 3,
     *               'I_{0.5,1.5,2.5}(3+i)')
C
      STOP
      END
C
C     ============================================================
C     Helper: call ZBESI and print results
C     ============================================================
C
      SUBROUTINE BESI_TEST(ZR, ZI, FNU, KODE, N, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU
      INTEGER KODE, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION CYR(10), CYI(10)
      INTEGER NZ, IERR, J
C
      CALL ZBESI(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zbesi ', LABEL
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
