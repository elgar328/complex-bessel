C     Generate reference values for Y Bessel function tests.
C     Calls ZBESY (upper interface) from TOMS 644 (zbsubs.f)
C     with IEEE 754 machine constants.
C
C     Output format: structured text parsed by txt_to_json.py
C     Each test case is delimited by BEGIN/END markers.
C
C     Usage: compile with zbsubs.f + machcon_ieee.f, then run.
C
      PROGRAM GEN_BESY_REFS
      IMPLICIT NONE
C
      WRITE(*,'(A)') '# Reference values for complex-bessel Y tests'
      WRITE(*,'(A)') '# Generator: TOMS 644 zbsubs.f (rev 930101)'
      WRITE(*,'(A)') '# Machine: IEEE 754 double precision'
      WRITE(*,'(A)') '#'
C
C     ============================================================
C     SECTION 1: Real axis tests
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbesy'
C
C     Y_0(1.0)
      CALL BESY_TEST(1.0D0, 0.0D0, 0.0D0, 1, 1, 'Y_0(1.0)')
C     Y_1(1.0)
      CALL BESY_TEST(1.0D0, 0.0D0, 1.0D0, 1, 1, 'Y_1(1.0)')
C     Y_0(0.5)
      CALL BESY_TEST(0.5D0, 0.0D0, 0.0D0, 1, 1, 'Y_0(0.5)')
C     Y_0(5.0)
      CALL BESY_TEST(5.0D0, 0.0D0, 0.0D0, 1, 1, 'Y_0(5.0)')
C     Y_0(10.0)
      CALL BESY_TEST(10.0D0, 0.0D0, 0.0D0, 1, 1, 'Y_0(10.0)')
C     Y_{0.5}(2.0), half-integer order
      CALL BESY_TEST(2.0D0, 0.0D0, 0.5D0, 1, 1, 'Y_0.5(2.0)')
C     Y_{2.5}(3.0), fractional order
      CALL BESY_TEST(3.0D0, 0.0D0, 2.5D0, 1, 1, 'Y_2.5(3.0)')
C
C     ============================================================
C     SECTION 2: Complex argument tests
C     ============================================================
C
C     Y_0(1+i)
      CALL BESY_TEST(1.0D0, 1.0D0, 0.0D0, 1, 1, 'Y_0(1+i)')
C     Y_0(1+2i)
      CALL BESY_TEST(1.0D0, 2.0D0, 0.0D0, 1, 1, 'Y_0(1+2i)')
C     Y_1(2+i)
      CALL BESY_TEST(2.0D0, 1.0D0, 1.0D0, 1, 1, 'Y_1(2+i)')
C     Y_0(0.1+0.1i), small z
      CALL BESY_TEST(0.1D0, 0.1D0, 0.0D0, 1, 1, 'Y_0(0.1+0.1i)')
C     Y_{0.25}(1.5+0.5i)
      CALL BESY_TEST(1.5D0, 0.5D0, 0.25D0, 1, 1, 'Y_0.25(1.5+0.5i)')
C
C     ============================================================
C     SECTION 3: Lower half-plane
C     ============================================================
C
C     Y_0(1-i)
      CALL BESY_TEST(1.0D0, -1.0D0, 0.0D0, 1, 1, 'Y_0(1-i)')
C     Y_1(2-i)
      CALL BESY_TEST(2.0D0, -1.0D0, 1.0D0, 1, 1, 'Y_1(2-i)')
C
C     ============================================================
C     SECTION 4: Scaled tests (KODE=2)
C     ============================================================
C
C     Y_0(5.0) scaled
      CALL BESY_TEST(5.0D0, 0.0D0, 0.0D0, 2, 1, 'Y_0(5.0) scaled')
C     Y_0(1+i) scaled
      CALL BESY_TEST(1.0D0, 1.0D0, 0.0D0, 2, 1, 'Y_0(1+i) scaled')
C
C     ============================================================
C     SECTION 5: Sequence (N>1) tests
C     ============================================================
C
C     Y_{0,1,2}(2.0)
      CALL BESY_TEST(2.0D0, 0.0D0, 0.0D0, 1, 3, 'Y_{0,1,2}(2.0)')
C     Y_{0,1,2}(1+i)
      CALL BESY_TEST(1.0D0, 1.0D0, 0.0D0, 1, 3, 'Y_{0,1,2}(1+i)')
C     Y_{0.5,1.5,2.5}(3+i)
      CALL BESY_TEST(3.0D0, 1.0D0, 0.5D0, 1, 3,
     *               'Y_{0.5,1.5,2.5}(3+i)')
C
      STOP
      END
C
C     ============================================================
C     Helper: call ZBESY and print results
C     ============================================================
C
      SUBROUTINE BESY_TEST(ZR, ZI, FNU, KODE, N, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU
      INTEGER KODE, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION CYR(10), CYI(10)
      DOUBLE PRECISION CWRKR(10), CWRKI(10)
      INTEGER NZ, IERR, J
C
      CALL ZBESY(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR, CWRKI,
     *           IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zbesy ', LABEL
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
