C     Generate reference values for Airy function tests.
C     Calls ZAIRY and ZBIRY from TOMS 644 (zbsubs.f)
C     with IEEE 754 machine constants.
C
C     Output format: structured text parsed by txt_to_json.py
C     Each test case is delimited by BEGIN/END markers.
C
C     Usage: compile with zbsubs.f + machcon_ieee.f, then run.
C
      PROGRAM GEN_AIRY_REFS
      IMPLICIT NONE
C
      WRITE(*,'(A)') '# Reference values for complex-bessel Airy tests'
      WRITE(*,'(A)') '# Generator: TOMS 644 zbsubs.f (rev 930101)'
      WRITE(*,'(A)') '# Machine: IEEE 754 double precision'
      WRITE(*,'(A)') '#'
C
C     ============================================================
C     SECTION 1: ZAIRY - Ai(z) value (ID=0)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zairy value'
C
C     Ai(0) - zero argument
      CALL AIRY_TEST(0.0D0, 0.0D0, 0, 1, 'Ai(0)')
C     Ai(0.5) - small real, power series
      CALL AIRY_TEST(0.5D0, 0.0D0, 0, 1, 'Ai(0.5)')
C     Ai(1.0) - boundary |z|=1
      CALL AIRY_TEST(1.0D0, 0.0D0, 0, 1, 'Ai(1.0)')
C     Ai(2.0) - real, K function path
      CALL AIRY_TEST(2.0D0, 0.0D0, 0, 1, 'Ai(2.0)')
C     Ai(5.0) - larger real
      CALL AIRY_TEST(5.0D0, 0.0D0, 0, 1, 'Ai(5.0)')
C     Ai(10.0) - large real
      CALL AIRY_TEST(10.0D0, 0.0D0, 0, 1, 'Ai(10.0)')
C     Ai(-1.0) - negative real, ZACAI path
      CALL AIRY_TEST(-1.0D0, 0.0D0, 0, 1, 'Ai(-1.0)')
C     Ai(-5.0) - larger negative, oscillatory
      CALL AIRY_TEST(-5.0D0, 0.0D0, 0, 1, 'Ai(-5.0)')
C     Ai(1+i) - complex
      CALL AIRY_TEST(1.0D0, 1.0D0, 0, 1, 'Ai(1+i)')
C     Ai(1+2i) - complex, larger imag
      CALL AIRY_TEST(1.0D0, 2.0D0, 0, 1, 'Ai(1+2i)')
C     Ai(-1+i) - left half plane
      CALL AIRY_TEST(-1.0D0, 1.0D0, 0, 1, 'Ai(-1+i)')
C     Ai(-2-i) - left half, lower
      CALL AIRY_TEST(-2.0D0, -1.0D0, 0, 1, 'Ai(-2-i)')
C     Ai(0.1+0.1i) - small complex
      CALL AIRY_TEST(0.1D0, 0.1D0, 0, 1, 'Ai(0.1+0.1i)')
C     Ai(3+2i)
      CALL AIRY_TEST(3.0D0, 2.0D0, 0, 1, 'Ai(3+2i)')
C
C     ============================================================
C     SECTION 2: ZAIRY - Ai'(z) derivative (ID=1)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zairy derivative'
C
C     Ai'(0)
      CALL AIRY_TEST(0.0D0, 0.0D0, 1, 1, 'Aip(0)')
C     Ai'(1.0)
      CALL AIRY_TEST(1.0D0, 0.0D0, 1, 1, 'Aip(1.0)')
C     Ai'(-1.0)
      CALL AIRY_TEST(-1.0D0, 0.0D0, 1, 1, 'Aip(-1.0)')
C     Ai'(1+i)
      CALL AIRY_TEST(1.0D0, 1.0D0, 1, 1, 'Aip(1+i)')
C     Ai'(2.0)
      CALL AIRY_TEST(2.0D0, 0.0D0, 1, 1, 'Aip(2.0)')
C
C     ============================================================
C     SECTION 3: ZAIRY - scaled (KODE=2)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zairy scaled'
C
C     Ai(2.0) scaled
      CALL AIRY_TEST(2.0D0, 0.0D0, 0, 2, 'Ai(2.0) scaled')
C     Ai(1+i) scaled
      CALL AIRY_TEST(1.0D0, 1.0D0, 0, 2, 'Ai(1+i) scaled')
C     Ai(-1.0) scaled
      CALL AIRY_TEST(-1.0D0, 0.0D0, 0, 2, 'Ai(-1.0) scaled')
C     Ai'(2.0) scaled
      CALL AIRY_TEST(2.0D0, 0.0D0, 1, 2, 'Aip(2.0) scaled')
C
C     ============================================================
C     SECTION 4: ZBIRY - Bi(z) value (ID=0)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbiry value'
C
C     Bi(0)
      CALL BIRY_TEST(0.0D0, 0.0D0, 0, 1, 'Bi(0)')
C     Bi(0.5)
      CALL BIRY_TEST(0.5D0, 0.0D0, 0, 1, 'Bi(0.5)')
C     Bi(1.0)
      CALL BIRY_TEST(1.0D0, 0.0D0, 0, 1, 'Bi(1.0)')
C     Bi(2.0)
      CALL BIRY_TEST(2.0D0, 0.0D0, 0, 1, 'Bi(2.0)')
C     Bi(5.0)
      CALL BIRY_TEST(5.0D0, 0.0D0, 0, 1, 'Bi(5.0)')
C     Bi(-1.0)
      CALL BIRY_TEST(-1.0D0, 0.0D0, 0, 1, 'Bi(-1.0)')
C     Bi(-5.0)
      CALL BIRY_TEST(-5.0D0, 0.0D0, 0, 1, 'Bi(-5.0)')
C     Bi(1+i)
      CALL BIRY_TEST(1.0D0, 1.0D0, 0, 1, 'Bi(1+i)')
C     Bi(1+2i)
      CALL BIRY_TEST(1.0D0, 2.0D0, 0, 1, 'Bi(1+2i)')
C     Bi(-1+i)
      CALL BIRY_TEST(-1.0D0, 1.0D0, 0, 1, 'Bi(-1+i)')
C     Bi(-2-i)
      CALL BIRY_TEST(-2.0D0, -1.0D0, 0, 1, 'Bi(-2-i)')
C     Bi(0.1+0.1i)
      CALL BIRY_TEST(0.1D0, 0.1D0, 0, 1, 'Bi(0.1+0.1i)')
C     Bi(3+2i)
      CALL BIRY_TEST(3.0D0, 2.0D0, 0, 1, 'Bi(3+2i)')
C
C     ============================================================
C     SECTION 5: ZBIRY - Bi'(z) derivative (ID=1)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbiry derivative'
C
C     Bi'(0)
      CALL BIRY_TEST(0.0D0, 0.0D0, 1, 1, 'Bip(0)')
C     Bi'(1.0)
      CALL BIRY_TEST(1.0D0, 0.0D0, 1, 1, 'Bip(1.0)')
C     Bi'(-1.0)
      CALL BIRY_TEST(-1.0D0, 0.0D0, 1, 1, 'Bip(-1.0)')
C     Bi'(1+i)
      CALL BIRY_TEST(1.0D0, 1.0D0, 1, 1, 'Bip(1+i)')
C     Bi'(2.0)
      CALL BIRY_TEST(2.0D0, 0.0D0, 1, 1, 'Bip(2.0)')
C
C     ============================================================
C     SECTION 6: ZBIRY - scaled (KODE=2)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: zbiry scaled'
C
C     Bi(2.0) scaled
      CALL BIRY_TEST(2.0D0, 0.0D0, 0, 2, 'Bi(2.0) scaled')
C     Bi(1+i) scaled
      CALL BIRY_TEST(1.0D0, 1.0D0, 0, 2, 'Bi(1+i) scaled')
C     Bi(-1.0) scaled
      CALL BIRY_TEST(-1.0D0, 0.0D0, 0, 2, 'Bi(-1.0) scaled')
C     Bi'(2.0) scaled
      CALL BIRY_TEST(2.0D0, 0.0D0, 1, 2, 'Bip(2.0) scaled')
C
      STOP
      END
C
C     ============================================================
C     Helper: call ZAIRY and print results
C     ============================================================
C
      SUBROUTINE AIRY_TEST(ZR, ZI, ID, KODE, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI
      INTEGER ID, KODE
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION AIR, AII
      INTEGER NZ, IERR
C
      CALL ZAIRY(ZR, ZI, ID, KODE, AIR, AII, NZ, IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zairy ', LABEL
      WRITE(*,'(A,E25.17)') 'z_re  ', ZR
      WRITE(*,'(A,E25.17)') 'z_im  ', ZI
      WRITE(*,'(A,I6)')     'id    ', ID
      WRITE(*,'(A,I6)')     'kode  ', KODE
      WRITE(*,'(A,I6)')     'nz    ', NZ
      WRITE(*,'(A,I6)')     'ierr  ', IERR
      WRITE(*,'(A,E25.17)') 'cy_re( 1) ', AIR
      WRITE(*,'(A,E25.17)') 'cy_im( 1) ', AII
      WRITE(*,'(A)') 'END'
C
      RETURN
      END
C
C     ============================================================
C     Helper: call ZBIRY and print results
C     ============================================================
C
      SUBROUTINE BIRY_TEST(ZR, ZI, ID, KODE, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI
      INTEGER ID, KODE
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION BIR, BII
      INTEGER IERR
C
      CALL ZBIRY(ZR, ZI, ID, KODE, BIR, BII, IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zbiry ', LABEL
      WRITE(*,'(A,E25.17)') 'z_re  ', ZR
      WRITE(*,'(A,E25.17)') 'z_im  ', ZI
      WRITE(*,'(A,I6)')     'id    ', ID
      WRITE(*,'(A,I6)')     'kode  ', KODE
      WRITE(*,'(A,I6)')     'ierr  ', IERR
      WRITE(*,'(A,E25.17)') 'cy_re( 1) ', BIR
      WRITE(*,'(A,E25.17)') 'cy_im( 1) ', BII
      WRITE(*,'(A)') 'END'
C
      RETURN
      END
