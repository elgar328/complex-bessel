C     Generate reference values for Phase 5b large-order tests.
C     Tests K, H, I, J functions at large orders (fnu > fnul ~ 86)
C     to exercise ZBUNK, ZBUNI, ZUNK2, ZUNI2, ZUNHJ paths.
C
C     Output format: structured text parsed by txt_to_json.py
C     Each test case is delimited by BEGIN/END markers.
C
C     Usage: compile with zbsubs.f + machcon_ieee.f, then run.
C
      PROGRAM GEN_PHASE5B_REFS
      IMPLICIT NONE
C
      WRITE(*,'(A)') '# Reference values for Phase 5b large-order'
      WRITE(*,'(A)') '# Generator: TOMS 644 zbsubs.f (rev 930101)'
      WRITE(*,'(A)') '# Machine: IEEE 754 double precision'
      WRITE(*,'(A)') '#'
C
C     ============================================================
C     SECTION 1: Large-order K function (ZBUNK path)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: large_order_K'
C
C     K_100(5+3i) — right half plane, region 1 (ZUNK1)
      CALL BESK_LO(5.0D0, 3.0D0, 100.0D0, 1, 1, 'K_100(5+3i)')
C     K_100(5+3i) scaled
      CALL BESK_LO(5.0D0, 3.0D0, 100.0D0, 2, 1, 'K_100(5+3i) sc')
C     K_90(10+5i)
      CALL BESK_LO(10.0D0, 5.0D0, 90.0D0, 1, 1, 'K_90(10+5i)')
C     K_100(1+10i) — region 2 (ZUNK2): |Im|>|Re|*sqrt(3)
      CALL BESK_LO(1.0D0, 10.0D0, 100.0D0, 1, 1, 'K_100(1+10i)')
C     K_100(1+10i) scaled
      CALL BESK_LO(1.0D0, 10.0D0, 100.0D0, 2, 1, 'K_100(1+10i) sc')
C     K_100(-5+3i) — left half plane: analytic continuation
      CALL BESK_LO(-5.0D0, 3.0D0, 100.0D0, 1, 1, 'K_100(-5+3i)')
C     K_100(-5-3i) — left half plane, Im < 0
      CALL BESK_LO(-5.0D0, -3.0D0, 100.0D0, 1, 1, 'K_100(-5-3i)')
C     K_90(-1+10i) — left half plane, region 2
      CALL BESK_LO(-1.0D0, 10.0D0, 90.0D0, 1, 1, 'K_90(-1+10i)')
C     K_{100,101}(5+3i) — sequence
      CALL BESK_LO(5.0D0, 3.0D0, 100.0D0, 1, 2, 'K_{100,101}(5+3i)')
C     K_100(50+30i) — large z
      CALL BESK_LO(50.0D0, 30.0D0, 100.0D0, 1, 1, 'K_100(50+30i)')
C
C     ============================================================
C     SECTION 2: Large-order H function (ZBUNK via ZBESH)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: large_order_H'
C
C     H^(1)_100(5+3i)
      CALL BESH_LO(5.0D0, 3.0D0, 100.0D0, 1, 1, 1, 'H1_100(5+3i)')
C     H^(2)_100(5-3i)
      CALL BESH_LO(5.0D0, -3.0D0, 100.0D0, 1, 2, 1, 'H2_100(5-3i)')
C     H^(1)_100(1+10i) — region 2
      CALL BESH_LO(1.0D0, 10.0D0, 100.0D0, 1, 1, 1,
     *             'H1_100(1+10i)')
C     H^(1)_100(5+3i) scaled
      CALL BESH_LO(5.0D0, 3.0D0, 100.0D0, 2, 1, 1,
     *             'H1_100(5+3i) sc')
C     H^(1)_{100,101}(5+3i) — sequence
      CALL BESH_LO(5.0D0, 3.0D0, 100.0D0, 1, 1, 2,
     *             'H1_{100,101}(5+3i)')
C     H^(1)_90(10+5i)
      CALL BESH_LO(10.0D0, 5.0D0, 90.0D0, 1, 1, 1, 'H1_90(10+5i)')
C
C     ============================================================
C     SECTION 3: Large-order I function (ZBUNI via ZBINU/ZBESI)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: large_order_I'
C
C     I_100(50+30i)
      CALL BESI_LO(50.0D0, 30.0D0, 100.0D0, 1, 1, 'I_100(50+30i)')
C     I_100(50+30i) scaled
      CALL BESI_LO(50.0D0, 30.0D0, 100.0D0, 2, 1,
     *             'I_100(50+30i) sc')
C     I_90(40+20i)
      CALL BESI_LO(40.0D0, 20.0D0, 90.0D0, 1, 1, 'I_90(40+20i)')
C     I_100(10+80i) — region 2: |Im|>|Re|*sqrt(3)
      CALL BESI_LO(10.0D0, 80.0D0, 100.0D0, 1, 1, 'I_100(10+80i)')
C     I_{100,101}(50+30i) — sequence
      CALL BESI_LO(50.0D0, 30.0D0, 100.0D0, 1, 2,
     *             'I_{100,101}(50+30i)')
C     I_100(-50+30i) — negative real
      CALL BESI_LO(-50.0D0, 30.0D0, 100.0D0, 1, 1,
     *             'I_100(-50+30i)')
C
C     ============================================================
C     SECTION 4: Large-order J function (ZBUNI via ZBINU/ZBESJ)
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: large_order_J'
C
C     J_100(50+30i)
      CALL BESJ_LO(50.0D0, 30.0D0, 100.0D0, 1, 1, 'J_100(50+30i)')
C     J_100(50+30i) scaled
      CALL BESJ_LO(50.0D0, 30.0D0, 100.0D0, 2, 1,
     *             'J_100(50+30i) sc')
C     J_90(40+20i)
      CALL BESJ_LO(40.0D0, 20.0D0, 90.0D0, 1, 1, 'J_90(40+20i)')
C     J_{100,101}(50+30i) — sequence
      CALL BESJ_LO(50.0D0, 30.0D0, 100.0D0, 1, 2,
     *             'J_{100,101}(50+30i)')
C
C     ============================================================
C     SECTION 5: Large-order Y function
C     ============================================================
C
      WRITE(*,'(A)') '# SECTION: large_order_Y'
C
C     Y_100(50+30i)
      CALL BESY_LO(50.0D0, 30.0D0, 100.0D0, 1, 1, 'Y_100(50+30i)')
C     Y_90(40+20i)
      CALL BESY_LO(40.0D0, 20.0D0, 90.0D0, 1, 1, 'Y_90(40+20i)')
C
      STOP
      END
C
C     ============================================================
C     Helper: ZBESK for large order
C     ============================================================
C
      SUBROUTINE BESK_LO(ZR, ZI, FNU, KODE, N, LABEL)
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
C     Helper: ZBESH for large order
C     ============================================================
C
      SUBROUTINE BESH_LO(ZR, ZI, FNU, KODE, M, N, LABEL)
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
C
C     ============================================================
C     Helper: ZBESI for large order
C     ============================================================
C
      SUBROUTINE BESI_LO(ZR, ZI, FNU, KODE, N, LABEL)
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
C
C     ============================================================
C     Helper: ZBESJ for large order
C     ============================================================
C
      SUBROUTINE BESJ_LO(ZR, ZI, FNU, KODE, N, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU
      INTEGER KODE, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION CYR(10), CYI(10)
      INTEGER NZ, IERR, J
C
      CALL ZBESJ(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, IERR)
C
      WRITE(*,'(A,A)') 'BEGIN zbesj ', LABEL
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
C     Helper: ZBESY for large order
C     ============================================================
C
      SUBROUTINE BESY_LO(ZR, ZI, FNU, KODE, N, LABEL)
      IMPLICIT NONE
      DOUBLE PRECISION ZR, ZI, FNU
      INTEGER KODE, N
      CHARACTER*(*) LABEL
C
      DOUBLE PRECISION CYR(10), CYI(10), CWRKR(10), CWRKI(10)
      INTEGER NZ, IERR, J
C
      CALL ZBESY(ZR, ZI, FNU, KODE, N, CYR, CYI, NZ, CWRKR,
     *           CWRKI, IERR)
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
