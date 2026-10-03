#!/usr/bin/env python3
"""Instrument a COPY of the Fortran hpkvd source for the per-peak harness (MATCHED_FORTRAN_2026-10.md, step 3).

For peaks with |x+107|,|y+107|,|z+107| < 12 Mpc/h (core of tile (2,2,2) in configs/matched_fortran.toml) get_homel
writes, per MPI rank, dbg_homel_<rank>.txt:
  PK  x y z R_f ZZon ir2min
  RNG mlow mupp m            (shell range of the strain kernel at the first tested shell m0)
  SH  m1 nshell rad SR11 SR22 SR33 Sshell_x   (per-shell sums entering the kernel)
  M0  m0 rad Fbar Frhoh e_v p_v zvir1p ;  M0S Frhoc Fbar(m0-1) iflag
  MP  mp rad Fbar Frhoh e_v p_v zvir1p iflag  (each inward step)
  FIN m0 RTHL zvir1 zvir1p ZZon dZvir ;  R  RTHL [Mpc/h]
The Julia side prints the same lines through a temporary hook in RadialShell.analyse_peak (not committed; see
matched_peak_dbg.jl). Never apply this to ~/work/peakpatch itself.
  usage: python3 fortran_dbg_patch.py <copy>/src/hpkvd     then rebuild hpkvd in <copy>/src
"""
import sys, os
d = sys.argv[1]
rd = lambda f: open(os.path.join(d, f)).read()
def wr(f, s): open(os.path.join(d, f), 'w').write(s)
def sub(s, old, new):
    assert old in s, old[:60]
    return s.replace(old, new, 1)

a = rd('arrays.f90')
wr('arrays.f90', sub(a, 'contains', '  logical :: dbgpk = .false.\n  integer :: dbgu = 77\ncontains'))

p = rd('peakvoidsubs.f90')
p = sub(p, """  if(zvir1p.ge.ZZon) then ! zvir1p should be -1""",
"""  if(dbgpk) then
     write(dbgu,'(a,3i5)') 'RNG ',mlow,mupp,m
     do m1=mlow,mupp
        write(dbgu,'(a,2i5,5g16.8)') 'SH ',m1,nshell(m1),rad(m1),SRshell(1,1,m1),SRshell(2,2,m1),SRshell(3,3,m1),Sshell(1,m1)
     enddo
  endif
  if(dbgpk) write(dbgu,'(a,i4,6g16.8)') 'M0 ',m0,rad(m0),Fbar(m0),Frhoh,e_v,p_v,zvir1p
  if(dbgpk) write(dbgu,'(a,2g16.8,i4)') 'M0S ',Frhoc,Fbar(m0-1),iflag
  if(zvir1p.ge.ZZon) then ! zvir1p should be -1""")
p = sub(p, """     m0=mp+1
     if(zvir1p.ge.ZZon) goto 300 ! if peak will collapse go to 300""",
"""     m0=mp+1
     if(dbgpk) write(dbgu,'(a,i4,6g16.8,i3)') 'MP ',mp,rad(mp),Fbar(mp),Frhoh,e_v,p_v,zvir1p,iflag
     if(zvir1p.ge.ZZon) goto 300 ! if peak will collapse go to 300""")
p = sub(p, """  if(RTHL.le.0.0) then
     RTHL=-1.0
     Srb=0.0""",
"""  if(dbgpk) write(dbgu,'(a,i4,5g16.8)') 'FIN ',m0,RTHL,zvir1,zvir1p,ZZon,dZvir
  if(RTHL.le.0.0) then
     RTHL=-1.0
     Srb=0.0""")
wr('peakvoidsubs.f90', p)

h = rd('hpkvd.f90')
h = sub(h, """        ! HOMOGENEOUS ELLIPSOID CALCULATION
        call get_homel(npart,jp,alatt,ir2min,&""",
"""        dbgpk = abs(xpk(jpp,ired)+107.0)<12.0 .and. abs(ypk(jpp,ired)+107.0)<12.0 .and. abs(zpk(jpp,ired)+107.0)<12.0
        if(dbgpk) then
           if(dbgu==77) then
              dbgu = 100+myid
              open(unit=dbgu,file='dbg_homel_'//trim(str(myid))//'.txt')
           endif
           write(dbgu,'(a,5g16.8,i4)') 'PK ',xpk(jpp,ired),ypk(jpp,ired),zpk(jpp,ired),Rfclv(ic),1+redshiftpk,ir2min
        endif
        ! HOMOGENEOUS ELLIPSOID CALCULATION
        call get_homel(npart,jp,alatt,ir2min,&""")
h = sub(h, """        RTHLv(jpp,ired)=RTHL*alatt
""", """        RTHLv(jpp,ired)=RTHL*alatt
        if(dbgpk) then
           write(dbgu,'(a,g16.8)') 'R ',RTHLv(jpp,ired)
           flush(dbgu)
        endif
        dbgpk = .false.
""")
wr('hpkvd.f90', h)
print("instrumented", d)
