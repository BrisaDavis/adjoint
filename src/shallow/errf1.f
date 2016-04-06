c
c --------------------------------------------------------------
c
      subroutine errf1(rctfine,nvar,rctcrse,mptr,mi2tot,mj2tot,
     2                 mitot,mjtot,rctflg,mibuff,mjbuff,auxfine,
     2                 naux,auxcrse)

      use amr_module
      use innerprod_module, only: calculate_max_innerproduct
      use adjoint_module, only: innerprod_index
      implicit none

      integer, intent(in) :: nvar, mptr, mi2tot, mj2tot
      integer, intent(in) :: mitot, mjtot,mibuff,mjbuff, naux
      real(kind=8), intent(in) :: rctfine(nvar,mitot,mjtot)
      real(kind=8), intent(inout) :: rctcrse(nvar,mi2tot,mj2tot)
      real(kind=8), intent(inout) :: auxfine(naux,mitot,mjtot)
      real(kind=8), intent(in) :: auxcrse(naux,mi2tot,mj2tot)
      real(kind=8), intent(inout) :: rctflg(mibuff,mjbuff)

      real(kind=8) :: qerr(3)

      logical :: interp2, interp3, interp4
      real(kind=8) :: errmax, err2, order, etaerr, terms
      integer :: i, j, ifine, jfine, jj, ii, levm, m
      real(kind=8) :: eta1, eta2, eta3, eta4, etacrse, est
      real(kind=8) :: term1, term2, term3, term4, aval
      real(kind=8) :: aux_a, aux_crse, aux_fine
      real(kind=8) :: hx, hy, dt, rflag, time
      real(kind=8) :: xofi, ybot, xleft, yofj
c
c
c ::::::::::::::::::::::::::::: ERRF1 ::::::::::::::::::::::::::::::::
c
c  Richardson error estimator:  Used when flag_richardson is .true.
c  Compare error estimates in rctfine, rctcrse, 
c  A point is flagged if the error estimate is greater than tol
c  later we check if its in a region where its allowed to be flagged
c  or alternatively required.
c
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

c
      time  = rnode(timemult, mptr)
      xleft = rnode(cornxlo,mptr)
      levm  = node(nestlevel, mptr)
      hx    = hxposs(levm)
      ybot  = rnode(cornylo,mptr)
      hy    = hyposs(levm)
      dt    = possk(levm)
 
      errmax = 0.0d0
      err2   = 0.0d0
c     order  = dt*dble(2**(iorder+1) - 2)
      order  = dble(2**(iorder+1) - 2)
c
c     Print values for debugging
      if (.not. (edebug)) go to 20
         write(outunit,107) mptr
 107     format(//,' coarsened grid values for grid ',i4)
         do 10 jj = nghost+1, mj2tot-nghost
            j = mj2tot + 1 - jj
            write(outunit,101) (rctcrse(1,i,j),
     .                          i = nghost+1, mi2tot-nghost)
10       continue
         write(outunit,108) mptr
 108     format(//, ' fine grid values for grid ',i4)
         do 15 jj = nghost+1, mjtot-nghost
            j = mjtot + 1 - jj
            write(outunit,101) (rctfine(1,i,j),i=nghost+1,mitot-nghost)
15       continue
101      format(' ',10e15.7)
c
c zero out the exterior locations so they don't affect err.est.
c
 20   continue
      jfine = nghost+1
      do 35  j = nghost+1, mj2tot-nghost
      yofj  = ybot + (dble(jfine) - .5d0)*hy
      ifine = nghost+1
c
      do 30  i  = nghost+1, mi2tot-nghost
          rflag = goodpt
          xofi  = xleft + (dble(ifine) - .5d0)*hx
c
c         Calculate error for each value of q and eta
c         if course and fine grids are in same wet/dry state
c
        aux_fine = auxfine(1,ifine,jfine)
        aux_crse = auxcrse(1,i,j)
        if(sign(aux_fine, aux_crse) .ne. sign(aux_fine, aux_fine))then
          go to 33
        else
          terms = 1.d0
          do 40 m = 1,nvar
              term1 = rctfine(m,ifine,jfine)

c             If the adjacent cell isn't in the same wet/dry
c             state, don't use it to calculate the error
              aux_a = auxfine(1,ifine+1,jfine)
              if(sign(aux_fine, aux_a) .ne.
     .                          sign(aux_fine, aux_fine))then
                  term2 = 0.d0
                  interp2 = .false.
              else
                  terms = terms + 1
                  term2 = rctfine(m,ifine+1,jfine)
                  interp2 = .true.
              endif

              aux_a = auxfine(1,ifine+1,jfine+1)
              if(sign(aux_fine, aux_a) .ne.
     .                          sign(aux_fine, aux_fine))then
                  term3 = 0.d0
                  interp3 = .false.
              else
                  terms = terms + 1
                  term3 = rctfine(m,ifine+1,jfine+1)
                  interp3 = .true.
              endif

              aux_a = auxfine(1,ifine,jfine+1)
              if(sign(aux_fine, aux_a) .ne.
     .                          sign(aux_fine, aux_fine))then
                  term4 = 0.d0
                  interp4 = .false.
              else
                  terms = terms + 1
                  term4 = rctfine(m,ifine,jfine+1)
                  interp4 = .true.
              endif

c             # divide by (aval*order) for relative error
              aval  = (term1+term2+term3+term4)/terms
              qerr(m)   =  dabs((aval-rctcrse(m,i,j))/ order)
40        continue

          est = qerr(1)
          if (est .gt. errmax) errmax = est
          err2 = err2 + est*est
c          write(outunit,102) i,j,est,rctcrse(1,i,j)
 102      format(' i,j,est ',2i5,2e15.7)
          write(outunit,104) term1,term2,term3,term4
 104      format('   ',4e15.7)
c         rctcrse(2,i,j) = est
c

c         Calculate error for eta
c
          eta1 = rctfine(1,ifine,jfine) + auxfine(1,ifine,jfine)
          if (abs(eta1) < 1d-90) then
              eta1 = 0.d0
          end if

          eta2 = 0.d0
          eta3 = 0.d0
          eta4 = 0.d0
          if (interp2) eta2 = rctfine(1,ifine+1,jfine) +
     .                          auxfine(1,ifine+1,jfine)
          if (abs(eta2) < 1d-90) then
              eta2 = 0.d0
          end if

          if (interp3) eta3 = rctfine(1,ifine+1,jfine+1) +
     .                          auxfine(1,ifine+1,jfine+1)
          if (abs(eta3) < 1d-90) then
              eta3 = 0.d0
          end if

          if (interp4) eta4 = rctfine(1,ifine,jfine+1) +
     .                          auxfine(1,ifine,jfine+1)
          if (abs(eta4) < 1d-90) then
              eta4 = 0.d0
          end if

c         # divide by (aval*order) for relative error
          aval  = (eta1+eta2+eta3+eta4)/terms
          etacrse = rctcrse(1,i,j) + auxcrse(1,i,j)
          if (abs(etacrse) < 1d-90) then
              etacrse = 0.d0
          end if
          etaerr   =  dabs((aval-etacrse)/ order)

c         Set innerproduct for fine grid
          auxfine(innerprod_index,ifine,jfine) =
     .           calculate_max_innerproduct(time,xofi,yofj,etaerr,
     .           qerr(2),qerr(3),auxfine(1,ifine,jfine))

          if (auxfine(innerprod_index,ifine,jfine) >= 0.05) then
          if ((xofi >= 200) .and. (xofi <=210)) then
          if ((yofj >= 10) .and. (yofj <=20)) then
        write(*,*) "TESTING:"
        write(*,*) xofi,yofj, etaerr
        write(*,*) eta1, eta2, eta3, eta4, etacrse
        write(*,*) aux_fine, aux_crse
        write(*,*) rctfine(1,ifine,jfine), rctcrse(1,i,j)
        write(*,*) "Innerprod: ", auxfine(innerprod_index,ifine,jfine)
          endif
          endif
          endif

          if (interp2) auxfine(innerprod_index,ifine+1,jfine)  =
     .                          auxfine(innerprod_index,ifine,jfine)
          if (interp3) auxfine(innerprod_index,ifine+1,jfine+1)  =
     .                          auxfine(innerprod_index,ifine,jfine)
          if (interp4) auxfine(innerprod_index,ifine,jfine+1)=
     .                          auxfine(innerprod_index,ifine,jfine)

          if (auxfine(innerprod_index,ifine,jfine) .ge. tol) then
             rflag  = badpt
          endif 
      rctcrse(1,i,j) = rflag
      endif
 33   ifine = ifine + 2
 30   continue
      jfine = jfine + 2
 35   continue
c
c  print out intermediate flagged rctcrse (for debugging)
c
      if (eprint) then
         err2 = dsqrt(err2/dble((mi2tot-2*nghost)*(mj2tot-2*nghost)))
         write(outunit,103) mptr, levm, errmax, err2
 103     format(' grid ',i4,' level ',i4,
     .          ' max. error = ',e15.7,' err2 = ',e15.7)
         if (edebug) then
           write(outunit,*) ' flagged points on coarsened grid ',
     .                      '(no ghost cells) for grid ',mptr
           do 45 jj = nghost+1, mj2tot-nghost
              j = mj2tot + 1 - jj
              write(outunit,106) (nint(rctcrse(1,i,j)),
     .                            i=nghost+1,mi2tot-nghost)
106           format(1h ,80i1)
45         continue
         endif
      endif
c
c     If coarse point is flagged, set point in flag array
      jfine   = nghost+1
      do 70 j = nghost+1, mj2tot-nghost
      ifine   = nghost+1
      do 60 i = nghost+1, mi2tot-nghost
         if (rctcrse(1,i,j) .eq. goodpt) go to 55
c           ## never set rctflg to good, since flag2refine may
c           ## have previously set it to bad
c           ## can only add bad pts in this routine
            rctflg(ifine,jfine)    = badpt
            rctflg(ifine+1,jfine)  = badpt
            rctflg(ifine,jfine+1)  = badpt
            rctflg(ifine+1,jfine+1)= badpt
 55       ifine   = ifine + 2
 60     continue
        jfine   = jfine + 2
 70   continue
c

      if (eprint) then
         write(outunit,118)
 118     format(' on fine grid (no ghost cells) flagged points are')
         if (edebug) then
          do 56 jj = nghost+1, mjtot-nghost
           j = mjtot + 1 - jj
           write(outunit,106)
     &      (nint(rctflg(i,j)),i=nghost+1,mitot-nghost)
 56       continue
        endif
      endif

      return
      end
