c
c-------------------------------------------------------------------------------------
c
       subroutine prepbigstep(nvar,naux,lcheck,mptr,nx,ny,midub,mjdub,
     .                       valbgc,auxbgc,mi2tot,mj2tot)

       use amr_module
       implicit double precision (a-h,o-z)

       double precision valdub(nvar,midub,mjdub)
       double precision auxdub(naux,midub,mjdub)
       double precision valbgc(nvar,mi2tot,mj2tot)
       double precision auxbgc(naux,mi2tot,mj2tot)
       dimension fp(nvar,mi2tot,mj2tot),gp(nvar,mi2tot,mj2tot)
       dimension fm(nvar,mi2tot,mj2tot),gm(nvar,mi2tot,mj2tot)

       !for setaux timing
       integer :: clock_start, clock_finish, clock_rate
       real(kind=8) :: cpu_start, cpu_finish

       write(*,*) "In prepbigstep"
          hx  = hxposs(lcheck)
          hy  = hyposs(lcheck)
          hx2 = 2.d0*hx
          hy2 = 2.d0*hy
          dt  = possk(lcheck)
          dt2 = 2. * dt
          time  = rnode(timemult,mptr)
          tpre  = time - dt

          mitot  = nx + 2*nghost
          mjtot  = ny + 2*nghost
          ng2    = 2*nghost
          locold = node(store2,mptr)
          xlow   = rnode(cornxlo,mptr) - nghost*hx2
          ylow   = rnode(cornylo,mptr) - nghost*hy2
c
       write(*,*) "About to call copysol"
c         # transfer soln. into grid with twice the ghost cells
          call copysol(valdub,alloc(locold),nvar,mitot,mjtot,
     1              nghost,midub,mjdub,ng2)

c
          if (naux .gt. 0) then
              xl     = rnode(cornxlo, mptr)
              yb     = rnode(cornylo, mptr)
              mx = midub - 4*nghost
              my = mjdub - 4*nghost
              auxdub = NEEDS_TO_BE_SET  ! signal that needs a val

c             # Generate aux for fine grid
              call system_clock(clock_start, clock_rate)
              call cpu_time(cpu_start)
              call setaux(2*nghost,mx,my,xl,yb,hx,hy,
     &                    naux,auxdub)
              call system_clock(clock_finish, clock_rate)
              call cpu_time(cpu_finish)
              timeSetaux = timeSetaux + clock_finish - clock_start
              timeSetauxCPU = timeSetauxCPU + cpu_finish - cpu_start

c             # Generate aux for coarse grid
              auxbgc = NEEDS_TO_BE_SET  ! signal that needs a val
              call setaux(2*nghost,mx/2,my/2,xl,yb,hx2,hy2,
     &                    naux,auxbgc)
          endif
          write(*,*) "Finished setting aux values"

c         # fill it - use enlarged (before coarsening) aux arrays
          call bound(tpre,nvar,ng2,valdub,midub,mjdub,mptr,
     1               auxdub,naux)
          write(*,*) "Finished calling bound"

c         coarsen by 2 in every direction
          call coarsen(valdub,midub,mjdub,auxdub,
     1                 valbgc,mi2tot,mj2tot,auxbgc,nvar,naux)
          write(*,*) "Finished calling coarsen"

          call stepgrid(valbgc,fm,fp,gm,gp,
     1                mi2tot,mj2tot,nghost,
     2                dt2,dtnew2,hx2,hy2,nvar,
     3                xlow,ylow,tpre,mptr,naux,auxbgc)
          write(*,*) "Finished stepping grid"

c         update counts for error estimation work
          evol = evol + (nx/2)*(ny/2)

          write(*,*) "Leaving prepbigstep"
           return
           end
