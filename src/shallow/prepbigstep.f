c
c-------------------------------------------------------------------------------------
c
       subroutine prepbigstep(nvar,naux,lcheck,mptr,nx,ny,midub,mjdub,
     .                       valbgc,auxbgc,mi2tot,mj2tot)

       use amr_module, only: rnode,cornylo,cornxlo,evol,hxposs,hyposs
       use amr_module, only: NEEDS_TO_BE_SET,nghost,store2,timemult
       use amr_module, only: timeSetaux,possk,alloc,node,timeSetauxCPU
       implicit none

       integer, intent(in) :: nvar, naux, lcheck, mptr, nx, ny
       integer, intent(in) :: midub, mjdub, mi2tot, mj2tot
       real(kind=8), intent(inout) :: valbgc(nvar,mi2tot,mj2tot)
       real(kind=8), intent(inout) :: auxbgc(naux,mi2tot,mj2tot)

c      # Local variables
       real(kind=8) :: valdub(nvar,midub,mjdub)
       real(kind=8) :: auxdub(naux,midub,mjdub)
       real(kind=8) :: fp(nvar,mi2tot,mj2tot),gp(nvar,mi2tot,mj2tot)
       real(kind=8) :: fm(nvar,mi2tot,mj2tot),gm(nvar,mi2tot,mj2tot)
       real(kind=8) :: hx,hy,hx2,hy2,dt,dt2,dtnew2,time,tpre
       integer :: mitot,mjtot,ng2,locold,mx,my
       real(kind=8) :: xlow,ylow,xl,yb

       !for setaux timing
       integer :: clock_start, clock_finish, clock_rate
       real(kind=8) :: cpu_start, cpu_finish

       write(*,*) "In prepbigstep. size of valbgc: ", size(valbgc)
       write(*,*) "In prepbigstep. size of auxbgc: ", size(auxbgc)
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

          write(*,*) "In prepbigstep 2. size of valbgc: ", size(valbgc)
          write(*,*) "In prepbigstep 2. size of auxbgc: ", size(auxbgc)
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

          write(*,*) "In prepbigstep 3. size of valbgc: ", size(valbgc)
          write(*,*) "In prepbigstep 3. size of auxbgc: ", size(auxbgc)

c         update counts for error estimation work
          evol = evol + (nx/2)*(ny/2)

          write(*,*) "Leaving prepbigstep"
           return
           end
