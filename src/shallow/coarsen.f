


c
c ---------------------------------------------------
c
       subroutine coarsen(valdub,midub,mjdub,auxdub,
     &                    valbgc,mi2tot,mj2tot,auxbgc,nvar,naux)


       use geoclaw_module, only: dry_tolerance
c     # modified from update.f to coarsen grid taking into
c     account wet and dry cells.
c
       use amr_module
       implicit double precision (a-h, o-z)

       dimension  valdub(nvar,midub, mjdub)
       dimension  auxdub(naux,midub, mjdub)
       dimension  valbgc(nvar,mi2tot,mj2tot)
       dimension  auxbgc(naux,mi2tot,mj2tot)

c :::::::::::::::::::::::: COARSEN ::::::::::::::::::::::::::::::::
c coarsen = coarsen the fine grid data (with double the usual
c           number of ghost cells to prepare coarsened
c           grid for error estimation.
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
       write(*,*) "In coarsen"
       do j = 1, mj2tot

         jfine = 2*(j-1) + 1

         do i = 1, mi2tot
           ifine = 2*(i-1) + 1

c          # Getting coarse capacity function and bathymetry
           if (mcapa .eq. 0) then
               capac=1.0d0
           else
               capac=auxbgc(2,i,j)
           endif

           bc = auxbgc(1,i,j)

c          # Calculating which values to use for coarsen
           etasum = 0.d0
           hsum = 0.d0
           husum = 0.d0
           hvsum = 0.d0

           nwet=0

           do jco = 1, 2
            do ico = 1, 2
               if (mcapa .eq. 0) then
                   capa=1.0d0
               else
                   capa=auxdub(2,ifine+ico-1,jfine+jco-1)
               endif

               hf =valdub(1,ifine+ico-1,jfine+jco-1)*capa
               huf=valdub(2,ifine+ico-1,jfine+jco-1)*capa
               hvf=valdub(3,ifine+ico-1,jfine+jco-1)*capa

               bf =auxdub(1,ifine+ico-1,jfine+jco-1)*capa

               if (hf > dry_tolerance) then
                   etaf = hf+bf
                   nwet=nwet+1
               else
                   etaf = 0.d0
                   huf=0.d0
                   hvf=0.d0
               endif

               hsum   = hsum + hf
               husum  = husum + huf
               hvsum  = hvsum + hvf
               etasum = etasum + etaf
             enddo
           enddo

c          # Getting coarse averages
           if (nwet.gt.0) then
               etaav=etasum/dble(nwet)
               hav= hsum/dble(nwet)
*              hc=max(etaav-bc*capac,0.d0) !tsunamiclaw method
               hc=min(hav,(max(etaav-bc*capac,0.d0)))
               huc=(min(hav,hc)/hsum)*husum
               hvc=(min(hav,hc)/hsum)*hvsum
           else
c               write(*,*) "Setting stuff to zero"
               hc=0.d0
               huc=0.d0
               hvc=0.d0
           endif

c          # set h on coarse grid based on surface, not conservative near shoreline
           valbgc(1,i,j) = hc / capac
           valbgc(2,i,j) = huc / capac
           valbgc(3,i,j) = hvc / capac
         enddo
       enddo

       write(*,*) "Leaving coarsen"
       return
       end


