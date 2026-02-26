c     SUBROUTINE EGLCAO
c     =================
c
c     Copyright 1997 by Peter Blaudeck, Dirk Porezag, Michael Haugk,
c                       Joachim Elsner
c
c *********************************************************************
c
c     PROGRAM CHARACTERISTICS
c     -----------------------
c
c eglcao calculates energy and gradient for dylcao as shown by Seifert. 
c The determination of the occupation numbers has been changed to be
c also valid for metallic systems.
c
c PARAMETERS:
c n3      i  3 * number of atoms
c telec   r  electronic temperature (for Fermi distribution)
c x       r  coordinates (n3)
c scftol  r  convergence criterion for scf (if negative, no scf)
c atomic  l  should there be an atomic energy calculation
c e       r  total energy
c eel     r  electronic energy
c gr      r  gradient (n3)
c miter   i  number of scf-iterations performed
c evector l  prints out the eigenvectors
c lgrad   l  should gradients be calculated?
c
C This program was modified by Yasuhito Ohta      8/1/2008
c *********************************************************************
c
       subroutine eglcao(n3,telec,x,scftol,atomic,e,eel,gr,
     & miter,evector,qmat,lgrad,lexc,nexc,nrem,entropy,timetab,iprint,
     & iMDstep)
c      =================
c
c maxima definitions
c
       implicit REAL*8 (A-H,O-Z)
       include 'maxima.inc'
c
       logical doscf, atomic, evector,constr,lgrad,lexc
       integer izp(NNDIM), ind(NNDIM+1),lmax(MAXTYP)
       integer indj,indk,indj1,indk1,nstart,nexc,nrem
       integer conat(NNDIM+1),cai
       integer numint(MAXTYP,MAXTYP)
       real*8 x(3,*), gr(3,*), hgrad(3,NNDIM), scftol
       real*8 coeff(6,MAXINT,MAXTYP,MAXTYP)
       real*8 xr(2,MAXINT,MAXTYP,MAXTYP)
       real*8 efkt(3,MAXTYP,MAXTYP),cutoff(MAXTYP,MAXTYP)
       real*8 eel, ecoul, e, efermi,convec(3,NNDIM)
       real*8 au(LDIM,LDIM), bu(LDIM,LDIM), nel
       real*8 a(MDIM,MDIM), b(MDIM,MDIM), spro
       real*8 gammamat(NNDIM,NNDIM),derivx(NNDIM,NNDIM)
       real*8 derivy(NNDIM,NNDIM),derivz(NNDIM,NNDIM)
c      real*8 h(MDIM), ev(MDIM), occ(MDIM)
c       modification for sgi-diagonalization
       real*8 h(34*MDIM), ev(MDIM), occ(MDIM)
       real*8 qmulli(MDIM), qmold(NNDIM), qmat(NNDIM)
       real*8 Qdiff(NNDIM),dipol(3),dipabs
       real*8 hamil(MDIM,MDIM), overl(MDIM,MDIM), shift(NNDIM)
       integer OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM
!------------
!     1: S,Hc
!     2: 2e
!     3: diag   
!     4: frc      
!     5:  total eglcao      
!     6:  io  
      real*8 timetab(10)
      
c       real*8 kspr
c      ME constraints like in CHARMM
c       integer lcc,jcc(nndim),atcc1(nndim),atcc2(nndim)
c       real*8 rccval,kccval,rccinc
c       common /cconstr/ lcc,jcc,atcc1,atcc2,rccval,kccval,rccinc
cccccccccccccccc  additional potential
c        real*8 ecorr
c        integer nco,nc(NNDIM,2)
cccccccccccccc CHARMM+ external charges
       real*8 xE(3,4*NNDIM),ZE(4*NNDIM),shiftE(NNDIM)
       integer nE
       character*2 EXT
ccccccccccccccc
       external skhpar,skspar
c     DISPERSION
       logical dispers
       real*8 Edis
       common /disper/ Edis, dispers
c     DISPERSION END
cccccccccccccc CHARMM
       common /extchr/ xE, ZE, nE, nEiter ,icycle, EXT
cccccccccc
cccc additional potential
       common /correct/  r0yz , r0x ,kspr , nco, nc
cccccccccccc
       common /muell/ a,b
c      only for not taking space from stack (FORTRAN problem)
       common /concgr/ convec, constr, conat
c      constrains
       common /lmax/ lmax
c              s : 1, p : 2, u.s.w. for all atom types
       common /izp/ izp, nel, nbeweg
c              atom types, number of electrons, number of movable atoms
       common /spltab/ coeff,xr,efkt,cutoff,numint
c              spline data for repulsiv force

       include 'commonbox.inc'

       common /machine/ racc
c              machine accuracy
       common /mcharge/ qzero(MAXTYP), uhubb(MAXTYP)
c
       common /varhelp/ Qdiff
c
c setup
c
c   parallel  trajectories execution 
       include 'mpiio.inc'     !ninfile, noutfile
       call cpu_time(rtmp0) 

c no gradients in mode 6
       if (atomic) lgrad= .true.  ! MOdified by Yasuhito Ohta: .false. --> .true.
       if (scftol .lt. 0.0) then
         doscf = .false. 
       else
         doscf = .true.
       endif
       nn = n3 / 3
       
c calculation of indices for matrices H and S 
c
       ind(1) = 0
       do j = 1,nn 
         !izpj = izp(j)
         ind(j+1) = ind(j)+lmax(izp(j))**2
       end do
c
c actual dimension of matrix
       ndim = ind(nn+1)      
       
       if (ndim .gt. MDIM) then
         write(istdout,*) ' Stop: eglcao: ndim > ', MDIM
         stop
       endif
       
c setup of charge-independent part of H and S
c
       do j = 1,nn 
         do k = 1,j 
         
           call slkmatrices(j,k,x,au,bu)
           
           do n = 1,ind(k+1)-ind(k) 
             do m = 1,ind(j+1)-ind(j)
               hamil(ind(j)+m,ind(k)+n) = au(m,n)
               hamil(ind(k)+n,ind(j)+m) = au(m,n)       
               overl(ind(j)+m,ind(k)+n) = bu(m,n)
               overl(ind(k)+n,ind(j)+m) = bu(m,n)
             end do
           end do
         end do
       end do
C Yasuhito Ohta: 2007/10/24
!       if(atomic) then
!         open(170,file='SMAT.DAT')
!         write(170,'(a7,I8)') 'ndim = ', ndim
!         write(170,*) 'Smat: overlap matrix elements  '
!         do i = 1,ndim
!           do j = 1,i
!             write(170,*) overl(i,j)
!           enddo
!         enddo
!         close(170)
!       endif
c for external charges : extra shift shiftE 
      if ((EXT.eq.'CH') .or. (EXT.eq.'EF')) then
       call externalshift(nn,x,izp,shiftE)
      endif
c setup for SCF cycle
c
       if (doscf ) then
!        if (atomic) then
!         open(70,file='scfsum',form='formatted',status='unknown') 
!        endif
c         close(70, status='delete')
         almix = 0.2

CYasuhito Ohta: maxiter = 70 -> 200
         maxiter = 70 
         eelold = 1.0d10
       else
          maxiter = 1
       endif
      call cpu_time(rtmp1) 
      timetab(9) =timetab(9)+rtmp1-rtmp0
c  *************************
c  SCF cycle starts here
c  *************************
       do niter = 1,maxiter 
         miter = niter
c save old charges
c
         if (doscf) then
           do i = 1,nn 
             qmold(i) = qmat(i)
           end do
         endif
       
       
c charge-independent part of H and S
c
       call cpu_time(r0)
!$OMP  PARALLEL SHARED(A,B,NCHUNK)
!$OMP+ PRIVATE(id,j,k,n,m,indj,indj1,indk,indk1,indjm,indkn)

CCC!$OMP+ PRIVATE(id,j,k,n,m,indj,indj1,indk,indk1)

       NTHREADS = OMP_GET_NUM_THREADS()
       NCHUNK = (nn+NTHREADS-1)/NTHREADS
!       ID = OMP_GET_THREAD_NUM()
!       write(istdout,*)'NTHREADS, NCHUNK,nn= ', NTHREADS, NCHUNK,nn
!       write(*,*)'NTHREADS, NCHUNK,nn= ', NTHREADS, NCHUNK,nn
       
!$OMP DO SCHEDULE(STATIC, NCHUNK)       
       do j = 1,nn 
         !write(*,*)' THREAD :',ID, '-->',j
         indj  = ind(j)
         indj1 = ind(j+1)
         do k = 1,j 
           indk  = ind(k)
           indk1 = ind(k+1)
           
           do n = 1,indk1-indk 
             do m = 1,indj1-indj 
             indkn = indk+n
               indjm = indj+m
               a(indjm,indkn) = hamil(indjm,indkn)
               b(indjm,indkn) = overl(indjm,indkn)
               !a(indj+m,indk+n) = hamil(indj+m,indk+n)
               !b(indj+m,indk+n) = overl(indj+m,indk+n)
             end do
           end do
         end do
         end do
       

!$OMP END PARALLEL       
       call cpu_time(r1)
       timetab(1) = timetab(1) +r1-r0
       
c zero shift 
c
         do i = 1,nn
           shift(i) = 0.0
         end do
       
c add charge dependent terms (Hubbard, etc. )
c start by defining shift 
c
         if (doscf) then
           call HAMILSHIFT(nn,x,izp,boxsiz,period,uhubb,qmat,qzero,
     &                     miter,gammamat,shift)
ccc add external charges! CHARMM
          if ((EXT.eq.'CH').or. (EXT.eq.'EF'))then 
           do i = 1,nn
           shift(i) = shift(i) - shiftE(i)
c minus, weil Z positiv sind!, positive DeltaQ meinen aber negative Ladung
           end do
          endif 
c
c update hamiltonian matrix 
c
       !write(istdout,*) 'nthreads = ',NTHREADS 
      call cpu_time(r2e0)
!$OMP  PARALLEL SHARED(A,NCHUNK)
!$OMP+ PRIVATE(id,i,j,li,lj)

       NTHREADS = OMP_GET_NUM_THREADS()
       NCHUNK = (nn+NTHREADS-1)/NTHREADS
       ID = OMP_GET_THREAD_NUM()
       !write(istdout,*)'eglcao: NTHREADS',NTHREADS

       
      !write(istdout,*)'eglcao: started TID : ',id
!$OMP DO SCHEDULE(STATIC, NCHUNK)       
      do i = 1,nn 
       do li = 1,lmax(izp(i))**2 
        do j = 1,i
         do lj = 1,lmax(izp(j))**2 
              a(ind(i)+li,ind(j)+lj) = a(ind(i)+li,ind(j)+lj) 
     &        +0.5*overl(ind(i)+li,ind(j)+lj)*(shift(i)+shift(j))
         end do
        end do
       end do
      end do
!$OMP END PARALLEL        
      call cpu_time(r2e1)
      timetab(2) = timetab(2) +r2e1-r2e0

      endif

       
c solution of eigenvalue problem 
c
        call cpu_time(rdia0)
         call ewevge(MDIM,MDIM,ndim,a,b,ev,h,1,-1,ier)
        call cpu_time(rdia1)
        timetab(3) = timetab(3) +rdia1-rdia0
      
         if (ier .ne. 0) then 
           write(istdout,*)' ewevge: ier =',ier,'niter=',niter
           stop
         endif
c Machine accuracy 
c
         dacc = 4*racc
c
c Calculate occupation (occ) and fermi energy (efermi)
c using fermi-distribution
c
         call cpu_time(rtmp0)
         
         call FERMI(nel,telec,ndim,ev,dacc,occ,efermi)
         
c sum of occupied eigenvalues
c
         call cpu_time(tt0)
         eel = 0.0
         do i = 1,ndim 
           if (occ(i) .ge. dacc) eel = eel + occ(i)*ev(i)
           !if (occ(i) .lt. dacc) goto 75
           !eel = eel + occ(i)*ev(i)
         end do


75       continue



         call cpu_time(tt1)
         


c Yasuhito Ohta  ---------------------
c Add entropy contribution,  -Te * entropy

         if( telec .ge. 0.1d0 ) then
           entropy = 0.0d0
           do i=1,ndim
            if( abs(occ(i)) .ge. 1.D-12 ) then
             entropy = entropy + entropy_orb(occ(i))
            endif
           enddo
           entropy = - 3.16683D-6  * 2.0d0 * entropy
c          eel = eel  - Telec * entropy
         endif
c ------------------------------------
c
c Lowest unoccupied level
c
         nstart = i
       
c
c determine mulliken charges , the charge of the whole system
c and the mulliken 
c
         if(doscf) then
          call cpu_time(rmull0)        
          call MULLIKEN(nn,x,izp,qmat,qzero,qmulli,qtot,ndim,dacc,occ,a,
     &                  overl,lmax,ind,dipol,dipabs)
          call cpu_time(rmull1)        
          timetab(6) =  timetab(6) + rmull1-rmull0
         endif
       
c output spektra, occupation numbers, mulliken charges 
c
c        if (atomic) then
c         call outspec(nn,nbeweg,ndim,ind,lmax,izp,ev,occ,efermi,
c     &                qmat,qmulli,qtot,dipol,dipabs)
c        endif 
c         if (evector) then
c           call outeigenvectors(a,ev,occ,ind,nn)
c         endif
        
       
c complete calculation of electronic energy:
c charge-dependent energy contribution
c warning: this will only lead to the right result if convergence
c has been reached
c
         if (doscf) then
           ecoul = 0.0
           eext = 0.0
           do i = 1,nn
            ecoul = ecoul + shift(i)*(qmat(i)+qzero(izp(i)))
            end do
ccccccc CHARMM
           if((EXT.eq.'CH').or. (EXT.eq.'EF'))then
           do i = 1,nn
             eext = eext + shiftE(i)*(qzero(izp(i))-qmat(i))
             end do
           endif
ccccccc

           eel = eel-0.5*ecoul + 0.5*eext
         endif
c remark: eel contains shiftE aready via ev,
c shift also contains -shiftE, i.e. ecoul also
c contains contributions from EXT
       
c print energy, check convergence, call broyden mixing
c
         if (doscf ) then
c           open(70,file='scfsum',form='formatted',status='unknown')
c     &         access='append')
!          if (atomic) then
!           write(70,*) 'iter: ',niter,', E= ',eel
!           !call flush(70)
!           write(71,*) 'iter: ',niter
!           write(71,'(I4,2F8.4,I4,2F8.4,I4,2F8.4,I4,2F8.4)') 
!     c       (i,occ(i),ev(i),i=1,ndim)
!           !call  flush(71)
!c           write(71,'(10F8.4)') (i,occ(i),ev(i),i=1,ndim)
!          endif
c           close(70) 
           if (abs(eel-eelold) .lt. scftol) goto 70
           eelold = eel
c
c   broyden mixing
c
           call broyden(niter,almix,nn,qmold,qmat)
           do i = 1,nn
             qmat(i)= qmold(i)
           end do
         end if
       end do
70     continue
       
       if (doscf) then
         call HAMILSHIFT(nn,x,izp,boxsiz,period,uhubb,qmat,qzero,
     &                   miter,gammamat,shift)
       endif
       
c end of SCF
c
c output spektra, occupation numbers, mulliken charges 
c        iprint =  (modulo(iMDstep,maxrun).eq.0)
c        if (atomic) then
         if (iprint.ge.1) call outspec(nn,nbeweg,ndim,ind,lmax,izp,ev,
     &         occ,efermi,qmat,qmulli,qtot,dipol,dipabs,iprint,iMDstep)
c        endif 
         if (evector) then
           call outeigenvectors(a,ev,occ,ind,nn)
         endif
        
ccc add EXTERNAL CHARGES CHARMM
          if ((EXT.eq.'CH').or. (EXT.eq.'EF'))then
           do i = 1,nn
           shift(i) = shift(i) - shiftE(i)
           enddo
          endif
c minus, weil Z positiv sind!, positive DeltaQ meinen aber negative Ladung  

!       if (atomic) then       
!       
c
c calculation of atomic energies
c 
!       call atomenerg(nn,ndim,ind,izp,shift,
!     &      x,hamil,a,occ,doscf,period,boxsiz,xinvbox,nlat)
!       endif
       
c start with repulsive energy
       erep = 0.0
       
c
c Calculate repulsive energy
c
       call repulsive(nn,izp,period,x,boxsiz,xinvbox,nlat,erep)
       
c total energy 
c
       e = eel + erep


       

c excitation energies via lin. response
c
c       if(lexc) then       
c         call lrespo(nn,ind,ndim,a,ev,overl,x,boxsiz,
c     &                         period,nexc,nrem,e)
c       endif

         call cpu_time(rtmp1)
         timetab(10) =  timetab(10)+ rtmp1-rtmp0
       
c
c start with gradient
       call cpu_time(rf0)        
       do j = 1,nn 
         do i = 1,3 
           gr(i,j) = 0.0
         end do
       end do
c DISERPSION
       if (dispers) then
         call  dispersion_egr(nn,x,hgrad,boxsiz,nlat,period)
        do j=1,nn
         do i=1,3
           gr(i,j) = gr(i,j) + hgrad(i,j)
         end do
        end do
        e = e + Edis
       endif
c END DISPERSION       
c BEGIN EXTERNAL CHARGES CHARMM
      if((EXT.eq.'CH').or. (EXT.eq.'EF'))then
       call externalchgrad(nn,x,izp,hgrad,qmat)
        do j=1,nn
         do i=1,3
           gr(i,j) = gr(i,j) + hgrad(i,j)
         end do
        end do
      endif
c END EXTERNAL CHRARGES
       if(lgrad) then
c
c Calculate gradient of the repulsive energy
c
       call repulsivegrd(nn,nbeweg,izp,period,x,boxsiz,
     &                   xinvbox,nlat,hgrad)
       do j=1,nn 
         do i=1,3 
           gr(i,j) = gr(i,j) + hgrad(i,j)
         end do
       end do
       
c gradient 
c
       call usualgrd(nn,nbeweg,ndim,izp,lmax,ind,period,doscf,
     &              x,ev,a,b,occ,shift,dacc,boxsiz,xinvbox,hgrad)
                    
       do j=1,nn 
        do i=1,3
         gr(i,j) = gr(i,j) + hgrad(i,j)
        end do
       end do
       
c here are the contributions due to gamma if in scf mode
       if (doscf) then
         call GAMMAGRAD(nn,nbeweg,x,izp,boxsiz,period,uhubb,
     &                  qmat,qzero,derivx,derivy,derivz,hgrad)
           do j=1,nn 
            do i=1,3 
              gr(i,j) = gr(i,j) + hgrad(i,j)
            end do
           end do
       endif
       
c
c constraints?
c
       if(constr) then
        do i=1,conat(1)
         cai=conat(i+1)
         spro=convec(1,i)*gr(1,cai) 
         spro=spro+convec(2,i)*gr(2,cai) 
         spro=spro+convec(3,i)*gr(3,cai) 
c         gr(1,cai)=gr(1,cai)-spro*convec(1,i)
c         gr(2,cai)=gr(2,cai)-spro*convec(2,i)
c         gr(3,cai)=gr(3,cai)-spro*convec(3,i)
c other constraint! fix x,y,z individually
         gr(1,cai)=gr(1,cai)*convec(1,i)
         gr(2,cai)=gr(2,cai)*convec(2,i)
         gr(3,cai)=gr(3,cai)*convec(3,i)
        end do
       endif
c
c End of gradient calculation , only done when not in mode 6
c
       call cpu_time(rf1)
       timetab(4)= timetab(4) + rf1-rf0
       endif
!       write(200,111) r1-r0,r2e1-r2e0,rdia1-rdia0,rf1-rf0
111    format(' Hc+S, 2e, diag  : ', 100(' ',F16.6))
       
       end
c**********************************************************************
      real*8 function entropy_orb(occ_num) 
c                                             Aug 6, 2007
c                                            Yasuhito Ohta
c**********************************************************************
      implicit none 
      real*8 occ_num, a, b, const1, const2

      const1 = 0.5d0 * occ_num
      const2 = 1.0d0 - 0.5d0 * occ_num 

  
      if( const1 .ge. 1.0D-12 ) then
        a =  const1 *log( const1 ) 
      else
        a = 0.0d0
      endif

      if( const2 .ge. 1.0D-12 ) then
        b = const2*log( const2 ) 
      else
        b = 0.0d0
      endif

      entropy_orb = a + b

      return
      end
