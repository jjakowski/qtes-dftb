C Output from Public domain Ratfor, version 1.0
      subroutine mulliken(nn,x,izp,qmat,qzero,qmulli,qtot,ndim,dacc,occ,
     * a,overl,lmax,ind,dipol,dipabs)
C
C maxima definitions
C
      implicit real*8 (a-h,o-z)
      include 'maxima.inc'
C
      integer nn,ndim,lmax(maxtyp),izp(nndim),ind(nndim+1)
      real*8 dacc,occ(mdim),a(mdim,mdim),overl(mdim,mdim)
      real*8 x(3,*),qzero(maxtyp),qtot,qmulli(mdim)
      real*8 qmat(nndim),dipol(3),dipabs,afoo(ndim,ndim)
      integer i,j,lj,m,n,izpj,adim
      real*8 sum,qhelp,conv

      call cpu_time(t0)
      do i = 1,ndim 
        if(occ(i) .lt. dacc)then
          goto 23001
        endif
      enddo
23001 continue

      adim = i-1

      call cpu_time(t1)
      call dsymm('L','U',ndim,adim,1.0d0,overl,mdim,a,mdim,
     &            0.0d0,afoo,ndim)
      call cpu_time(t2)

      do n = 1,ndim 
        qmulli(n) = 0.0
        do i = 1,adim 
          qmulli(n) = qmulli(n) + occ(i)*afoo(n,i)*a(n,i)
        enddo
      enddo

      do j = 1,nn 
        qtot = 0.0
        do lj = 1,lmax(izp(j)) 
          jofn = ind(j)+(lj-1)**2
          qhelp = 0.0

          do mj = 1,2*lj-1 
           qhelp = qhelp + qmulli(jofn+mj)
          enddo

          qtot = qtot + qhelp
        enddo
        qmat(j) = qtot
      enddo

C
C calculation of total charge 
C
      qtot = 0.0
      do j = 1,nn 
        qtot = qtot+qmat(j)
      enddo

C
C calculation of dipol moment
C
      do i = 1,3 
        dipol(i) = 0.d0
        do j = 1,nn 
         !izpj = izp(j)
         !qhelp = qzero(izpj) - qmat(j)
         qhelp = qzero(izp(j)) - qmat(j)
         dipol(i) = dipol(i) + qhelp*x(i,j)
        enddo
C conversion to debye
        dipol(i) = dipol(i)*2.541765d0
      enddo

C
C norm of dipol moment
C
      dipabs = 0.d0
      do i = 1,3 
        dipabs = dipabs + dipol(i)**2
      enddo

      dipabs = sqrt(dipabs)

      call cpu_time(t3)
!      write(istdout,111)t1-t0,t2-t1,t3-t2
111   format('bef, dsymm, aft:',3(' ', F14.8))      
      end
