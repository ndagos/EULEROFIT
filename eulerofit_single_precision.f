C --------------------------------------------
c
c    EULEROFIT
c
c  Programma per il calcolo del polo euleriano
c    da un file velocita' GMT
c
c   Legge un file input di velocita' (velocity.in) ed usa
c    i siti geodetici elencati nel file sites.fix per
c    calcolare ai minimi quadrati il polo di rotazione
c    euleriano
c
c  Output files:
c    
c    velmod.out
c       velocita' e sigma previsti dal polo di rotazione calcolato per i siti in fixsta
c
c    velcalc.out
c       velocita' calcolate dal polo di rotazione per tutti i siti in input. 
c       Errori calcolati propagando l'errore del polo di rotazione (solo se imode=2)
c   
c    rescalc.out
c       residui velocita' di tutte le stazioni.
c       Gli errori calcolati includono l'errore iniziale piu' l'errore
c       associato al polo di rotazione (se imode=2 e iprop=1)
c
c  Nicola D'Agostino apr. 2003  
c
c  10/08/2004:
c     aggiunta conversione latitudine geodetica->geocentrica in input
c
c  12/10/2004:
c     data resolution matrix  (Menke)
c
c  01/09/2006:
c     inversione similtanea di +1 placche 
c
c  22/04/2011
c     subroutines LAPACK per calcolo eigenvalues/eigenvectors     
c
      program eulerofit
      implicit real*4 (a-h,o-z)
      parameter(nmax=8000,nbig=80000,np=50)
      real*4 cgeo(nmax),vel(nmax),covvel(nmax,nmax)
      real*4 velmod(nmax),covmod(nmax,nmax),covvelinv(nmax,nmax)
      real*4 velcalc(nmax),covcalc(nmax,nmax)
      real*4 rescalc(nmax),addcalc(nmax)
      real*4 cgeotot(nmax),veltot(nmax),covveltot(nmax,nmax)
c
      real*4 rmat(2,3),amat(3,3),tmp1(2,3),tmat(nmax,np)
      real*4 veltmp(nmax),tmpmat(nmax,np),covtmp(nmax,nmax)
      real*4 ctmp1(nmax,nmax),ctmp2(nmax,nmax),rotmp(np,np)
      real*4 wp(np),covwp(np,np)
      real*4 wp1(3),covwp1(3,3),rotmat(3,3)
      real*4 wprot(3),covwprot(3,3),cov2d(2,2),eigval(2),eigvect(2,2)
      real*4 wres(nmax),res(nmax),psum(nmax)
      real*4 work(nbig)
      character*4 sta,fixsta(nmax),stafix(nmax),statot(nmax)
      character*2 pl(nmax),plfix(nmax),plid(20),plate
      integer ipv(nmax),imp(nmax),ipl(nmax)
      character*35 filevel,filefix,filepole
c
      write(*,*) " ============================================="
      write(*,*) "         EULEROFIT  "
c
      data bigr/6371008.E0/
      pi=acos(-1.E0)
      pi2=pi*pi
      rad2dg=180.E0/pi
      dg2rad=pi/180.E0
      ff=1.E0/298.2572235630D0
      e2=2*ff-ff**2
c
c      write(*,*) 'e2 =',e2
c -----------------------------------------------
c   Read drive file
c
c  drive_mode 
c     0 = pole estimate (needs fixed points in filefix)
c     1 = read pole and rotate velocities (needs an eulerian pole in filepole)
c           format: lat,lon,rot                         (deg/Myr)
c     2 = read pole and rotate velocities (needs an eulerian pole in filepole)
c           format: rotX, rotY, rotZ                    (10^-3 rad/Myr)
c                   covXX,covXY,covXZ,covYY,covYZ,covZZ (10^-6 rad/Myr)
c  filevel
c  filepole
c  filefix (only if drive_mode = 0)
c
      open(11,file="eulerofit.drv")
      read(11,*) imode, iprop
      read(11,'(35a)') filevel
      read(11,'(a)') filepole
      read(11,'(a)') filefix
      close(11)
c ------------------------------------------------
c   Read fixed sites 
c
      if(imode.eq.0) then
      open(10,file=filefix)
      nfx=1
   10 read(10,*,end=20) fixsta(nfx), plfix(nfx)
      nfx=nfx+1
      goto 10
   20 nfx=nfx-1
      write(*,*) 'Read fixed sites =',nfx
      close(10)
      endif
c ------------------------------------------------
c  Read input pole file
c
      if(imode.eq.2) then
       npl=1
       plid(1)="RF"
       open(11,file=filepole)
       read(11,*) eulat,eulon,rot
       rot=rot*dg2rad
       wp(1)=rot*cos(eulat*dg2rad)*cos(eulon*dg2rad)*1.E-6
       wp(2)=rot*cos(eulat*dg2rad)*sin(eulon*dg2rad)*1.E-6
       wp(3)=rot*sin(eulat*dg2rad)*1.E-6
        do i=1,3
           covwp(i,i)=0.0
        enddo
      elseif(imode.eq.1) then
       npl=1
       plid(1)="RF"
       open(11,file=filepole)
       read(11,*) wp(1),wp(2),wp(3)
       read(11,*) covwp(1,1),covwp(2,1),covwp(3,1),
     & covwp(2,2),covwp(2,3),covwp(3,3)
       covwp(3,2)=covwp(2,3)
       covwp(1,2)=covwp(2,1)
c
c --- riscala la rotazione da 10^-3 rads/Myr a rad/yr ---
c
       covwp(1,3)=covwp(3,1)
        do i=1,3
          wp(i)=wp(i)*1.E-9
c
c --- e la covarianza da rad**2/Myr**2 a rad**2/yr**2
c
         do j=1,3
          covwp(i,j)=covwp(i,j)*1.E-18
         enddo
        enddo
       close(11)
       endif
c ------------------------------------------------------------
c   Input file in formato GMT psvelo -Se  (velocita' in mm/yr)
c     [Lon, Lat, Vel_E, Vel_N, Sig_E, Sig_N, Corr, Sta]
c
      open(20,file=filevel) 
      i=1
      itot=1
   30 read(20,*,end=40) c1,c2,v1,v2,se,sn,corr,sta 
c      write(*,'(2F12.4,5F10.2,A6)') c1,c2,v1,v2,se,sn,corr,sta
c         if(se.eq.0.0) se=0.5
c         if(sn.eq.0.0) sn=0.5
         cgeotot(itot)=c1
         cgeotot(itot+1)=c2
         veltot(itot)=v1
         veltot(itot+1)=v2
         covveltot(itot,itot)=se**2
         covveltot(itot+1,itot+1)=sn**2
         covveltot(itot,itot+1)=corr*sn*se
         covveltot(itot+1,itot)=covveltot(itot,itot+1)
         statot(itot)=sta
         statot(itot+1)=sta
         itot=itot+2
      if(imode.eq.0) then
      do ii=1,nfx
       if(sta.eq.fixsta(ii)) then
         kk=(i+1)/2
         pl(kk)=plfix(ii)
         stafix(i)=sta
         stafix(i+1)=sta
         cgeo(i)=c1
         cgeo(i+1)=c2
         vel(i)=v1
         vel(i+1)=v2
         covvel(i,i)=se**2
         covvel(i+1,i+1)=sn**2
         covvel(i,i+1)=corr*sn*se
         covvel(i+1,i)=covvel(i,i+1)
c         write(*,'(2F12.4,5F10.2,A6,A4)')  
c     & cgeo(i),cgeo(i+1),vel(i),vel(i+1),se,sn,corr,sta,plfix(kk) 
         i=i+2
        endif
      enddo
      endif
      goto 30
   40 write(*,*) 'End of filein '
       nv=(i-1)/2
       ntot=(itot-1)/2
      write(*,*) 'Fixed sites found in velocity file = ',nv
      write(*,*) 'Total sites found in velocity file = ',ntot
c
      close(20)
c
c ------ conta quante placche ed associa plid -----------
c
      if(imode.eq.0) then
       npl=1
       plid(1)=pl(1)
       do i=1,nv
        iflg=0
        do j=1,npl
         if(pl(i).eq.plid(j)) then
           iflg=1
           ipl(i)=j
         endif
        enddo
         if(iflg.ne.1) then
          npl=npl+1
          plid(npl)=pl(i)
          ipl(i)=j
         endif
       enddo
       write(*,*) 'found ', npl, ' plates'
      endif
c -------------------------------------------------------------------
c   Forma la matrice per l'inversione minimi quadrati in XYZ
c
      if(imode.eq.0) then
      do i=1,nv
      kk=ipl(i)*3-2
      ii=i*2-1
      iii=i*3-2
      alon=cgeo(ii)
      alat=cgeo(ii+1)
c
      cx=bigr*cos(alat*dg2rad)*cos(alon*dg2rad)
      cy=bigr*cos(alat*dg2rad)*sin(alon*dg2rad)
      cz=bigr*sin(alat*dg2rad)
c
      dist=sqrt(cx**2+cy**2+cz**2)
c
      amat(1,2)=cz*1.E3
      amat(1,3)=-cy*1.E3
      amat(2,1)=-cz*1.E3
      amat(2,3)=cx*1.E3
      amat(3,1)=cy*1.E3
      amat(3,2)=-cx*1.E3
c ----------------------------------------------------------------
c   Prepara la matrice di rotazione R di trasformazione Vxyz => Ven
c    
       rmat(1,1)=-sin(dg2rad*alon)
       rmat(1,2)=cos(dg2rad*alon)
       rmat(1,3)=0.D0
       rmat(2,1)=-sin(dg2rad*alat)*cos(dg2rad*alon)
       rmat(2,2)=-sin(dg2rad*alat)*sin(dg2rad*alon)
       rmat(2,3)=cos(dg2rad*alat)
c -------------------------------------------------------------------
c  Moltiplica le due matrici per formare la design matrix finale Tmat
c      tmat = rmat * amat
c
       call sgemm ('N','N',2,3,3,1.0,rmat,2,amat,3,0.0,tmp1,2)
         tmat(ii,kk)    =tmp1(1,1)
         tmat(ii+1,kk)  =tmp1(2,1)
         tmat(ii,kk+1)  =tmp1(1,2)
         tmat(ii+1,kk+1)=tmp1(2,2)
         tmat(ii,kk+2)  =tmp1(1,3)
         tmat(ii+1,kk+2)=tmp1(2,3)
      enddo
c ----------------------------------------------------------------
c  Inversione ai minimi quadrati del polo di rotazione euleriano
c 
c   tmat*wp = vel                                   (Model equation)
c   P = covvel^-1                                   (Weigth matrix)
c   wp = [(tmat * P * tmat')^-1 ] [tmat * P * vel]  (Least squares solution)
c     
      m1=nv*2
      do i=1,m1
       do j=1,3*npl
        tmpmat(i,j)=tmat(i,j)
       enddo
       do k=1,m1
        covtmp(i,k)=covvel(i,k)
       enddo
        veltmp(i)=vel(i)
      enddo
c     
      call sggglm(m1,npl*3,m1,tmpmat,nmax,covtmp,nmax,veltmp,wp,wres,
     & work,nbig,info)
      write(*,*) 'sggglm info = ',info
      write(*,*) 'sggglm Optimal dimension for work (nbig) = ',work(1)

c --- calcola velocita' modello ----
c
      call sgemm('N','N',m1,1,3*npl,1.0,tmat,nmax,wp,np,0.0,velmod,
     & nmax)
       do i=1,nv
        ii=i*2-1
        write(*,'(2F12.2,A6)') velmod(ii),velmod(ii+1),stafix(ii)
c        write(*,'(4F12.2,A6)') vel(ii),vel(ii+1),
c     &                    sqrt(covvel(ii,ii)),sqrt(covvel(ii+1,ii+1)),stafix(ii)
c        write(*,*) '--------------------------------------'
       enddo
c
c --- calcola residui --------
c
      do ii=1,m1
       res(ii)=vel(ii)-velmod(ii)
      enddo
c
c --- calcola inversa matrice covarianza osservazioni P e chi2/dof ----
c
      do i=1,m1
        do j=1,m1
          covtmp(i,j)=covvel(i,j)
        enddo
      enddo
c
      call sgetrf(m1,m1,covtmp,nmax,ipv,info)
       write(*,*) 'sgetrf info = ', info
      call sgetri(m1,covtmp,nmax,ipv,work,-1,info)
       write(*,*) 'sgetri info = ', info
       write(*,*) 'sgetri optimal dimension for work=',work(1)
      call sgetri(m1,covtmp,nmax,ipv,work,nbig,info)
       write(*,*) 'sgetri info = ', info
       write(*,*) 'sgetri optimal dimension for work=',work(1)
c
      do i=1,m1
       do j=1,m1
        covvelinv(i,j)=covtmp(i,j)
       enddo
      enddo
c
c ---- calcola res^T*covtmp*res -----
c
      call sgemm('T','N',1,m1,m1,1.0,res,nmax,covtmp,nmax,0.0,
     & ctmp1,nmax)
      call sgemm('N','N',1,1,m1,1.0,ctmp1,nmax,res,nmax,0.0,
     & ctmp2,nmax)
       chi2=ctmp2(1,1)
       idof=m1-npl*3
       chi2dof=chi2/float(idof)
       write(*,*) '----------------------------------'
       write(*,'("  CHI2      = ",F10.3)') chi2
       write(*,'("  dof       = ",I10)') idof
       write(*,'("  CHI^2/dof = ",F10.3)') chi2dof
       write(*,*) '----------------------------------'
c
c ---- calcola matrice covarianza parametri -----
c
      do i=1,nmax
       do j=1,3*npl
        tmpmat(i,j)=tmat(i,j)
       enddo
c       do k=1,nmax
c        ctmp1(i,k)=0.D0
c        ctmp2(i,k)=0.D0
c       enddo
      enddo
c
      call sgemm('T','N',3*npl,m1,m1,1.0,tmpmat,nmax,covtmp,nmax,
     & 0.0,ctmp1,nmax)
      call sgemm('N','N',3*npl,3*npl,m1,1.0,ctmp1,nmax,tmpmat,nmax,
     & 0.0,covwp,np)
c
c --- inverte tmat^T*P*tmat  ----
c
c      write(*,*) 'inverte matrice tmat^T*P*tmat'
      call sgetrf(3*npl,3*npl,covwp,np,ipv,info)
      write(*,*) 'sgetrf info = ', info
      call sgetri(3*npl,covwp,np,ipv,work,nbig,info)
      write(*,*) 'sgetri info = ', info
      write(*,*) 'sgetri optimal dimension for work=',work(1)
c
c --- calcola generalizzata inversa G^g [tmat^T*P*tmat]^-1*tmat^T*P -------
c
      call sgemm('N','T',3*npl,m1,3*npl,1.0,covwp,np,tmat,nmax,0.0,
     & ctmp1,nmax)
      call sgemm('N','N',3*npl,m1,m1,1.0,ctmp1,nmax,covvelinv,nmax,
     & 0.0,ctmp2,nmax)
c
c --- Data resolution matrix -------
c
      call sgemm('N','N',m1,m1,3*npl,1.0,tmat,nmax,ctmp2,nmax,0.0,
     & ctmp1,nmax)
c
c --- Calcola data importance per ogni stazione e ordina il vettore psum ------
c
      sum=0.0
        do i=1,nv
         ii=2*i-1
         psum(i)=ctmp1(ii,ii)+ctmp1(ii+1,ii+1)
         sum=sum+ctmp1(ii,ii)+ctmp1(ii+1,ii+1)
         imp(i)=i
        enddo
c
      do j=2,nv
       a=psum(j)
       n=imp(j)
        do i=j-1,1,-1
          if(psum(i).ge.a) goto 55
          psum(i+1)=psum(i)
          imp(i+1)=imp(i)
        enddo
          i=0
   55 psum(i+1)=a
      imp(i+1)=n
      enddo
c
      write(*,*) 'Data Importance Table'
      write(*,56)
   56 format('      Site      East     North      Tot       ',
     & ' %        %cum')
      pcentcum=0.0
       do i=1,nv
        ii=2*imp(i)-1
         ppsum=ctmp1(ii,ii)+ctmp1(ii+1,ii+1)
         pcent=100.0*ppsum/sum
         pcentcum=pcentcum+pcent
         write(*,60) stafix(ii),ctmp1(ii,ii),ctmp1(ii+1,ii+1),
     &               ppsum,pcent,pcentcum,pl(imp(i))
   60    format(A10,2F10.2,F10.2,2F10.1,A6,'  Data Importance Table')
         enddo
        write(*,*) 'sum = ', sum
c
c ---- scala la matrice di covarianza dei parametri ------
c
c        do i=1,3*npl
c         do j=1,3*npl
c          covwp(i,j)=chi2dof*covwp(i,j)
c         enddo
c        enddo
c
c --- calcola matrice covarianza velocita' modello tmat * covwp * tmat^T --------
c
      call sgemm('N','N',m1,3*npl,3*npl,1.0,tmat,nmax,covwp,np,
     & 0.0,ctmp1,nmax)
      call sgemm('N','T',m1,m1,3*npl,1.0,ctmp1,nmax,tmat,nmax,0.0,
     & covmod,nmax)
c
c --- tavola dei residui --------
c
      aae=0.0
      bbe=0.0
      aan=0.0
      bbn=0.0
      aat=0.0
      bbt=0.0
      rmse=0.0
      rmsn=0.0
c
      write(*,*) '----------------------------------'
      write(*,"('    %      DCHI2 (2dof)')")
      write(*,"('   68.3     2.30 ')")
      write(*,"('   90.0     4.61 ')")
      write(*,"('   95.4     6.17 ')")
      write(*,"('   99.0     9.21 ')")
      write(*,"('   99.99   18.40 ')")
      write(*,*) '-----------------------------------------------------'
c
      do j=1,npl      
       write(*,*) 'PLATE ',plid(j)
      write(*,*)  'Residuals:'
      write(*,50)
   50 format('      E(mm/yr)    N(mm/yr)      sigE     ',
     & '   sigN      Fix_Sta       dchi2')

       npp=0
       aae=0.
       aan=0.
       bbe=0.
       bbn=0.
       aat=0.
       bbt=0.
      do i=1,nv
       if(pl(i).eq.plid(j)) then
       ii=i*2-1
       cov2d(1,1)=covvel(ii,ii)+covmod(ii,ii)
       cov2d(2,1)=covvel(ii+1,ii)+covmod(ii+1,ii)
       cov2d(2,2)=covvel(ii+1,ii+1)+covmod(ii+1,ii+1)
       cov2d(1,2)=cov2d(2,1)
       call sgetrf(2,2,cov2d,2,ipv,info)
c         write(*,*) 'info = ',info
       call sgetri(2,cov2d,2,ipv,work,nbig,info)
c         write(*,*) 'info = ',info
       fchi2=cov2d(1,1)*res(ii)**2+
     &       cov2d(2,2)*res(ii+1)**2+
     &     2*cov2d(2,1)*res(ii)*res(ii+1)

       aae=aae+(res(ii)**2/covvel(ii,ii))
       bbe=bbe+(1.0/covvel(ii,ii))
c
       aan=aan+(res(ii+1)**2/covvel(ii+1,ii+1))
       bbn=bbn+(1.0/covvel(ii+1,ii+1))
c
       aat=aat+(res(ii)**2/covvel(ii,ii))+
     &         (res(ii+1)**2/covvel(ii+1,ii+1))
       bbt=bbt+(1.0/covvel(ii,ii))+(1.0/covvel(ii+1,ii+1))
       rmse=rmse+res(ii)**2
       rmsn=rmsn+res(ii+1)**2
c
c      out=1.D0
       out1=sqrt(covvel(ii,ii))
       out2=sqrt(covvel(ii+1,ii+1))
       write(*,"(4F12.2,A12,2X,F10.2,A6)")
     &      res(ii),res(ii+1),out1,out2,stafix(ii),fchi2,pl(i)
       npp=npp+1
       endif
      enddo
      wrmse=sqrt(aae/bbe)
      wrmsn=sqrt(aan/bbn)
      wrmst=sqrt(aat/bbt)
      rmse=sqrt(rmse/float(npp))
      rmsn=sqrt(rmsn/float(npp))
      cdf=aat/float(2*npp-3)
      write(*,*) 'WRMS (mm/yr): East    North     Tot'
      write(*,'(7X,3F10.3)') wrmse , wrmsn, wrmst
      write(*,*) ' RMS (mm/yr): East    North'
      write(*,'(7X,2F10.3)') rmse , rmsn
c
      write(*,'("  CHI2      = ",F10.3)') aat
      write(*,'("  dof       = ",I10)') 2*npp-3
      write(*,'("  CHI^2/dof = ",F10.3)') cdf
      write(*,*) '==================================================='
      enddo
c
      write(*,'("Multi-Plate Statistics")')
      write(*,'("  CHI2      = ",F10.3)') chi2
      write(*,'("  dof       = ",I10)') idof
      write(*,'("  CHI^2/dof = ",F10.3)') chi2dof
      write(*,*) '----------------------------------'
c
c ================= fine ciclo su imod=0 =======================================
c
      endif
c
c --- ciclo su npl per scrivere residui/calcolati/poli_rotazione
c
      open(40,file='rescalc.out')
      open(45,file='addcalc.out')
      open(50,file='velcalc.out')
      open(60,file='enupole.out')
      open(70,file='xyzpole.out')
c -------------------------------
      m1=ntot*2
      do k=1,npl
       write(*,*) 'Scrive residui/poli per ', plid(k)
       write(*,*) '--------------------------------------'
C
       do ii=1,3                
         wp1(ii)=wp(k*3-3+ii)
         plate=plid(k)
        do jj=1,3
         covwp1(ii,jj)=covwp(k*3-3+ii,k*3-3+jj)
        enddo
       enddo
c
c      do ii=1,nmax
c       do jj=1,nmax
c        tmat(ii,jj)=0.D0
c       enddo
c      enddo
c -------------------------------------------------------------------
c   Forma la matrice per il calcolo velocita teoriche
c
      do i=1,ntot
      ii=i*2-1
      iii=i*3-2
      alon=cgeotot(ii)
      alat=cgeotot(ii+1)
c
      cx=bigr*cos(alat*dg2rad)*cos(alon*dg2rad)
      cy=bigr*cos(alat*dg2rad)*sin(alon*dg2rad)
      cz=bigr*sin(alat*dg2rad)
c
      dist=sqrt(cx**2+cy**2+cz**2)
c
      amat(1,2)=cz*1.E3
      amat(1,3)=-cy*1.E3
      amat(2,1)=-cz*1.E3
      amat(2,3)=cx*1.E3
      amat(3,1)=cy*1.E3
      amat(3,2)=-cx*1.E3
c ----------------------------------------------------------------
c   Prepara la matrice di rotazione R di trasformazione Vxyz => Ven
c
       rmat(1,1)=-sin(dg2rad*alon)
       rmat(1,2)=cos(dg2rad*alon)
       rmat(1,3)=0.0
       rmat(2,1)=-sin(dg2rad*alat)*cos(dg2rad*alon)
       rmat(2,2)=-sin(dg2rad*alat)*sin(dg2rad*alon)
       rmat(2,3)=cos(dg2rad*alat)
c -------------------------------------------------------------------
c  Moltiplica le due matrici per formare la design matrix finale Tmat
c      tmat = rmat * amat
c
       call sgemm ('N','N',2,3,3,1.0,rmat,2,amat,3,0.0,tmp1,2)
         tmat(ii,  1)=tmp1(1,1)
         tmat(ii+1,1)=tmp1(2,1)
         tmat(ii,  2)=tmp1(1,2)
         tmat(ii+1,2)=tmp1(2,2)
         tmat(ii,  3)=tmp1(1,3)
         tmat(ii+1,3)=tmp1(2,3)
      enddo
c
      call sgemm('N','N',m1,1,3,1.0,tmat,nmax,wp1,3,0.0,velcalc,nmax)
      do i=1,m1
       rescalc(i)=veltot(i)-velcalc(i)
       addcalc(i)=veltot(i)+velcalc(i)
      enddo
c
c --- calcola matrice covarianza velocita' modello ---------
c
      do i=1,nmax
       do j=1,nmax
        ctmp1(i,j)=0.0
        ctmp2(i,j)=0.0
        covcalc(i,j)=0.0
       enddo
      enddo
      call sgemm('N','N',m1,3,3,1.0,tmat,nmax,covwp1,3,0.0,ctmp1,nmax)
      call sgemm('N','T',m1,m1,3,1.0,ctmp1,nmax,tmat,nmax,0.0,
     & covcalc,nmax)
c
c ------ scrive velocita' in output ---------
c
      do i=1,ntot
       ii=i*2-1
       se=sqrt(covcalc(ii,ii))
       sn=sqrt(covcalc(ii+1,ii+1))
       corr=covcalc(ii+1,ii)/(sn*se)
       write(50,'(2F10.4,4F9.2,F8.2,A8,A4)') cgeotot(ii),cgeotot(ii+1),
     & velcalc(ii),velcalc(ii+1),se,sn,corr,statot(ii),plate
c
      if(iprop.eq.1)then
        stote=covcalc(ii,ii)+covveltot(ii,ii)
        stotn=covcalc(ii+1,ii+1)+covveltot(ii+1,ii+1)
        stoten=covcalc(ii,ii+1)+covveltot(ii,ii+1)
        se=sqrt(stote)
        sn=sqrt(stotn)
        corr=stoten/(se*sn)
      else
        stote=covveltot(ii,ii)
        stotn=covveltot(ii+1,ii+1)
        stoten=covveltot(ii,ii+1)
        se=sqrt(stote)
        sn=sqrt(stotn)
        corr=stoten/(se*sn)
      endif
      write(40,'(2F10.4,4F9.2,F8.2,A8,A4)') cgeotot(ii),cgeotot(ii+1),
     & rescalc(ii),rescalc(ii+1),se,sn,corr,statot(ii),plate
      write(45,'(2F10.4,4F9.2,F8.2,A8,A4)') cgeotot(ii),cgeotot(ii+1),
     & addcalc(ii),addcalc(ii+1),se,sn,corr,statot(ii),plate
      enddo
c
      eulon=atan2(wp1(2),wp1(1))*rad2dg
      aa=sqrt(wp1(1)**2+wp1(2)**2)
      eulat=atan2(wp1(3),aa)*rad2dg
      rot=sqrt(wp1(1)**2+wp1(2)**2+wp1(3)**2)
c
c --- converte rot da rad/yr a deg/my
c
      rot=rot*rad2dg*1.E6
c
      if(eulon.gt.360.0) eulon=eulon-360.0
c
c ----- ruota le componenti del vettore di rotazione wp..... in riferimento ENU  ----
c
      rotmat(1,1)=-sin(dg2rad*eulon)
      rotmat(1,2)=cos(dg2rad*eulon)
      rotmat(1,3)=0.0
c
      rotmat(2,1)=-sin(dg2rad*eulat)*cos(dg2rad*eulon)
      rotmat(2,2)=-sin(dg2rad*eulat)*sin(dg2rad*eulon)
      rotmat(2,3)=cos(dg2rad*eulat)
c
      rotmat(3,1)=-cos(dg2rad*eulat)*cos(dg2rad*eulon)
      rotmat(3,2)=-cos(dg2rad*eulat)*sin(dg2rad*eulon)
      rotmat(3,3)=-sin(dg2rad*eulat)
      call sgemm('N','N',3,1,3,1.0,rotmat,3,wp1,3,0.0,wprot,3)
c
c ---- .... e la sua matrice di covarianza in coordinate ENU -----
c
      call sgemm ('N','N',3,3,3,1.0,rotmat,3,covwp1,3,0.0,rotmp,3)
      call sgemm ('N','T',3,3,3,1.0,rotmp,3,rotmat,3,0.0,covwprot,3)
c
      do i=1,2
       do j=1,2
         cov2d(i,j)=covwprot(i,j)
       enddo
      enddo
c
c ---- Calcolo eigenvalues/vectors with Numerical Recipes --
c      (deprecated: not compiling well with gfortran 21 apr 2011)
c
c      call jacobi(cov2d,2,2,eigval,eigvect,nrot)
c      call eigsrt(eigval,eigvect,2,2)
c      az=datan2(eigvect(1,1),eigvect(2,1))*rad2dg
c      aa=dabs(wprot(3))
c      smaj=(sqrt(eigval(1))/aa)*rad2dg
c      smin=(sqrt(eigval(2))/aa)*rad2dg
c
c ---- Use LAPACK subroutines for eigenvalues/vectors 
c
      call SSYEV('V','U',2,cov2d,2,eigval,work,nmax,info)
      write(*,*) 'ssyef info =',info
      az=atan2(cov2d(1,2),cov2d(2,2))*rad2dg
      aa=abs(wprot(3))
      smaj=(sqrt(eigval(2))/aa)*rad2dg*sqrt(2.0)
      smin=(sqrt(eigval(1))/aa)*rad2dg*sqrt(2.0)
c
      if(az.gt.90.0) az=az-180.0
      if(az.lt.-90.0) az=az+180.0
c
      rote=wprot(1)*rad2dg*1.E6
      srote=sqrt(covwprot(1,1))*rad2dg*1.E6
      rotn=wprot(2)*rad2dg*1.E6
      srotn=sqrt(covwprot(2,2))*rad2dg*1.E6
      rotu=wprot(3)*rad2dg*1.E6
      srotu=sqrt(covwprot(3,3))*rad2dg*1.E6
c
        if(eulon.ge.180.0) eulon=eulon-360.0
        if(eulon.le.-180.0) eulon=eulon+360.0
        write(60,165) eulat,eulon,rot,smaj,smin,az,srotu,plate
        eulon=eulon-180.0
        if(eulon.ge.180.0) eulon=eulon-360.0
        if(eulon.le.180.0) eulon=eulon+360.0
        eulat=-eulat
        rot=-rot
        az=-az
        write(60,165) eulat,eulon,rot,smaj,smin,az,srotu,plate
  165 format(2F10.3,F8.3,3F8.1,F8.3,2X,A2,4X,'POLE')
c
      write(70,170) wp1(1)*1.E9,wp1(2)*1.D9,wp1(3)*1.E9,plate
  170 format(3F14.6,2X,A2,4X,'ROT')
      write(70,180) covwp1(1,1)*1.E18,covwp1(1,2)*1.E18,
     &              covwp1(3,1)*1.E18,covwp1(2,2)*1.E18,
     &              covwp1(2,3)*1.E18,covwp1(3,3)*1.E18,
     &              plate
  180 format(6F14.6,2X,A2,4X,'COVAR')
      enddo
c ===================================================================
c scrive i poli di rotazione relativi
c
      if(npl.gt.1) then
       write(*,*) 'scrive poli rotazione relativi'
c
      do i=1,npl
       do j=1,npl
c
      do ii=1,3
         wp1(ii)=0.0
       do  jj=1,3
         covwp1(ii,jj)=0.0
       enddo
      enddo
c
      do ii=1,nmax
       do jj=1,nmax
         ctmp1(ii,jj)=0.0
         ctmp2(ii,jj)=0.0
         covtmp(ii,jj)=0.0
       enddo
      enddo
c
       if(i.ne.j) then
        ii=(i-1)*3
        jj=(j-1)*3
c        wp1(1)=wp(ii+1)-wp(jj+1)
c        wp1(2)=wp(ii+2)-wp(jj+2)
c        wp1(3)=wp(ii+3)-wp(jj+3)
          ctmp1(1,ii+1)=1.0
          ctmp1(2,ii+2)=1.0
          ctmp1(3,ii+3)=1.0
          ctmp1(1,jj+1)=-1.0
          ctmp1(2,jj+2)=-1.0
          ctmp1(3,jj+3)=-1.0
      call sgemm('N','N',3,1,3*npl,1.0,ctmp1,nmax,
     &           wp,np,0.0,wp1,3)
      call sgemm('N','N',3,3*npl,3*npl,1.0,ctmp1,nmax,
     &           covwp,np  ,0.0,ctmp2 ,nmax)
      call sgemm('N','T',3,3    ,3*npl,1.0,ctmp2,nmax,
     &           ctmp1,nmax,0.0,covwp1,3   )
c      write(*,*) covwp1(1,1),covwp1(1,2),
c     &           covwp1(3,1),covwp1(2,2),
c     &           covwp1(2,3),covwp1(3,3)
      eulon=atan2(wp1(2),wp1(1))*rad2dg
      aa=sqrt(wp1(1)**2+wp1(2)**2)
      eulat=atan2(wp1(3),aa)*rad2dg
      rot=sqrt(wp1(1)**2+wp1(2)**2+wp1(3)**2)
c
c --- converte rot da rad/yr a deg/my
c
      rot=rot*rad2dg*1.E6
      if(eulon.gt.360.0) eulon=eulon-360.0
c
c ----- ruota le componenti del vettore di rotazione wp..... in riferimento ENU  ----
c
      rotmat(1,1)=-sin(dg2rad*eulon)
      rotmat(1,2)= cos(dg2rad*eulon)
      rotmat(1,3)=0.0
c
      rotmat(2,1)=-sin(dg2rad*eulat)*cos(dg2rad*eulon)
      rotmat(2,2)=-sin(dg2rad*eulat)*sin(dg2rad*eulon)
      rotmat(2,3)= cos(dg2rad*eulat)
c
      rotmat(3,1)=-cos(dg2rad*eulat)*cos(dg2rad*eulon)
      rotmat(3,2)=-cos(dg2rad*eulat)*sin(dg2rad*eulon)
      rotmat(3,3)=-sin(dg2rad*eulat)
      call sgemm('N','N',3,1,3,1.0,rotmat,3,wp1,3,0.0,wprot,3)
c
c ---- .... e la sua matrice di covarianza in coordinate ENU -----
c
      call sgemm ('N','N',3,3,3,1.0,rotmat,3,covwp1,3,0.0,rotmp,3)
      call sgemm ('N','T',3,3,3,1.0,rotmp,3,rotmat,3,0.0,covwprot,3)
c
      do ii=1,2
       do jj=1,2
         cov2d(ii,jj)=covwprot(ii,jj)
       enddo
      enddo
c
c ---- Calcolo eigenvalues/vectors with Numerical Recipes --
c      (deprecated: not compiling well with gfortran 21 apr 2011)
c
c      call jacobi(cov2d,2,2,eigval,eigvect,nrot)
c      call eigsrt(eigval,eigvect,2,2)
c      az=datan2(eigvect(1,1),eigvect(2,1))*rad2dg
c      aa=dabs(wprot(3))
c      smaj=(sqrt(eigval(1))/aa)*rad2dg*sqrt(2.D0)
c      smin=(sqrt(eigval(2))/aa)*rad2dg*sqrt(2.D0)
c
c ---- Uses LAPACK subroutines for eigenvalues/vectors
c
      call SSYEV('V','U',2,cov2d,2,eigval,work,nmax,info)
      write(*,*) 'dsyef info =',info
      az=atan2(cov2d(1,2),cov2d(2,2))*rad2dg
      aa=abs(wprot(3))
      smaj=(sqrt(eigval(2))/aa)*rad2dg*sqrt(2.0)
      smin=(sqrt(eigval(1))/aa)*rad2dg*sqrt(2.0)
c --------------------------------------------------
      if(az.gt.90.0) az=az-180.0
      if(az.lt.-90.0) az=az+180.0
c
      rote=wprot(1)*rad2dg*1.E6
      srote=sqrt(covwprot(1,1))*rad2dg*1.E6
      rotn=wprot(2)*rad2dg*1.E6
      srotn=sqrt(covwprot(2,2))*rad2dg*1.E6
      rotu=wprot(3)*rad2dg*1.E6
      srotu=sqrt(covwprot(3,3))*rad2dg*1.E6
c
        if(eulon.ge.180.0) eulon=eulon-360.0
        if(eulon.le.-180.0) eulon=eulon+360.0
        write(60,185) eulat,eulon,rot,smaj,smin,az,srotu,
     &                plid(i),plid(j)
        eulon=eulon-180.0
        if(eulon.ge.180.0) eulon=eulon-360.0
        if(eulon.le.180.0) eulon=eulon+360.0
        eulat=-eulat
        rot=-rot
        az=-az
        write(60,185) eulat,eulon,rot,smaj,smin,az,srotu,
     &                plid(i),plid(j)
  185 format(2F10.3,F8.3,3F8.1,F8.3,2X,2A2,2X,'POLE')
c
      write(70,190) wp1(1)*1.E9,wp1(2)*1.E9,wp1(3)*1.E9,
     &              plid(i),plid(j)
  190 format(3F14.6,2X,2A2,2X,'ROT')
      write(70,200) covwp1(1,1)*1.E18,covwp1(1,2)*1.E18,
     &              covwp1(3,1)*1.E18,covwp1(2,2)*1.E18,
     &              covwp1(2,3)*1.E18,covwp1(3,3)*1.E18,
     &              plid(i),plid(j)
  200 format(6F14.6,2X,2A2,2X,'COVAR')
c
      endif
      enddo
      enddo
      endif
c
      stop
      end
c
c ---------------------------------------------------------------
c
