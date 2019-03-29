      Program cell2cell
      use fun
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="SnTe"
      integer,parameter::nkpath=5, ndk=100, ncite=4
      real*8,parameter::ef=0.4117865695*13.6 !7.434771770
!------------------------------------------------------
      character(len=30)::klabel(nkpath)
      character(len=80) hamil_file,nnkp,win,line,hfile,new_nnkp,w2winup,w2windn
      integer*4,parameter::nk=(nkpath-1)*ndk+1,natom_max=200
      integer*4,parameter::n1=7,n2=7,n3=7,npp=(2*n1+1)*(2*n2+1)*(2*n3+1)
      integer*4 i,j,k,nr,i1,i2,nb,lwork,info,r0,ib,jb,iib,jjb,mb,p,np,n,natom
      integer*4 icite,jcite,rr(3,1),rt0(3),rexist(-n1:n1,-n2:n2,-n3:n3)
      real*8,parameter::third=1d0/3d0, pi2=4.0d0*atan(1.0d0)*2.0d0
      real*8 phase,jk,a,b,rf(3,natom_max),acell2bcell(3,3),bcell2acell(3,3)
      real*8 klist(3,1:nk),xk(nk),kpath(3,ndk),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath)
      real*8 avec(3,3),atranvec(3,3),btranvec(3,3),rrr(3,1),rcite(3,ncite),rrt(3),r1(3,1),r2(3,1)
      real*8,allocatable::ene(:,:),rwork(:)
      integer*4,allocatable:: rvec(:,:),rrvec(:,:),rv(:,:),ndeg(:),atype(:)
      integer*4,allocatable:: ibind(:,:),jbind(:,:),rt(:,:,:)
      complex*16,allocatable:: Hk(:,:),Hamr(:,:,:),work(:),Hr(:,:,:),Hrt(:,:,:)
      integer*4 lmn(-30:30,-30:30,-30:30)
      character(len=4)label(natom_max)
!----------------taarget lattice--------------------------------------
      data atranvec(:,1)   /   3.1550000d0, 0.0000000d0, 0.0000000d0 /
      data atranvec(:,2)   /   0.0000000d0, 3.1550000d0, 0.0000000d0 /
      data atranvec(:,3)   /   0.0000000d0, 0.0000000d0, 3.1550000d0 /
           atranvec=atranvec*2.d0

      data rcite(:,1)      /   0d0        , 0d0        , 0d0         /
      data rcite(:,2)      /   1d0        , 0d0        , 0d0         /
      data rcite(:,3)      /   0d0        , 1d0        , 0d0         /
      data rcite(:,4)      /   0d0        , 0d0        , 1d0         /

!----------------------kpath------------------------------------------
!      data kpath(:,1) /     0.0d0,    0d0,    0.0d0/  !G
!      data kpath(:,2) /     0.5d0,  0.0d0,    0.0d0/  !X

      data kpath(:,1) /     0.5d0,  0.5d0,    0.5d0/  !R
      data kpath(:,2) /     0.0d0,    0d0,    0.0d0/  !G
      data kpath(:,3) /     0.5d0,  0.0d0,    0.0d0/  !X
      data kpath(:,4) /     0.5d0,  0.5d0,    0.0d0/  !M
      data kpath(:,5) /     0.0d0,    0d0,    0.0d0/  !G

!      data klabel     /'{/Symbol \107}','X'/
      data klabel     /'R','{/Symbol \107}','X','M','{/Symbol \107}'/

!---------------------------------------------------------------------
      write(hamil_file,'(a,a)')trim(adjustl(prefix)),"_hr.dat"
      write(hfile,'(a,a)')     trim(adjustl(prefix)),"_new_hr.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"
      write(win,'(a,a)')       trim(adjustl(prefix)),".win"
      write(w2winup,'(a,a)')       trim(adjustl(prefix)),".w2winup"
      write(w2windn,'(a,a)')       trim(adjustl(prefix)),".w2windn"
      write(new_nnkp,'(a,a)')  trim(adjustl(prefix)),"_new.nnkp"

!--------------- real and  reciprocal vectors-------------------------
      open(98,file=trim(adjustl(nnkp)),err=333)

109   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin real_lattice") goto 109
      read(98,*)avec
      close(98)

      bvec = reciprocal(atranvec)
      acell2bcell=matmul(matinv3(atranvec),avec)
      bcell2acell=matmul(matinv3(avec),atranvec)
      
      open(198,file=trim(adjustl(new_nnkp)))

      write(198,'(a)') 'Created by M. S. Bahramy'
      write(198,*)

      write(198,'(a)') 'calc_only_A  :  F'
      write(198,*)

      write(198,'(a)')'begin real_lattice'
      write(198,'(3f12.7)')atranvec
      write(198,'(a)')'end real_lattice'
      write(198,*)

      write(198,'(a)')'begin recip_lattice'
      write(198,'(3f12.7)')bvec
      write(198,'(a)')'end recip_lattice'
      write(198,*)
      close(198)
!------------------ atomic positions in origianl cell ----------------
      open(89,file=trim(adjustl(win)),err=334)
      natom=0
110   read(89,'(a)')line
      if(trim(adjustl(line)).ne."begin atoms_frac") goto 110
      do while(natom.le.natom_max)
      read(89,'(a)')line
      if(trim(adjustl(line)).eq."end atoms_frac") exit
      natom=natom+1
      read(line,*) label(natom),rf(1:3,natom)
      enddo
      close(89)
      if(trim(adjustl(line)).ne."end atoms_frac") goto 335
!-------------------kmesh---------------------------------------------
      ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+&
                (kpath(2,1)-kpath(2,2))*bvec(:,2)+&
                (kpath(3,1)-kpath(3,2))*bvec(:,3)

      xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
      xkl(1)=xk(1)
      
      k=0
      ktemp1=0d0
      do i=1,nkpath-1
       do j=1,ndk
        k=k+1
        jk=dfloat(j-1)/dfloat(ndk)
        klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
        ktemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
        if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
        if(j.eq.1) xkl(i)=xk(k)
        ktemp1=ktemp2
       enddo
      enddo
      klist(:,nk)=kpath(:,nkpath)
      ktemp2=klist(1,nk)*bvec(:,1)+klist(2,nk)*bvec(:,2)+klist(3,nk)*bvec(:,3)
      xk(nk)=xk(nk-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
      xkl(nkpath)=xk(nk)
      klist=klist*pi2

!------read H(R)
      open(99,file=trim(adjustl(hamil_file)),err=444)

      lmn=0
      read(99,*)
      read(99,*)nb,nr
      mb=nb*ncite
      allocate(rvec(3,nr) ,Hamr(nb,nb,nr),ndeg(nr))
      allocate(rrvec(3,npp),  Hrt(mb,mb,npp),Hk(mb,mb),ene(mb,nk))
      read(99,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)rvec(1,k),rvec(2,k),rvec(3,k),i1,i2,a,b
               hamr(i1,i2,k)=dcmplx(a,b)
            enddo
         enddo
           lmn(rvec(1,k),rvec(2,k),rvec(3,k))=k
      enddo
      close(99)
!------------------ find wannier centers in the original cell --------
      open(90,file=trim(adjustl(w2winup)),err=336)
      read(90,'(a)')line
      read(90,'(a)')line
      read(90,*) i, i
      if(i.ne.nb) goto  346
      allocate(atype(nb))
      do i=1,nb
         read(90,*)i1
         do j=1,i1
            read(90,*)atype(i)
         enddo
      enddo
      close(90)
      open(90,file=trim(adjustl(w2windn)),err=337)
      read(90,'(a)')line
      read(90,'(a)')line
      read(90,*) i, i
      if(i.ne.nb) goto  347
      do i=1,nb
         read(90,*)i1
         do j=1,i1
            read(90,*)atype(i)
         enddo
      enddo
      close(90)
      r1(1,1)=0.5d0
      r1(2,1)=-0.5d0
      r1(3,1)=-0.5d0
      write(*,*)matmul(bcell2acell,r1)
!-------------- find the appropriate atoms residing inside the new unitcell and index them
      allocate(ibind(mb,mb),jbind(mb,mb),rt(3,mb,mb))
      rt=0
      ibind=0
      jbind=0
      do icite =1,ncite
         do jcite=1,ncite
            do i1=1,nb
               do i2=1,nb
!--------------------------transform
                  r1(:,1)=rcite(:,icite)+rf(:,atype(i1))
                  r2(:,1)=rcite(:,jcite)+rf(:,atype(i2))
                  r1     =matmul(acell2bcell,r1)
                  r2     =matmul(acell2bcell,r2)
                  r1     =mod(r1+1000d0,1d0)        ! the offset 1000 added to ensure all r vectors are positive
                  r2     =mod(r2+1000d0,1d0)       
                  r1     =matmul(bcell2acell,r1)
                  r2     =matmul(bcell2acell,r2)
                  r1(:,1)=r1(:,1)-rf(:,atype(i1))
                  r2(:,1)=r2(:,1)-rf(:,atype(i2))
                  rrr    =r2-r1
!--------------------------upup
                  if    (i1.le.nb/2.and.i2.le.nb/2) then
                           iib=i1+(      icite-1)*nb/2
                           jjb=i2+(      jcite-1)*nb/2
!--------------------------updn
                  elseif(i1.le.nb/2.and.i2.gt.nb/2) then
                           iib=i1+(      icite-1)*nb/2
                           jjb=i2+(ncite+jcite-2)*nb/2
!--------------------------dnup
                  elseif(i1.gt.nb/2.and.i2.le.nb/2) then
                           iib=i1+(ncite+icite-2)*nb/2
                           jjb=i2+(      jcite-1)*nb/2
!--------------------------dndn
                  elseif(i1.gt.nb/2.and.i2.gt.nb/2) then
                           iib=i1+(ncite+icite-2)*nb/2
                           jjb=i2+(ncite+jcite-2)*nb/2
                  endif
!--------------------------finally index
                  ibind(iib,jjb)=i1
                  jbind(iib,jjb)=i2
                  rt(:,iib,jjb) =int(rrr(:,1))
               enddo
            enddo
         enddo
      enddo
!-------construct new H(R)
      p=0
      Hrt=dcmplx(0d0,0d0)
      rexist=0
      do i=-n1,n1
         do j=-n2,n2
            do k=-n3,n3
               p=p+1
               rrvec(1,p)=i
               rrvec(2,p)=j
               rrvec(3,p)=k
               rr(1,1)   =i
               rr(2,1)   =j
               rr(3,1)   =k
               rrr=matmul(bcell2acell,rr)
               do i1=1,mb
                  do i2=1,mb
                     rt0=rrr(:,1)+rt(:,i1,i2)
                     r0=lmn(rt0(1),rt0(2),rt0(3))
                     if(r0.eq.0)cycle
                     rexist(i,j,k)=1
                     ib=ibind(i1,i2)
                     jb=jbind(i1,i2)
                     Hrt(i1,i2,p)=Hamr(ib,jb,r0)
                  enddo
               enddo
            enddo
         enddo
      enddo 
!-------------------------clean up Hr
      np=sum(rexist)
      open(199,file=trim(adjustl(hfile)))
      write(199,'(a)') 'Created by M. S. Bahramy'
      write(199,'(i12)')mb,np
      write(199,'(15i5)') (i/i,i=1,np)
      allocate(Hr(mb,mb,np),rv(3,np))
      p=0
      n=0
      Hr=dcmplx(0d0,0d0)
      do i=-n1,n1
         do j=-n2,n2
            do k=-n3,n3
               p=p+1
               if(rexist(i,j,k).eq.0)cycle
               n=n+1
               Hr(:,:,n)=Hrt(:,:,p)
               rv(:,n)  =rrvec(:,p)
               do jb=1,mb
                  do ib=1,mb
                      write(199,"(5i5,2f12.8)")rv(1:3,n),ib,jb,Hr(ib,jb,n)
                  enddo
               enddo
            enddo
         enddo
      enddo


!---- Fourrier transform H(R) to H(k)
     lwork=max(1,2*mb-1)
     allocate(work(max(1,lwork)),rwork(max(1,3*mb-2)))

      ene=0d0
      do k=1,nk
         HK=(0d0,0d0)
         do j=1,np

            phase=0.0d0
            do i=1,3
               phase=phase+klist(i,k)*rv(i,j)
            enddo

            do i1=1,mb
               do i2=1,mb
                  Hk(i1,i2)=Hk(i1,i2)+Hr(i1,i2,j)* &
                  dcmplx(cos(phase),sin(phase))!/float(ndeg(j))
               enddo
            enddo

         enddo
         call zheev('V','U',mb,Hk,mb,ene(:,k),work,lwork,rwork,info)
      enddo

!-------export E(k)     
      open(100,file='band.dat')
      do i=1,mb
         do k=1,nk
           write(100,'(2(x,f12.6))') xk(k),ene(i,k)
         enddo
           write(100,*)
      enddo
      close(100)
      call write_plt_rgb(nkpath,xkl,klabel,ef)
      stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
334   write(*,'(3a)')'ERROR: input file "',trim(adjustl(win)),' not found'
      stop
335   write(*,'(2a)')'ERROR: reached the end of "',trim(adjustl(win))
      stop
336   write(*,'(3a)')'ERROR: input file "',trim(adjustl(w2winup)),' not found'
      stop
337   write(*,'(3a)')'ERROR: input file "',trim(adjustl(w2windn)),' not found'
      stop
346   write(*,'(2a)')'ERROR: incompatible "',trim(adjustl(w2winup))
      stop
347   write(*,'(2a)')'ERROR: incompatible "',trim(adjustl(w2windn))
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file)),' not found'
      stop
      end
!--------------------------------------------------------------------------------

     subroutine write_plt_rgb(nkp,xkl,kl,ef)
     implicit none
     integer nkp,i
     real*8 xkl(nkp),ef,kz
     character(len=30)kl(nkp)
     
     open(99,file='band.plt')
     write(99,'(a,f12.8)')'ef=',ef
     write(99,'(a)') 'set xtics ( \'
     do i=1,nkp
        if(trim(adjustl(kl(i))).eq.'g'.or.trim(adjustl(kl(i))).eq.'G')kl(i)="{/Symbol \107}"
        if(i.ne.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i),", \"
        if(i.eq.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i)," )"
     enddo
     write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [',xkl(1),':',xkl(nkp),']'
!     write(99,'(a,f4.2,a)')'set title "k_z=',kz,'{/Symbol \160}/c"'
     write(99,'(a)') &
          'set terminal pdfcairo enhanced font "Helvetica"  transparent fontscale 1 size 5.00in, 5.00in'
     write(99,'(a)')&
          'set output "band.pdf"'
     write(99,'(14(a,/),a)') &
          'set style line 10 lt 1 lc rgb "black" lw 2',&
          'set border ls 10',&
          'set tics textcolor rgb "black"',&
          'set encoding iso_8859_1',&
          'set size ratio 0 1.0,1.0',&
          'set ylabel "E-E_{F} (eV)"',&
          'set yrange [ -6.0 : 4 ]',&
          'unset key',&
          'set ytics 1.0 scale 1',&
          'set mytics 5',&
          'set parametric',&
          'set trange [-10:10]',&
          'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)',&
          'plot \',&
          '"band.dat" using 1:($2-ef) with l lt 1 lw 4 lc 3, \'

    do i=2,nkp-1
      write(99,'(f12.6,a)') xkl(i),',t with l ls 10,\'
    enddo
    write(99,'(a)') 't,0 with l ls 10'
    end subroutine write_plt_rgb
