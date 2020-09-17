           program nfnc
           implicit none
           include 'Glob_decl'
           real*8 ef
           external theta,sem
           include 'Glob_cons'
           print*,'Enter the value of ef'
           read*, ef
           N=50000
           dw=0.0001

           open(unit=20,file='ef.dat',status='unknown')
           write(20,*) ef
           close(20)
!*********Read input***********************************
          open(unit=21,file='par.dat',status='unknown')
          read(21,*)ec,ufc,v2
          read(21,*)nc,nf
          read(21,*)uff,m
          close(21)
          
          call  hartree(ef,1)
!********write output**********************************          
          open(unit=21,file='par.dat',status='unknown')
          write(21,*)ec,ufc,v2
          write(21,*)nc,nf
          write(21,*)uff,m
          write(21,*)'ec,ufc,v2'
          write(21,*)'nc,nf'
          write(21,*)'uff,m'
          close(21)
          open(unit=22,file='outpar.dat',status='unknown')
          write(22,*)nc,nf,m
          close(22)
          open(unit=23,file='Gc.dat',status='unknown')
          do i=-N,N
          w=float(i)*dw 
          write(23,*)w,dc(i)
          enddo
          close(23)
          open(unit=24,file='Gf.dat',status='unknown')
          do i=-N,N
          w=float(i)*dw
          write(24,*)w,df(i)
          enddo
          close(24)
          end program nfnc

!**************************************************
! HILBER TRANSFORMATION
!**************************************************
        complex*16 function sem(z)
        
!       This block calculates Hilbert transform for a semicircular
!       DOS with t^*=t
!       In general sem=8*(z-sqrt(z**2-D**2/4))/D**2
        complex*16 z,z1,z2,sem1,sem2,sz2
        real*8 dband
       
        t=1.0d0

        
        dband=t
        z1=dcmplx(dband,0.d0)

        z2=z**2-z1**2

        sz2=zsqrt(z2)

        sem1=2.d0/(z+sz2)

        sem2=2.d0/(z-sz2)


        if(dimag(sem1).le.0.d0) then
             sem=sem1
        else if(dimag(sem2).le.0.d0) then
             sem=sem2
        else
             write(6,*) 'no causal root found'
        end if

        return
        end

!**************************************************
! STEP FUNCTION
!**************************************************
        real*8 function theta(x)
        real*8 x

        if(x.lt.0.d0) then
             theta=0.d0
        else if(x.eq.0.d0) then
             theta=0.5d0
        else
             theta=1.d0
        end if
        return
        end function theta

!**************************************************
!SUBROUTINE FOR HARTREE_FOCK
!**************************************************
        subroutine hartree(ef,flag)
        include 'Glob_decl'
        integer miter,flag
        real*8 ef,dop
        real*8 conv,mconv,x
        real*8 pnf,pnc,pmag,dctot,dftot
        real*8 r,r1,r2,r3,r4,r5
        real*8 f,f1,f2
        real*8 a1,a2,alpha,beta
        real*8 c1,c2,c3
        real*8 theta
        real*8 dcup(-Nm:Nm),dcdn(-Nm:Nm),dfup(-Nm:Nm),dfdn(-Nm:Nm)
        complex*16 sem
        complex*16 z,phyb,nhyb
        complex*16 gamma,gup,gdn,g0,hyb(-Nm:Nm) 
        complex*16 sbar(-Nm:Nm)
        complex*16 gc(-Nm:Nm),gcup(-Nm:Nm),gcdn(-Nm:Nm)
        complex*16 gf(-Nm:Nm),gfup(-Nm:Nm),gfdn(-Nm:Nm)
        external sem,theta
        include 'Glob_cons'
        dop=0.99d0
        N=50000
        dw=0.0001
        v=dsqrt(v2) 
        a1=t**2/4.d0
        a2=t**4/8.d0
        alpha=a1
        beta=a2-a1**2
        
        mconv=1.d0
        miter=1
        pmag=m
        pnf=nf
        pnc=nc
        
       do while (mconv.gt.1.D-5.and.miter.le.50)
            r=0.d0
           r1=0.d0
           r2=0.d0
           r3=0.0d0
           r4=0.0d0
           r5=0.0d0
        do i=-N,N
          w=dfloat(i)*dw
          z=w+ii*1.D-7 
          hyb(i)=(t**2/4.d0)*sem(z)
          phyb=hyb(i)
         gup=z-(ec+ufc*nf)-v**2/(z-(ef+uff*nf/2.d0  &
                 +ufc*nc-uff*m/2.d0))
         gdn=z-(ec+ufc*nf)-v**2/(z-(ef+uff*nf/2.d0 & 
                   +ufc*nc+uff*m/2.d0))
         g0=z-(ec+ufc*nf)
         conv=1.d0
        iter=1
        do while (conv.ge.1.D-10.and.iter.le.100000)
        phyb=hyb(i)
!        gamma=(2.d0*gup*gdn-hyb(i)*(gup+gdn))/(gup   &
!            +gdn-2.d0*hyb(i))
        
           gamma=hyb(i)+(2.d0*(g0-hyb(i))*(gup*gdn  & 
                   -hyb(i)*(gup+gdn)+hyb(i)**2))/  &
                 ((1.d0-dop)*(g0-hyb(i))* &
                  (gup+gdn-2.d0*hyb(i))+2.d0*dop* &
                (gup*gdn-hyb(i)*(gup+gdn)+hyb(i)**2))
!         gamma=(-hyb(i)*gup*g0-hyb(i)*gdn*g0+(i)**2*gup &
!	     +hyb(i)**2*gdn-2.0d0*gup*gdn*hyb(i)+2.0d0*g0*gup*gdn &
!             +dop*(-hyb(i)*gup*g0 &
!	     -hyb(i)*gdn*g0-hyb(i)**2*gup-hyb(i)**2*gdn &
!	     +2.0d0*hyb(i)**2*g0+2.0d0*hyb(i)*gup*gdn))/(gup*g0 &
!	     +gdn*g0-2.0d0*hyb(i)*g0-hyb(i)*gup-hyb(i)*gdn &
!             +2.0d0*hyb(i)**2+dop*(-gup*g0 &
!             -gdn*g0+2.0d0*hyb(i)*g0+2*gup*gdn-hyb(i)*gup-hyb(i)*gdn))

        if(zabs(gamma).le.1.D3)then
        nhyb=gamma-1.d0/sem(gamma)
        else
        nhyb=alpha/gamma+abeta/gamma**3
        endif
        hyb(i)=0.5d0*(phyb+nhyb)
        conv=zabs(hyb(i)-phyb)
        iter=iter+1
        enddo
        sbar(i)=gamma
        gc(i)=1.d0/(sbar(i)-hyb(i))
        gfup(i)=1.d0/(z-ef-uff*nf/2.d0-ufc*nc+uff*m/2.d0  &
                -v**2/(z-(ec+ufc*nf)-hyb(i)))

        gfdn(i)=1.d0/(z-ef-uff*nf/2.d0-ufc*nc-uff*m/2.d0  &
                   -v**2/(z-(ec+ufc*nf)-hyb(i)))
        gf(i)=0.5d0*(gfup(i)+gfdn(i))

!***************Spectral Functions****************************
       dc(i)=-dimag(gc(i))/pi
       dfup(i)=-dimag(gfup(i))/pi 
       dfdn(i)=-dimag(gfdn(i))/pi
       df(i)=-dimag(gf(i))/pi
       r=r+2.d0*dc(i)*dw*theta(-w)
       r1=r1+dfup(i)*dw*theta(-w)
       r2=r2+dfdn(i)*dw*theta(-w)
       r3=r3+dc(i)*dw
       r4=r4+df(i)*dw
       enddo
       f=r
       f1=r1
       f2=r2
       nc=f
       nf=f1+f2
       m=f1-f2
       dctot=r3
       dftot=r4
!       c1=dabs(m-pmag)
       c2=dabs(nf-pnf)
!       c3=dabs(nc-pnc)
       mconv=c2
       miter=miter+1
       pmag=m
       pnf=nf
       pnc=nc
       enddo
       if(flag.eq.1)write(6,*) 'mconv=',mconv
       if(flag.eq.1) write(6,*)'----------------------------'
       if(flag.eq.1) write(6,*)'spectral weight'
       if(flag.eq.1)write(6,*)'nc=',nc,'nf=',nf
       if(flag.eq.1) write(6,*)dctot,dftot

       return
       end 
