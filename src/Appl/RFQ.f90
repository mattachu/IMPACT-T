!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! RFQclass: RFQ cavity class in Lattice module of APPLICATION layer.
! Version:  2.0
! Author:   Zhi Wang, Peking University
! Description: This class defines the linear transfer map and field
!              for the RFQ cavity.
! Comments:
!----------------------------------------------------------------
      module RFQclass
        use PhysConstclass
        use Dataclass
        integer, private, parameter :: Nparam = 26

        type RFQ
          !Itype = 21
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
            end type RFQ
          ! Param(1) : zedge
          !      (2) : scale
          !      (3) : frequency
          !      (4) : voltage
          !      (5) : min-a
          !      (6) : modulation
          !      (7) : synchronous-fais

          !      (8) : multipole field1 -A10
          !      (9) : multipole field2 -A0(A01)
          !      (10) : multipole field1 -A12
          !      (11) : multipole field2 -A1(A03)
          !      (12) : multipole field1 -A30
          !      (13) : multipole field2 -A21
          !      (14) : multipole field2 -A32
          !      (15) : multipole field1 -A23

          !   below are former cell coefficients
          !      (16) : multipole field1 -A10
          !      (17) : multipole field2 -A0(A01)
          !      (18) : multipole field1 -A12
          !      (19) : multipole field2 -A1(A03)
          !      (20) : multipole field1 -A30
          !      (21) : multipole field2 -A21
          !      (22) : multipole field2 -A32
          !      (23) : multipole field1 -A23


        interface getparam_RFQ
        module procedure getparam1_RFQ,  &
                          getparam2_RFQ,   &
                          getparam3_RFQ
         end interface
        interface setparam_RFQ
        module procedure setparam1_RFQ,  &
                          setparam2_RFQ, setparam3_RFQ
        end interface
        contains
        subroutine construct_RFQ(this,numseg,nmpstp,type,blength)
        implicit none
        type (RFQ), intent(out) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength

        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0

        end subroutine construct_RFQ

        subroutine setparam1_RFQ(this,i,value)
        implicit none
        type (RFQ), intent(inout) :: this
        integer, intent(in) :: i
        double precision, intent(in) :: value

        this%Param(i) = value

        end subroutine setparam1_RFQ

        subroutine setparam2_RFQ(this,values)
        implicit none
        type (RFQ), intent(inout) :: this
        double precision, dimension(:), intent(in) :: values

        this%Param = values

        end subroutine setparam2_RFQ

        subroutine setparam3_RFQ(this,numseg,nmpstp,type,blength)
        implicit none
        type (RFQ), intent(inout) :: this
        integer, intent(in) :: numseg,nmpstp,type
        double precision, intent(in) :: blength

        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        end subroutine setparam3_RFQ


        subroutine getparam1_RFQ(this,i,blparam)
        implicit none
        type (RFQ), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_RFQ

        subroutine getparam2_RFQ(this,blparams)
        implicit none
        type (RFQ), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_RFQ

        subroutine getparam3_RFQ(this,blength,bnseg,bmapstp,&
                                       btype)
        implicit none
        type (RFQ), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype
        end subroutine getparam3_RFQ

        subroutine  getfld_RFQ(pos,momentum,extfld,this)
        implicit none
        include 'mpif.h'
        integer:: bnseg,N,i,j,ii,jj,nsg,cellNum,numRMS
        type (RFQ), intent(in) :: this
        double precision, dimension(4), intent(in) :: pos
        double precision, dimension(3), intent(in) :: momentum
        double precision, dimension(6), intent(out) :: extfld
        double precision:: zz,zzm,bgrad,zedge,interV,modM,synP,escale,freQ,betax,betay,betaz
        double precision:: exj,eyj,ezj,aj,bj,cj,dj,aaj,bbj,ccj,ddj,rsqq,kb,kbeta
        double precision:: ww,tt,apx,apy,apz,apx1,apx2,apx3,apx4,apy1,apy2,apy3,apy4,apz1,apz2,apz3,apz4,apz5
        double precision:: tmpsin,tmpcos,rsq,theta,tmpex,tmpey,tmpez,tmpbx,tmpby,tmpbz
        double precision:: am,bm,blength,blength2,fai0,fais,fai,arg,rbsq,xI0
      double precision:: k,kk,X,Ex,Ey,Ez,Ipd0,Ipd2,Ipd4,Ipd6,Ipd0p,Ipd2p,Ipd4p,Ipd6p
      double precision:: a01,a03,a10,a12,a21,a23,a30,a32,minA,midA,betaS
        double precision:: a01f,a03f,a10f,a12f,a21f,a23f,a30f,a32f
        double precision:: a01k,a03k,a10k,a12k,a21k,a23k,a30k,a32k
        double precision:: Focu,Acc,Kfocu,Kacc,Ka,Kf,km,ki,apxte,apyte,exte,eyte
        double precision:: a10kk,a01kk,a03kk,a30kk,a21kk,a12kk,a23kk,a32kk
        double precision:: sin1t,cos1t,cos2t,sin2t,sin4t,cos4t,cos6t,sin6t
        double precision::apx1t,apx2t,apx3t,apx4t,apy1t,apy2t,apy3t,apz1t,apz2t,apz3t,xa,apxr,apyr,apzr
        double precision::Lr,Lf,EXT,EYT,EZT,apxr1,apxr2,apxr3,apxr4,apyr1,apyr2,apyr3,apyr4,apzr1,apzr2,apzr3,apzr4
         double precision::kt,lt
         double precision ::zt
         double precision :: ml,al
         double precision :: ddk
          integer :: ierr, myid
         double precision :: T10k,T30k,m,A10T,A30T,betaS2,betaSi
         double precision :: kka
        ! length unit is meter
        ! 2           3             4        6       7          8        9         10      11        12       13
        ! Z           a             m       thi     fir       four      sec       sev     five       eig     six
        ! Length(cm)  Aper(cm)    modu      A10   A0(A01)    A12*I4   A1(A03)    A30*I0   A21*I2   A32*I4   A23*I6

        bnseg = this%Nseg
        blength = this%Length
        zedge = this%Param(1)
        blength2 =this%Param(2)
       ! betaS2=this%Param(2)
        cellNum = this%Param(3)
        Lr=this%Param(4)
        minA= 0.01*this%Param(5)
        midA= 0.01*this%Param(6)
        interV=this%Param(7)
        betaS=this%Param(8)
        a10 = this%Param(9)   !acc term
        a10f  = this%Param(10)
        a01 = this%Param(11)  !focus term
        a01f  = this%Param(12)
        !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 2-term potentinal expression
        a12  = this%Param(13)
        a12f  = this%Param(14)
        a03  = this%Param(15)
        a03f  = this%Param(16)
        a30  = this%Param(17)
        a30f  = this%Param(18)
        a21  = this%Param(19)
        a21f  = this%Param(20)
        a32  = this%Param(21)
        a32f  = this%Param(22)
        a23  = this%Param(23)
        a23f  = this%Param(24)

        ml=this%Param(25)
        al=this%Param(26)

        extfld(1) = 0.0
        extfld(2) = 0.0
        extfld(3) = 0.0
        extfld(4) = 0.0
        extfld(5) = 0.0
        extfld(6) = 0.0

        xa=0.002
        ww = 2*Pi*scfreq
        zz=pos(3)
        zzm=pos(3)-zedge
        betax=momentum(1)
        betay=momentum(2)
        betaz=momentum(3)

          !LHP
        ki=pi/blength
        kb=pi/blength2
        kk=ki+zzm*(kb-ki)/blength
        ddk=(kb-ki)/blength

        if (this%MAPSTP.le.0) then
          kk=ki
       !   kb=pi/2*blength2
       !   kk=ki+zzm*(kb-ki)/blength
        endif

       !! kbeta=2*pi/(betaz*clight/scfreq)

       ! Lt=zedge
     !  betaSi=zzm*(betaS-betaS2)/blength+betaS
    !  kk=2*pi/(betaSi*clight/scfreq)
      ! kk=ki

        rsq=sqrt(pos(1)*pos(1)+pos(2)*pos(2))
        rbsq=pos(1)*pos(1)+pos(2)*pos(2)
        tt = pos(4)
        fai0=0 !for beam center at -Lc/2
        fais=(ww*tt+fai0)

        tmpcos = cos(fais)
        tmpsin = sin(fais)

        a01k=(a01-a01f)/blength
        a03k=(a03-a03f)/blength
        a10k=(a10-a10f)/blength
        a12k=(a12-a12f)/blength
        a21k=(a21-a21f)/blength
        a23k=(a23-a23f)/blength
        a30k=(a30-a30f)/blength
        a32k=(a32-a32f)/blength
             !2015-12-4
        A10KK=A10f+A10K*ZZm
        a01kk=a01f+a01k*zzm
        A12KK=A12f+A12K*ZZm
        a03kk=a03f+a03k*zzm
        A21KK=A21f+A21K*ZZm
        a23kk=a23f+a23k*zzm
        A30KK=A30f+A30K*ZZm
        a32kk=a32f+a32k*zzm

        if (rsq .ne.0) then
        cos1t=pos(1)/rsq
        sin1t=pos(2)/rsq
        else
            cos1t=atan(1.0)
            sin1t=atan(1.0)
        endif

        sin2t=2*sin1t*cos1t
        cos2t=2*cos1t*cos1t-1
        sin4t=-4*sin1t*cos1t*(2*sin1t*sin1t-1)
        cos4t=1+8*(cos1t**4-cos1t*cos1t)
        sin6t=2*sin1t*cos1t*(2*sin1t+1)*(2*sin1t-1)*(-3+4*sin1t*sin1t)
        cos6t=(2*cos1t*cos1t-1)*(16*cos1t**4-16*cos1t*cos1t+1)

        apx1=0
        apx2=0
        apx3=0
        apx4=0
        apz1=0
        apz2=0
        apz3=0
        apz4=0
        apz5=0

        if (xa.eq.0)then
        apxr=0
        apyr=0
         endif

        if(bnseg.eq.-1)then
       !Lr=0
    !   midA=0.4421*0.01
     !  midA=0.557600*0.01
      ! a01=1.0
       km=pi/(2*Lr)
       kk=km
        apx=   6*pos(1)*a01*(cos(km*(zz-Lr))/4+cos(3*km*(zz-Lr))/12)/(midA**2)+&
        34560*A03*rsq**5*(cos(km*(zz-Lr))+cos(3*km*(zz-Lr))/3)*(cos6t*cos1t+sin6t*sin1t)/(7680*midA**6)
        apy=   6*pos(2)*a01*(cos(km*(zz-Lr))/4+cos(3*km*(zz-Lr))/12)/(midA**2)+&
        34560*A03*rsq**5*(cos(km*(zz-Lr))+cos(3*km*(zz-Lr))/3)*(-cos6t*sin1t+sin6t*cos1t)/(7680*midA**6)



       ! apz=6*a01*km*cos2t*rbsq/8.0*(sin(km*(zz-Lr))+sin(3*km*(zz-Lr)))/(midA**2)
        apz=0



         elseif(bnseg.eq.1) then
            !for test
        Lr=0


        !r direction field
         if (rsq.eq.0) then
             apx=0.0
             apy=0.0
         else


        apx1= 2*pos(1)*a01kk/(midA**2)+a10kk*kk*I0p(kk*rsq)*cos1t*cos(ki*(zzm-Lr))*(-1)**(cellnum+1)
       ! apx1= 2*pos(1)*a01kk/(midA**2)+a10kk*ki*I0p(kk*rsq)*cos1t*cos(ki*(zzm-Lr))*(-1)**(cellnum)
        apy1= 2*pos(2)*a01kk/(midA**2)-a10kk*kk*I0p(kk*rsq)*sin1t*cos(ki*(zzm-Lr))*(-1)**(cellnum+1)

       ! apx1= 2*pos(1)/(midA**2)+a10kk*kk*I0p(kk*rsq)*cos1t*cos(ki*(zzm-Lr))*(-1)**(cellnum)
       ! apx1= 2*pos(1)*a01kk/(midA**2)+a10kk*ki*I0p(kk*rsq)*cos1t*cos(ki*(zzm-Lr))*(-1)**(cellnum)
       ! apy1= 2*pos(2)/(midA**2)-a10kk*kk*I0p(kk*rsq)*sin1t*cos(ki*(zzm-Lr))*(-1)**(cellnum)

        apx2= a03kk*6*rsq**5*(cos1t*cos6t+sin1t*sin6t)/(midA**6)+&
        cos(ki*(zzm-Lr))*a12kk*(kk*I4p(kk*rsq)*cos4t*cos1t+4*i4(kk*rsq)*sin4t*sin1t/RSQ)*(-1)**(cellnum+1)
        apy2= a03kk*6*rsq**5*(-sin1t*cos6t+cos1t*sin6t)/(midA**6)+&
        cos(ki*(zzm-Lr))*a12kk*(-kk*I4p(kk*rsq)*cos4t*sin1t+4*i4(kk*rsq)*sin4t*cos1t/RSQ)*(-1)**(cellnum+1)


        apx3=(cos(2*ki*(zzm-Lr))*cos1t*((a21kk*2*kk*I2p(2*kk*rsq)*cos2t+a23kk*2*kk*I6p(2*kk*rsq)*cos6t))+&
        cos(2*ki*(zzm-Lr))*sin1t*(a21kk*2*I2(2*kk*rsq)*sin2t+a23kk*6*I6(2*kk*rsq)*sin6t)/RSQ)
        apy3=(cos(2*ki*(zzm-Lr))*sin1t*((-a21kk*2*kk*I2p(2*kk*rsq)*cos2t-a23kk*2*kk*I6p(2*kk*rsq)*cos6t))+&
        cos(2*ki*(zzm-Lr))*cos1t*(a21kk*2*I2(2*kk*rsq)*sin2t+a23kk*6*I6(2*kk*rsq)*sin6t)/RSQ)

        apx4=(cos(3*ki*(zzm-Lr))*cos1t*((a30kk*3*kk*I0p(3*kk*rsq)+a32kk*3*kk*I4p(3*kk*rsq)*cos4t))+&
        cos(3*ki*(zzm-Lr))*sin1t*a32kk*4*I4(3*kk*rsq)*sin4t/rsq)*(-1)**(cellnum+1)
        apy4=(cos(3*ki*(zzm-Lr))*sin1t*(-a30kk*3*kk*I0p(3*kk*rsq)-a32kk*3*kk*I4p(3*kk*rsq)*cos4t)+&
        cos(3*ki*(zzm-Lr))*cos1t*a32kk*4*I4(3*kk*rsq)*sin4t/rsq)*(-1)**(cellnum+1)

       ! apx2=0
       ! apx3=0
       ! apx4=0
       ! apy2=0
       ! apy3=0
       ! apy4=0

        apx=apx1+apx2+apx3+apx4
        apy=apy1+apy2+apy3+apy4
         endif


        apz1=ki*sin(ki*(zzm-Lr))*(A10KK*I0(kk*rsq)+a12kk*I4(kk*rsq)*cos4t)*(-1)**(cellnum+1)
      !  apz1=ki*sin(ki*(zzm-Lr))*A10KK*I0(kk*rsq)*(-1)**(cellnum+1)

        apz2=-a01k*rsq*rsq*cos2t/midA**2-cos6t*a03k*rsq**6/midA**6  ! ???
        apz3=2*ki*sin(2*ki*(zzm-Lr))*(a21kk*I2(2*kk*rsq)*cos2t+a23kk*I6(2*kk*rsq)*cos6t)
        apz4=3*ki*sin(3*ki*(zzm-Lr))*(a30kk*I0(3*kk*rsq)+a32kk*I4(3*kk*rsq)*cos4t)*(-1)**(cellnum+1)
        apz5=-cos(ki*(zzm-Lr))*(A10k*I0(kk*rsq)+A10kk*I0p(kk*rsq)*ddk*rsq)*(-1)**(cellnum+1)
    !    apz5=100
      !  apz2=0
      !  apz3=0
      !  apz4=0

      !  apz=((-1)**(cellnum+1)*(apz1+apz4)+apz2+apz3)
        apz=apz1+apz4+apz2+apz3+apz5

     !  apz=(ki*sin(ki*(zzm-Lr))*A10KK*I0(kk*rsq))*(-1)**(cellnum+1)
         elseif(bnseg.eq.-2) then
        kt=pi/(blength*2)
        kb=pi/0.053093

        lr=zedge
      !  kk=kb+(kb-kt)/blength*(pos(3)-lr)
      kk=kt
     minA =al/100
     T10k =ml**2*I0(kk*minA)+I0(ml*kk*minA)
     T30k=ml**2*I0(3*kk*minA)+I0(3*ml*kk*minA)
     A10=(ml**2-1)/(T10k+T30k/3*(I0(kk*midA)/I0(3*kk*midA)))
     A30=A10*I0(kk*midA)/(3*I0(3*kk*midA))
     A10=A10f+(A10-A10f)/blength*(pos(3)-lr)
     A30=A30f+(A30-A30f)/blength*(pos(3)-lr)
     !ki=kk-kk/blength*(pos(3)-lr)

        !apx=( cos1t*(2*rsq*cos2t/midA**2-kk*A10*I0p(kk*rsq)*cos(kk*(zz-Lr))-3*kk*A30*I0p(3*kk*rsq)*cos(3*kk*(zz-Lr)))+&
         !     2*sin1t*sin2t*rsq/midA**2)
         !apy= -(sin1t*(2*rsq*cos2t/midA**2-kk*A10*I0p(kk*rsq)*cos(kk*(zz-Lr))-3*kk*A30*I0p(3*kk*rsq)*cos(3*kk*(zz-Lr)))-&
          !    2*cos1t*sin2t*rsq/midA**2)
         !- should be taken
         !apz=-( -kk*A10*I0(kk*rsq)*sin(kk*(zz-Lr))+3*kk*A30*I0(3*kk*rsq)*sin(3*kk*(zz-Lr)))

         apx=2*A01kk*pos(1)/midA**2-cos1t*(kk*A10*I0p(kk*rsq)*cos(kk*(zz-Lr))+3*kk*A30*I0p(3*kk*rsq)&
             *cos(3*kk*(zz-Lr)))*(-1)**(cellnum+1)+A03*6*rsq**5/midA**6*(cos6t*cos1t+sin6t*sin1t)
     !   apx3=(cos(2*kk*(zz-Lr))*cos1t*((a21kk*2*kk*I2p(2*kk*rsq)*cos2t+a23kk*2*kk*I6p(2*kk*rsq)*cos6t))+&
     !   cos(2*kk*(zz-Lr))*sin1t*(a21kk*2*I2(2*kk*rsq)*sin2t+a23kk*6*I6(2*kk*rsq)*sin6t)/RSQ)
     !   apx2=-cos(kk*(zz-Lr))*a12kk*(kk*I4p(kk*rsq)*cos4t*cos1t+4*I4(kk*rsq)*sin4t*sin1t/RSQ)*(-1)**(cellnum+1)
     !    apx4=(cos(3*kk*(zz-Lr))*cos1t*(a32kk*3*kk*I4p(3*kk*rsq)*cos4t)+&
     !   cos(3*kk*(zz-Lr))*sin1t*a32kk*4*I4(3*kk*rsq)*sin4t/rsq)*(-1)*(-1)**(cellnum+1)
     !   apx=apx1+apx2+apx3+apx4


        apy=2*A01*pos(2)/midA**2-sin1t*(kk*A10*I0p(kk*rsq)*cos(kk*(zz-Lr))+3*kk*A30*I0p(3*kk*rsq)&
            *cos(3*kk*(zz-Lr)))*(-1)**(cellnum+1)+A03*6*rsq**5/midA**6*(-cos6t*sin1t+sin6t*cos1t)
         !- should be taken
        apz=-((kk*A10*I0(kk*rsq)*sin(kk*(zz-Lr))+3*kk*A30*I0(3*kk*rsq)*sin(3*kk*(zz-Lr))))!*(-1)**(cellnum+1)!&
          !   +(-12.90009*I0(kk*rsq)*cos(kk*(zz-Lr))+0.13294*I0(3*kk*rsq)*cos(kk*(zz-Lr)))
        !apz=(-1)**(cellnum+1)*(apz1)

         !- should be taken
      !  apzr=(-kk*A10*I0(kk*xa)*sin(kk*(zz-Lr))+3*kk*A30*I0(3*kk*xa)*sin(3*kk*(zz-Lr)))

       ! open(102,file='Lr.txt',status='unknown')
       ! write(102,*) zz,kt,kk,lr,lt,Lr,zz-Lr
        !close(102)

        elseif(bnseg.eq.-3)then
       ! A03=0.033870
        apx= pos(1)/(midA**2)*2*A01+A03*6*rsq**5/midA**6*(cos6t*cos1t+sin6t*sin1t)
        apy= pos(2)/(midA**2)*2*A01+A03*6*rsq**5/midA**6*(cos6t*cos1t+sin6t*sin1t)
        APZ=0

        apxr=xa/(midA**2)
        apyr=xa/(midA**2)
        apzr=0


        elseif(bnseg.eq.-4)then
         !    A03=0.033870
             lf=zedge
            kk=pi/(2*blength)
        apx=6*pos(1)*a01f*(cos(kk*(zz-Lf))/4+cos(3*kk*(zz-Lf))/12)/(midA**2)+&
         34560*A03f*rsq**5*(cos(kk*(zz-Lf))+cos(3*kk*(zz-Lf))/3)*(cos6t*cos1t+sin6t*sin1t)/(7680*midA**6)
        apy=6*pos(2)*a01f*(cos(kk*(zz-Lf))/4+cos(3*kk*(zz-Lf))/12)/(midA**2)+&
        34560*A03f*rsq**5*(cos(kk*(zz-Lf))+cos(3*kk*(zz-Lf))/3)*(-cos6t*sin1t+sin6t*cos1t)/(7680*midA**6)

       ! apx1=6*A01/kk**2/midA**2*(cos2t*cos1t*(kk*I2p(kk*rsq)*cos(kk*(zz-lf))+1/9*kk*I2p(3*kk*rsq)*cos(3*kk*(zz-lf)))+&
       !                           sin2t*sin1t*2/rsq*(I2(kk*rsq)*cos(kk*(zz-lf))+1/27*I2(3*kk*rsq)*cos(3*kk*(zz-lf))))
       ! apx2=34560*A03/kk**6/midA**6*(cos6t*cos1t*(kk*I6p(kk*rsq)*cos(kk*(zz-lf))+1/3**6*kk*I6p(3*kk*rsq)*cos(3*kk*(zz-lf)))+&
       !                           sin6t*sin1t*6/rsq*(I6(kk*rsq)*cos(kk*(zz-lf))+1/3**7*I6(3*kk*rsq)*cos(3*kk*(zz-lf))))
       ! apy1=6*A01/kk**2/midA**2*(-cos2t*cos1t*(kk*I2p(kk*rsq)*cos(kk*(zz-lf))-1/9*kk*I2p(3*kk*rsq)*cos(3*kk*(zz-lf)))-&
       !                           sin2t*cos1t*2/rsq*(I2(kk*rsq)*cos(kk*(zz-lf))-1/27*I2(3*kk*rsq)*cos(3*kk*(zz-lf))))
       ! apy2=34560*A03/kk**6/midA**6*(-cos6t*sin1t*(kk*I6p(kk*rsq)*cos(kk*(zz-lf))-1/3**6*kk*I6p(3*kk*rsq)*cos(3*kk*(zz-lf)))-&
       !                           sin6t*cos1t*6/rsq*(I6(kk*rsq)*cos(kk*(zz-lf))-1/3**7*I6(3*kk*rsq)*cos(3*kk*(zz-lf))))
       ! apx2=0
       ! apy2=0
       ! apx=(apx1+apx2)*(-1)**(cellnum+1)
       ! apy=(apy1+apy2)*(-1)**(cellnum+1)

        apz=0
        apzr=0
        endif

        if (zz .le. 0) then
            apx=0.0
            apy=0.0
          !  apz=0.0
        endif
        Ex= -interV*1000*apx/2
        Ey= interV*1000*apy/2
        Ez=  interV*1000*apz/2

      !  Ext= -interV*1000*apxR/2
      !  Eyt=  interV*1000*apyR/2
      !  Ezt=  interV*1000*apzR/2

        tmpex=ex*tmpcos
        tmpey=ey*tmpcos
        tmpez=ez*tmpcos
        tmpbx = 0.0
        tmpby = 0.0
        tmpbz = 0.0

        extfld(1) = extfld(1) + tmpex
        extfld(2) = extfld(2) + tmpey
        extfld(3) = extfld(3) + tmpez
        extfld(4) = 0
        extfld(5) = 0
        extfld(6) = 0

        call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
        if (myid.eq.0) then
          if (betax.eq.(-2)) then
            write(9,100)zz,ex,ey,ez,tmpex,tmpey,tmpez !,bnseg
       ! if (bnseg.eq.-2) then
        ! write(10,100)zz,apx,apx1,apx2,apx3,apx4,2*pos(1)*a01kk   !a10kk*kk*I0p(kk*rsq)*cos1t*cos(ki*(zzm-Lr))
         ! kka=kk*rsq
         !write(33,100)kka,I2(kka),I4p(kka),I6p(kka),I0p(kka),I2p(kka)
          endif
        endif

100     format(8(1x,e12.6),/)

        end subroutine getfld_RFQ


        double precision FUNCTION I0(x)
        implicit none
        double precision::i6,i5,i4,i3,i2,i1
        double precision,intent(in)::x
        if (x.eq.0) then
        i0=1
        else
        i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
        i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
        i4=10*i5/x+i6
        i3=8*i4/x+i5
        i2=6*i3/x+i4
        i1=4*i2/x+i3
        i0=2*i1/x+i2
        endif
        end function I0

       double precision   FUNCTION I0p(x)
       implicit none
       double precision::i6,i5,i4,i3,i2,i1
       double precision,intent(in)::x
       if(x.eq.0) then
       i0p=1
       else
       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
       i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
       i4=10*i5/x+i6
       i3=8*i4/x+i5
       i2=6*i3/x+i4
       i0p=4*i2/x+i3
       endif
       end function I0p

       double precision FUNCTION I2p(x)
       implicit none
       double precision::i6,i5,i4,i3,i2,i1
       double precision,intent(in)::x
       double precision:: elap
       if(x.eq.0) then
       i2p=1
       else
       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
       i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
       i4=10*i5/x+i6
       i3=8*i4/x+i5
       i2=6*i3/x+i4
       i1=4*i2/x+i3
       I2p=i1-2*i2/x
       endif
       end function I2p

       double precision  FUNCTION I4p(x)
       implicit none
       double precision::i6,i5,i4,i3,i2,i1
       double precision,intent(in)::x
       double precision:: elap
       if(x.eq.0) then
       i4p=1
       else
       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
       i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
       i4=10*i5/x+i6
       i3=8*i4/x+i5
       i2=6*i3/x+i4
       i1=4*i2/x+i3
       I4p=i3-4*i4/x
       endif
       end function I4p

       double precision  FUNCTION I6p(x)
       implicit none
       double precision::i6,i5,i4,i3,i2,i1
       double precision,intent(in)::x
       if(x.eq.0)then
       i6p=1
       else
       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
       i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
       i6p=i5-6*i6/x

       endif
       end function I6p

       double precision  FUNCTION I6(x)
       implicit none
       double precision,intent(in)::x

       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720

       end function I6

       double precision FUNCTION I4(x)
       implicit none
       double precision::i6,i5
       double precision,intent(in)::x
       if(x.eq.0) then
       i4=1
       else
       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
       i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
       i4=10*i5/x+i6
       endif
       end function I4

       double precision FUNCTION I2(x)
       implicit none
       double precision::i6,i5,i4,i3
       double precision,intent(in)::x
       if (x.eq.0) then
       i2=1
       else
       i6=(X/2)**6*(1+(X/2)**2/7+(X/2)**4/112)/720
       i5=(x/2)**5*(1+(x/2)**2/6+(x/2)**4/84)/120
       i4=10*i5/x+i6
       i3=8*i4/x+i5
       i2=6*i3/x+i4
       endif
       end function I2

       end module RFQclass
