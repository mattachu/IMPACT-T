!----------------------------------------------------------------
! (c) Copyright, 2018 by the Regents of the University of California.
! Fieldmapclass: 3D field map class in Field module of APPLICATION
!                 layer.
! Version: 2.0
! Author: Ji Qiang
! Description: This class reads in a field map
! Comments:
!----------------------------------------------------------------
    module Fieldmapclass
     use PhysConstclass
     use Dataclass
     integer, private, parameter :: Nparam = 24
     type Fieldmap
          !Itype=202
          integer :: Nseg,Mapstp,Itype
          double precision :: Length
          double precision, dimension(Nparam) :: Param
     end type Fieldmap
     !Param(1):zmin/edge
     !!!Param(2):Frequency
     !Param(2): Field ID
     !Param(3):Input phase
     !Param(4):xmin
     !Param(5):xmax
     !Param(6):ymin
     !Param(7):ymax
     !Param(8):Nx
     !Param(9):Ny
     !Param(10):Nz
     !Param(11):ke
     !Param(12):kb
     !Param(13):Aperture
     interface getparam_Fieldmap
        module procedure getparam1_Fieldmap, &
                                      getparam2_Fieldmap, &
                                      getparam3_Fieldmap
     end interface
     interface setparam_Fieldmap
        module procedure setparam1_Fieldmap, &
                                      setparam2_Fieldmap, &
                                      setparam3_Fieldmap
     end interface

     interface setparam_fielddata
        module procedure setparam1_fielddata, &
                                        setparam2_fielddata
     end interface
    contains

        subroutine construct_Fieldmap(this,numseg,nmpstp,type,blength)
        implicit none
        type(Fieldmap),intent(out) :: this
        integer,intent(in) :: numseg,nmpstp, type
        double precision, intent(in) :: blength
        this%Nseg = numseg
        this%Mapstp = nmpstp
        this%Itype = type
        this%Length = blength

        this%Param = 0.0
        end subroutine construct_Fieldmap

        subroutine setparam1_Fieldmap(this,i,value)
          implicit none
          type(Fieldmap),intent(inout) :: this
           integer, intent(in) :: i
           double precision, intent(in) :: value

           this%Param(i) =value

        end subroutine setparam1_Fieldmap

        subroutine setparam2_Fieldmap(this, values)
         implicit none
          type(Fieldmap),intent(inout) :: this
           double precision, dimension(:), intent(in) :: values

           this%Param=values

        end subroutine setparam2_Fieldmap

        subroutine setparam3_Fieldmap(this,numseg,nmpstp, type, blength)
        implicit none
          type(Fieldmap),intent(inout) :: this
           integer, intent(in) :: numseg,nmpstp,type
           double precision, intent(in) :: blength

           this%Nseg = numseg
           this%Mapstp =nmpstp
           this%Itype = type
           this%Length = blength

        end subroutine setparam3_Fieldmap

        subroutine getparam1_Fieldmap(this,i,blparam)
         implicit none
        type (Fieldmap), intent(in) :: this
        integer, intent(in) :: i
        double precision, intent(out) :: blparam

        blparam = this%Param(i)

        end subroutine getparam1_Fieldmap

        subroutine getparam2_Fieldmap(this, blparams)
        implicit none
        type (Fieldmap), intent(in) :: this
        double precision, dimension(:), intent(out) :: blparams

        blparams = this%Param

        end subroutine getparam2_Fieldmap

        subroutine getparam3_Fieldmap(this,blength, bnseg, bmapstp, btype)
          implicit none
        type (Fieldmap), intent(in) :: this
        double precision, intent(out) :: blength
        integer, intent(out) :: bnseg,bmapstp,btype

        blength = this%Length
        bnseg = this%Nseg
        bmapstp = this%Mapstp
        btype = this%Itype

        end subroutine getparam3_Fieldmap

        subroutine setparam1_fielddata(this,n)
        implicit none
        type(fielddata),intent(inout) :: this
        !character, intent(in) :: n
        double precision, intent(in) :: n
        integer :: row, xcol, ycol, zcol
        integer :: Nx, Ny, Nz
        integer :: linenum
        double precision :: xstep, ystep, zstep, xmin, ymin, zmin ,xmax,ymax,zmax
        double precision :: xc, yc, zc, ex, ey, ez, bx, by, bz
        !double precision :: xstep, ystep, zstep
        character(len=4) :: FieldID
        !FieldID=(int(n))
        write(FieldID,'(i3.3)') int(n)
        print*, "reading fieldmap: " ,FieldID
    !    print*,  Trim(FieldID)
         zmin=this%ZminRfgt
         zmax=this%ZmaxRfgt
         xmin=this%XminRfgt
         xmax=this%XmaxRfgt
         ymin=this%YminRfgt
         ymax=this%YmaxRfgt
         Nx=this%NxIntvRfgt
         Ny=this%NyIntvRfgt
         Nz=this%NzIntvRfgt
         xstep=(xmax-xmin)/Nx
        ystep=(ymax-ymin)/Ny
        zstep=(zmax-zmin)/Nz

        row=0
        linenum = 0
        open(5,file="E"//Trim(FieldID)//".txt",status="old")            !read data
10      continue
            read(5,*,end=100) ex,ey, ez, bx, by, bz
            linenum = linenum +1

        !    xc=xc/1000.0
        !    yc=yc/1000.0
        !    zc=zc/1000.0
            bx=bx/10000000.0
            by=by/10000000.0
            bz=bz/10000000.0

        !    row=row+1
        !    xcol=anint((xc-xmin)/xstep)+1
        !    ycol=anint((yc-ymin)/ystep)+1
        !    zcol=anint((zc-zmin)/zstep)+1
            zcol=(linenum-1)/((Nx+1)*(Ny+1))+1
            ycol=(linenum-1-(zcol-1)*(Nx+1)*(Ny+1))/(Nx+1)+1
            xcol=linenum-(zcol-1)*(Nx+1)*(Ny+1)-(ycol-1)*(Nx+1)
            this%Exgridt(xcol,ycol,zcol)=ex
            this%Eygridt(xcol,ycol,zcol)=ey
           ! if (zcol.le.15) then
            !    this%Ezgridt(xcol,ycol,zcol)=0.0
            !else
            this%Ezgridt(xcol,ycol,zcol)=ez
            !endif
            this%Bxgridt(xcol,ycol,zcol)=bx
            this%Bygridt(xcol,ycol,zcol)=by
            this%Bzgridt(xcol,ycol,zcol)=bz
         goto 10
100      continue
         close(5)
        end subroutine setparam1_fielddata

        subroutine setparam2_fielddata(this,fldmp,blength)
        implicit none
        type(fielddata),intent(out) :: this
        type(Fieldmap), intent(in) :: fldmp
        double precision, intent(in) :: blength

        this%NxIntvRfgt=fldmp%param(8)
        this%NyIntvRfgt=fldmp%param(9)
        this%NzIntvRfgt=fldmp%param(10)
        this%XminRfgt=fldmp%param(4)
        this%YminRfgt=fldmp%param(6)
        this%ZminRfgt=fldmp%param(1)
        this%XmaxRfgt=fldmp%param(5)
        this%YmaxRfgt=fldmp%param(7)
        this%ZmaxRfgt=fldmp%param(1)+blength


        allocate(this%Exgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Eygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Ezgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bxgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bygridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))
        allocate(this%Bzgridt(this%NxIntvRfgt+1,this%NyIntvRfgt+1,&
                 this%NzIntvRfgt+1))

        end subroutine setparam2_fielddata

        subroutine getfld_Fieldmap (pos,momentum,extfld,this,fldata)
        implicit none

        double precision, dimension(4), intent(in) :: pos
        double precision,dimension(3), intent(in) :: momentum
        type (Fieldmap),intent(in) :: this
        type (fielddata), intent(in) :: fldata
        double precision, dimension(6), intent(out) :: extfld

        double precision :: xmin, ymin, zmin, xmax, ymax, zmax, xstep, ystep, zstep
        double precision :: xc,yc,zc,tt, betax, betay, betaz
        double precision :: dx1, dy1, dz1,dx2,dy2,dz2
        double precision :: ke, kb
        double precision :: fai0, fais, ww,freq
     !   double precision :: Ex1,Ex2,Ex3, Ex4, Ex5, Ex6, Ex7, Ex8
    !  double precision :: Ey1,Ey2,Ey3, Ey4, Ey5, Ey6, Ey7, Ey8
     !   double precision :: Ez1,Ez2,Ez3, Ez4, Ez5, Ez6, Ez7, Ez8
      !  double precision :: Bx1, Bx2, Bx3, Bx4, Bx5, Bx6, Bx7, Bx8
      !  double precision :: By1, By2, By3, By4, By5, By6, By7, By8
       ! double precision :: Bz1, Bz2, Bz3, Bz4, Bz5, Bz6, Bz7, Bz8
       double precision :: Ex,Ey,Ez,Bx,By,Bz
        double precision :: tmpex,tmpey,tmpez, tmpbx, tmpby, tmpbz
        integer :: xcol, ycol, zcol    !, xcol2,ycol2,zcol2
        integer :: i, Nx, Ny, Nz


         zmin=fldata%ZminRfgt
         zmax=fldata%ZmaxRfgt
         xmin=fldata%XminRfgt
         xmax=fldata%XmaxRfgt
         ymin=fldata%YminRfgt
         ymax=fldata%YmaxRfgt
         Nx=fldata%NxIntvRfgt
         Ny=fldata%NyIntvRfgt
         Nz=fldata%NzIntvRfgt

         !freq= this%Param(2)
         fai0=this%Param(3)/ Rad2deg
         ke=this%Param(11)
         kb=this%Param(12)

        xstep=(xmax-xmin)/Nx
        ystep=(ymax-ymin)/Ny
        zstep=(zmax-zmin)/Nz

        xc=pos(1)
        yc=pos(2)
        zc=pos(3)
        tt=pos(4)
        betax=momentum(1)
        betay=momentum(2)
        betaz=momentum(3)

        ww=2*Pi*scfreq
        fais=fai0+ww*tt

         extfld(1)=0
         extfld(2)=0
         extfld(3)=0
         extfld(4)=0
         extfld(5)=0
         extfld(6)=0

        if ((zc.le.zmin).or.(zc.ge.zmax)) then
            Ex=0.0
            Ey=0.0
            Ez=0.0
            Bx=0.0
            By=0.0
            Bz=0.0
        else

        xcol=floor((xc-xmin)/xstep)+1
        ycol=floor((yc-ymin)/ystep)+1
        zcol=floor((zc-zmin)/zstep)+1
       ! xcol2=xcol+1
       ! ycol2=ycol+1
       ! zcol2=zcol+1
        if(xcol.le.0 .or. ycol.le.0 .or. xcol.gt.Nx .or. ycol.gt.Ny .or. zcol.le.0 .or. zcol.gt.Nz) then
            Ex=0.0
            Ey=0.0
            Ez=0.0
            Bx=0.0
            By=0.0
            Bz=0.0
        else
        dx1=xc-xstep*(xcol-1)-xmin
        dy1=yc-ystep*(ycol-1)-ymin
        dz1=zc-zstep*(zcol-1)-zmin
        !dx2=xstep-dx1
        !dy2=ystep-dy1
        !dz2=zstep-dz1

        !Interpolation(this,xx,yy,zz,xcol,ycol,zcol,dx1,dy1,dz1,dx,dy,dz,ff)

        call Interpolation(fldata%Exgridt,xc,yc,zc,xcol,ycol,zcol,dx1,dy1,dz1,xstep,ystep,zstep,Ex)
        call Interpolation(fldata%Eygridt,xc,yc,zc,xcol,ycol,zcol,dx1,dy1,dz1,xstep,ystep,zstep,Ey)
        call Interpolation(fldata%Ezgridt,xc,yc,zc,xcol,ycol,zcol,dx1,dy1,dz1,xstep,ystep,zstep,Ez)
        call Interpolation(fldata%Bxgridt,xc,yc,zc,xcol,ycol,zcol,dx1,dy1,dz1,xstep,ystep,zstep,Bx)
        call Interpolation(fldata%Bygridt,xc,yc,zc,xcol,ycol,zcol,dx1,dy1,dz1,xstep,ystep,zstep,By)
        call Interpolation(fldata%Bzgridt,xc,yc,zc,xcol,ycol,zcol,dx1,dy1,dz1,xstep,ystep,zstep,Bz)

        endif
        endif

        tmpex=Ex*cos(fais)*ke
        tmpey=Ey*cos(fais)*ke
        tmpez=Ez*cos(fais)*ke
        tmpbx=Bx*sin(fais)*kb
        tmpby=By*sin(fais)*kb
        tmpbz=Bz*sin(fais)*kb

        extfld(1)=extfld(1)+tmpex
        extfld(2)=extfld(2)+tmpey
        extfld(3)=extfld(3)+tmpez
        extfld(4)=extfld(4)+tmpbx
        extfld(5)=extfld(5)+tmpby
        extfld(6)=extfld(6)+tmpbz

   !    write(3,120)xc,yc,zc,betax,betay,betaz,Ex,Ey,Ez,tmpex,cos(fais),tmpez
!120     format(13(1x,e12.6),/)
         if (betax.eq.(-2)) then
           write(9,110)zc,Ex*cos(fai0),Ey*cos(fai0),Ez*cos(fai0),tmpex,tmpey,tmpez !,betax,betay,betaz
110     format(8(1x,e12.6),/)
       endif
        end subroutine getfld_Fieldmap
       ! test field
        subroutine getField_fldmp(xx,yy,this)
        type (fielddata), intent(in) :: this
        double precision, intent(in) :: xx, yy

        double precision :: zmin,zmax,xmin,xmax,ymin, ymax
        double precision :: xstep, ystep, zstep
        integer:: Nx,Ny,Nz,xcol,xcol2,ycol,ycol2,zcol,zcol2,i
        double precision :: ex,ey,ez
        double precision :: dx, dy, dz, Ex1,Ex2,Ey1,Ey2,Ez1,Ez2
        double precision :: Bx1,Bx2,By1,By2,Bz1, Bz2
        double precision :: tmpex,tmpyey,tmpez
        double precision :: fai0,fais
        double precision :: ke, kb

        ke = 1.95
        kb = 1.95

        zmin=this%ZminRfgt
         zmax=this%ZmaxRfgt
         xmin=this%XminRfgt
         xmax=this%XmaxRfgt
         ymin=this%YminRfgt
         ymax=this%YmaxRfgt
         Nx=this%NxIntvRfgt
         Ny=this%NyIntvRfgt
         Nz=this%NzIntvRfgt
         xstep=(xmax-xmin)/Nx
        ystep=(ymax-ymin)/Ny
        zstep=(zmax-zmin)/Nz


        xcol=floor((xx-xmin)/xstep)+1
        ycol=floor((yy-ymin)/ystep)+1

        do i=1,Nz+1
        zcol=i
        xcol2=xcol+1
        ycol2=ycol+1
       ! zcol2=zcol+1
        zz=(i-1)*zstep+zmin
        dx=xx-(xstep*(xcol-1)+xmin)
        dy=yy-(ystep*(ycol-1)+ymin)
      !  dz=zc-zstep*zcol

        Ex1=this%Exgridt(xcol,ycol,zcol)
        Ex2=this%Exgridt(xcol2,ycol,zcol)
        Ey1=this%Eygridt(xcol,ycol,zcol)
        Ey2=this%Eygridt(xcol,ycol2,zcol)
       ! Ez1=fldata%Ezgridt(xcol,ycol,zcol)
       ! Ez2=fldata%Ezgridt(xcol2,ycol2,zcol2)
        Bx1=this%Bxgridt(xcol,ycol,zcol)
        Bx2=this%Bxgridt(xcol2,ycol,zcol)
        By1=this%Bygridt(xcol,ycol,zcol)
        By2=this%Bygridt(xcol,ycol2,zcol)

        Ex=(Ex2-Ex1)/xstep*dx+Ex1
        Ey=(Ey2-Ey1)/ystep*dy+Ey1
        Bx=(Bx2-Bx1)/xstep*dx+Bx1
        By=(By2-By1)/ystep*dy+By1
        Ez=this%Ezgridt(xcol,ycol,zcol)
        Bz=this%Bzgridt(xcol,ycol,zcol)
         write(2,130)zz,Ex*ke,Ey*ke,Ez*ke,Bx*kb,By*kb,Bz*kb
130     format(8(1x,e12.6),/)
        enddo

        end subroutine getField_fldmp

        subroutine Interpolation(this,xx,yy,zz,xcol,ycol,zcol,dx1,dy1,dz1,dx,dy,dz,ff)

        double complex, pointer,dimension(:,:,:),intent(in)::this
        double precision, intent(in):: xx,yy,zz
        double precision, intent(in):: dx1,dy1,dz1
        double precision, intent(in):: dx,dy,dz
        double precision,intent(out):: ff
        integer, intent(in) :: xcol,ycol,zcol

        integer :: xcol2,ycol2,zcol2
        double precision :: dx2, dy2, dz2
        double precision :: f1,f2,f3,f4,f5,f6,f7,f8
        double precision :: v1,v2,v3,v4,v5,v6,v7,v8
        double precision :: vv

        dx2=dx-dx1
        dy2=dy-dy1
        dz2=dz-dz1
        xcol2=xcol+1
        ycol2=ycol+1
        zcol2=zcol+1

        f1=this(xcol,ycol,zcol)
        f2=this(xcol2,ycol,zcol)
        f3=this(xcol2,ycol,zcol2)
        f4=this(xcol,ycol,zcol2)
        f5=this(xcol,ycol2,zcol)
        f6=this(xcol2,ycol2,zcol)
        f7=this(xcol2,ycol2,zcol2)
        f8=this(xcol,ycol2,zcol2)

        vv=dx*dy*dz
        v1=dx2*dy2*dz2
        v2=dx1*dy2*dz2
        v3=dx1*dy2*dz1
        v4=dx2*dy2*dz1
        v5=dx2*dy1*dz2
        v6=dx1*dy1*dz2
        v7=dx1*dy1*dz1
        v8=dx2*dy1*dz1

        ff=(v1*f1+v2*f2+v3*f3+v4*f4+v5*f5+v6*f6+v7*f7+v8*f8)/vv

    end subroutine Interpolation



    end module Fieldmapclass
