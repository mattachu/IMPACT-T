!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! Inputclass: Input class in I/O module of CONTROL layer. 
! 
! MODULE  : ... Inputclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
!
! DESCRIPTION: 
!> This class defines functions to input the global
!> beam and computational parameters and the lattice input
!> parameters in the accelerator.
! Comments: J.Q modified the source code so that the user can put comment
!           lines starting with "!" for each number line in the input file
!           "ImpactT.in".
!----------------------------------------------------------------
      module Inputclass
        !use mpistub

        ! Parameters for interactions input file
        CHARACTER(15), PARAMETER :: INTERACTION_INPUT_FILE = "Interactions.in"
        INTEGER :: INTERACTION_FILEUNIT ! set when file is opened

        interface in_Input
          module procedure in1_Input, in2_Input, in3_Input
        end interface
      contains
        !> Start MPI
        subroutine init_Input(time)
        implicit none
        include 'mpif.h'
        double precision, intent(out) :: time
        integer :: ierr

        ! start MPI library.
!        call MPI_INIT(ierr)
        time = MPI_WTIME()
        ! for measurement of memory usage.
        !call system_stats()

        end subroutine init_Input
 
        !> Input all parameters except beam line element parameters.
        subroutine in1_Input(odim,onp,onx,ony,onz,oflagbc,oflagdist, &
        orstartflg,oflagmap,distparam,nparam,obcurr,obkenergy,obmass,&
        obcharge,obfreq,oxrad,oyrad,operdlen,onblem,onpcol,onprow,oflagerr,&
        oflagdiag,oflagsbstp,ophsini,odt,ontstep,onbunch,oflagimg,&
        onemission,otemission,ozimage,oflaginteract)

        implicit none
        include 'mpif.h'
        integer, intent(out) :: odim,onp,onx,ony,onz,oflagbc,oflagdist
        integer, intent(out) :: orstartflg,oflagmap,onblem,onpcol,onprow 
        integer, intent(out) :: oflagerr,oflagdiag,oflagsbstp,ontstep,&
                                onbunch,oflagimg,onemission,oflaginteract
        integer, intent(in) :: nparam
        double precision, dimension(nparam), intent(out) :: distparam
        double precision, intent(out) :: obcurr,obkenergy,obmass,odt
        double precision, intent(out) :: obcharge,obfreq,operdlen,&
                                         oxrad,oyrad,ophsini,&
                                         otemission,ozimage
        double precision :: xjunk
        integer :: my_rank,nproc,ierr,np,itot,njunk1,njunk2,njunk3
        character*1 comst
        integer :: ii,jj,i

        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

!        onpcol = 32
!        onprow = 16

        if(my_rank.eq.0) then
!          print*,"input the column and row number of processors:"
!          read(*,*)onpcol,onprow

          open(unit=13,file='ImpactT.in',status='old')
!          open(unit=24,file='fort.24')
!          open(unit=25,file='fort.25')
!          open(unit=26,file='fort.26')

          ii = 0
          jj = 0
10        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 10
          else
            backspace(13,err=789)
            read(13,*)onpcol,onprow
            ii = ii+1
          endif
20        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 20
          else
            backspace(13,err=789)
            read(13,*)odt,ontstep,onbunch
            ii = ii+1
          endif
30        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 30
          else
            backspace(13,err=789)
            read(13,*)odim,onp,oflagmap,oflagerr,oflagdiag,oflagimg,ozimage 
            CALL readbonusflags_Input(13, oflaginteract)
            ii = ii+1
          endif
40        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 40
          else
            backspace(13,err=789)
            read(13,*)onx,ony,onz,oflagbc,oxrad,oyrad,operdlen
            ii = ii+1
          endif
50        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 50
          else
            backspace(13,err=789)
            read(13,*)oflagdist,orstartflg,oflagsbstp,onemission,otemission
            ii = ii+1
          endif

          distparam = 0.0
80        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 80
          else
            backspace(13,err=789)
            read(13,*)distparam(1),distparam(2),distparam(3),distparam(4),&
                    distparam(5),distparam(6),distparam(7)
            ii = ii+1
          endif
90        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 90
          else
            backspace(13,err=789)
            read(13,*)distparam(8),distparam(9),distparam(10),&
                    distparam(11),distparam(12),distparam(13),distparam(14)
            ii = ii+1
          endif
101       continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 101
          else
            backspace(13,err=789)
            read(13,*)distparam(15),distparam(16),&
                    distparam(17),distparam(18),distparam(19),distparam(20),distparam(21)
            ii = ii+1
          endif
102       continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 102
          else
            backspace(13,err=789)
            read(13,*)obcurr,obkenergy,obmass,obcharge,obfreq,ophsini
            ii = ii+1
          endif
!          print*,"# of lines before beam line elements: ",ii,jj

!          read(13,*)onblem

          !count the # of beam line elements.
          itot=0
          njunk3 = 0
123       continue
            read(13,*,end=789)comst
            if(comst.ne."!") then
              backspace(13,err=789)
              read(13,*,end=789)xjunk,njunk1,njunk2,njunk3
              itot = itot + 1
              if(njunk3.eq.-99)then
                goto 789
              endif
            endif
          goto 123
789       continue
          onblem=itot

!          write(6,*)'onblem = ',onblem
          rewind(13)

          do i = 1, jj
            read(13,*)comst
          enddo
        endif

        call MPI_BCAST(onpcol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onprow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ontstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(odim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagmap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagerr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagdiag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagimg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onemission,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onbunch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ony,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagbc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(odt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(otemission,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ozimage,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oxrad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(oyrad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(operdlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(oflagdist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(orstartflg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagsbstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onblem,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(distparam(1),nparam,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                       ierr)
        call MPI_BCAST(obcurr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obkenergy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obmass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obcharge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obfreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(ophsini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(oflaginteract,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

        end subroutine in1_Input

        !> Input beam line element parameters.
!        subroutine in2_Input(onblem,operd,oblength,obnseg,obmpstp,&
!                             obtype,value1,value2,value3,value4,value5)
        subroutine in2_Input(onblem,oblength,obnseg,obmpstp,&
        obtype,value0,value1,value2,value3,value4,value5,value6,&
        value7,value8,&
        value9,value10,value11,value12,value13,value14,value15,value16,&
        value17,value18,value19,value20,value21,value22,value23,value24,value25)
        implicit none
        include 'mpif.h'
!        integer,intent(in) :: onblem,operd
        integer,intent(in) :: onblem
        integer,intent(out) :: obnseg(onblem)
        integer,intent(out) :: obmpstp(onblem)
        integer,intent(out) :: obtype(onblem)
        double precision,intent(out) :: oblength(onblem)
        double precision,dimension(onblem),intent(out) :: value0,&
        value1,value2,&
        value3,value4,value5,value6,value7,value8,value9,value10,value11,&
        value12,value13,value14,value15,value16,value17,value18,value19,&
        value20,value21,value22,value23,value24,value25
        integer :: i,irf
        integer :: myrank,ierr
        character*1 comst

        call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

        value0 = 0.0
        value1 = 0.0
        value2 = 0.0
        value3 = 0.0
        value4 = 0.0
        value5 = 0.0
        value6 = 0.0
        value7 = 0.0
        value8 = 0.0
        value9 = 0.0
        value10 = 0.0
        value11 = 0.0
        value12 = 0.0
        value13 = 0.0
        value14 = 0.0
        value15 = 0.0
        value16 = 0.0
        value17 = 0.0
        value18 = 0.0
        value19 = 0.0
        value20 = 0.0
        value21 = 0.0
        value22 = 0.0
        value23 = 0.0
        value24 = 0.0
        value25 = 0.0

        if(myrank.eq.0) then

          i=0
123       continue
            read(13,*,end=789)comst
            if(comst.ne."!") then
              backspace(13,err=789)
              i = i + 1
              read(13,*)oblength(i),obnseg(i),obmpstp(i),obtype(i),&
              value0(i),value1(i),value2(i),value3(i),value4(i),value5(i),&
              value6(i),&
              value7(i),value8(i),value9(i),value10(i),value11(i),value12(i),&
              value13(i),value14(i),value15(i),value16(i),value17(i),value18(i),&
              value19(i),value20(i),value21(i),value22(i),value23(i),value24(i),value25(i)
              if(obtype(i).eq.-99)then
                goto 789
              endif
            endif
          goto 123
789       continue

          print*,"nblem: ",i,onblem

          close(13)
        endif

        call MPI_BCAST(oblength,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value0,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value1,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value2,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value3,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value4,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value5,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value6,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value7,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value8,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value9,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value10,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value11,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value12,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value13,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value14,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value15,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value16,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value17,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value18,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value19,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value20,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value21,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value22,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value23,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(value24,onblem,MPI_DOUBLE_PRECISION,0,&
                       MPI_COMM_WORLD,ierr)
        call MPI_BCAST(obnseg,onblem,MPI_INTEGER,0,MPI_COMM_WORLD,&
                       ierr)
        call MPI_BCAST(obmpstp,onblem,MPI_INTEGER,0,MPI_COMM_WORLD,&
                       ierr)
        call MPI_BCAST(obtype,onblem,MPI_INTEGER,0,MPI_COMM_WORLD,&
                       ierr)

        end subroutine in2_Input

        !> Input all parameters except beam line element parameters.
        subroutine in3_Input(odim,onp,onx,ony,onz,oflagbc,oflagdist, &
        orstartflg,oflagmap,distparam,nparam,obcurr,obkenergy,obmass,&
        obcharge,obfreq,oxrad,oyrad,operdlen,onblem,onpcol,onprow,oflagerr,&
        oflagdiag,oflagsbstp,ophsini,odt,ontstep,onbunch,oflagimg,ib)

        implicit none
        include 'mpif.h'
        integer, intent(out) :: odim,onp,onx,ony,onz,oflagbc,oflagdist
        integer, intent(out) :: orstartflg,oflagmap,onblem,onpcol,onprow 
        integer, intent(out) :: oflagerr,oflagdiag,oflagsbstp,ontstep,&
                                onbunch,oflagimg
        integer, intent(in) :: nparam,ib
        double precision, dimension(nparam), intent(out) :: distparam
        double precision, intent(out) :: obcurr,obkenergy,obmass,odt
        double precision, intent(out) :: obcharge,obfreq,operdlen,&
                                         oxrad,oyrad,ophsini
        double precision :: xjunk
        integer :: my_rank,nproc,ierr,np,itot,njunk1,njunk2,njunk3
        character*11 name1
        character*12 name2
        character*13 name3
        integer :: i,j,k,l
        integer :: onemission
        double precision :: otemission, ozimage
        character*1 comst
        integer :: ii,jj
 
        name1 = 'ImpactTx.in'
        name2 = 'ImpactTxx.in'
        name3 = 'ImpactTxxx.in'


        call MPI_COMM_RANK(MPI_COMM_WORLD,my_rank,ierr)
        call MPI_COMM_SIZE(MPI_COMM_WORLD,np,ierr)

!        onpcol = 32
!        onprow = 16

        if(my_rank.eq.0) then
          if(ib.lt.10) then
            name1(8:8) = char(ib+48)
            open(unit=13,file=name1,status='old')
          else if(ib.lt.100) then
            i = ib/10
            j = ib - 10*i
            name2(8:8) = char(i+48)
            name2(9:9) = char(j+48)
            open(unit=13,file=name2,status='old')
          else if(ib.lt.1000) then
            i = ib/100
            j = ib - 100*i
            k = j/10
            l = j - 10*k
            name3(8:8) = char(i+48)
            name3(9:9) = char(k+48)
            name3(10:10) = char(l+48)
            open(unit=13,file=name3,status='old')
          else
            print*,"over maximum # of input files:...."
            stop
          endif
          print*,"file name: ",name1,name2,name3

          ii = 0
          jj = 0
10        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 10
          else
            backspace(13,err=789)
            read(13,*)onpcol,onprow
            ii = ii+1
          endif
20        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 20
          else
            backspace(13,err=789)
            read(13,*)odt,ontstep,onbunch
            ii = ii+1
          endif
30        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 30
          else
            backspace(13,err=789)
            read(13,*)odim,onp,oflagmap,oflagerr,oflagdiag,oflagimg,ozimage 
            ii = ii+1
          endif
40        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 40
          else
            backspace(13,err=789)
            read(13,*)onx,ony,onz,oflagbc,oxrad,oyrad,operdlen
            ii = ii+1
          endif
50        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 50
          else
            backspace(13,err=789)
            read(13,*)oflagdist,orstartflg,oflagsbstp,onemission,otemission
            ii = ii+1
          endif

          distparam = 0.0
80        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 80
          else
            backspace(13,err=789)
            read(13,*)distparam(1),distparam(2),distparam(3),distparam(4),&
                    distparam(5),distparam(6),distparam(7)
            ii = ii+1
          endif
90        continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 90
          else
            backspace(13,err=789)
            read(13,*)distparam(8),distparam(9),distparam(10),&
                    distparam(11),distparam(12),distparam(13),distparam(14)
            ii = ii+1
          endif
101       continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 101
          else
            backspace(13,err=789)
            read(13,*)distparam(15),distparam(16),&
                    distparam(17),distparam(18),distparam(19),distparam(20),distparam(21)
            ii = ii+1
          endif
102       continue
          read(13,*)comst
          jj = jj + 1
          if(comst.eq."!") then
            goto 102
          else
            backspace(13,err=789)
            read(13,*)obcurr,obkenergy,obmass,obcharge,obfreq,ophsini
            ii = ii+1
          endif
          print*,"# of lines before beam line elements: ",ii,jj

!          read(13,*)onblem
          !count the # of beam line elements.
          itot=0
          njunk3 = 0
123       continue
            read(13,*,end=789)comst
            if(comst.ne."!") then
              backspace(13,err=789)
              read(13,*,end=789)xjunk,njunk1,njunk2,njunk3
              itot = itot + 1
              if(njunk3.eq.-99)then
                goto 789
              endif
            endif
          goto 123
  789     continue
          onblem=itot
          write(6,*)'onblem = ',onblem

          rewind(13)

          do i = 1, jj
            read(13,*)comst
          enddo

          close(13)
        endif

        call MPI_BCAST(onpcol,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onprow,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ontstep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(odim,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagmap,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagimg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagerr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagdiag,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onbunch,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(ony,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagbc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(odt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oxrad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(oyrad,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(operdlen,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(oflagdist,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(orstartflg,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(oflagsbstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(onblem,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(distparam(1),nparam,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                       ierr)
        call MPI_BCAST(obcurr,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obkenergy,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obmass,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obcharge,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(obfreq,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)
        call MPI_BCAST(ophsini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,&
                         ierr)

        end subroutine in3_Input



        ! Read in any bonus flags for additional modules that are
        ! not part of the standard input file
        !  - this allows adding extra module flags without needing to
        !    change all existing input files, assisting backwards compatibility
        !  - this is called *after* the standard reading of the line with the
        !    main flags on, so starts by backing up the file and reading the
        !    line again
        !  - at the end of the routine, the file should be ready to continue
        !    reading as normal
        SUBROUTINE readbonusflags_Input(fileunit, flag_interact)
          ! Parameters
          INTEGER, INTENT(IN)  :: fileunit
          INTEGER, INTENT(OUT) :: flag_interact
          ! Variables
          INTEGER          :: dim_in, np_in, flagmap_in, flagerr_in, &
                              flagdiag_in, flagimg_in
          INTEGER          :: dim_chk, np_chk, flagmap_chk
          DOUBLE PRECISION :: zimage_in
          INTEGER          :: ierr

          ! Read the last input line again with the bonus flags included
          BACKSPACE(fileunit)
          READ(fileunit, *, IOSTAT=ierr) dim_in, np_in, flagmap_in, flagerr_in, &
                                         flagdiag_in, flagimg_in, zimage_in, &
                                         flag_interact

          ! If we hit an error at this point, it's probably because
          !  there are no bonus flags
          IF(ierr /= 0) THEN
            flag_interact = 0
            BACKSPACE(fileunit)
            RETURN
          ENDIF

          ! Read in the first few standard flags again and compare
          !  - this is to check that we haven't moved to the next line
          BACKSPACE(fileunit)
          READ(fileunit, *, IOSTAT=ierr) dim_chk, np_chk, flagmap_chk
          IF((ierr /= 0).OR.(dim_chk /= dim_in).OR.(np_chk /= np_in) &
                        .OR.(flagmap_chk /= flagmap_in)) THEN
            flag_interact = 0
            BACKSPACE(fileunit)
            RETURN
          ENDIF

        END SUBROUTINE readbonusflags_Input



        ! Open the interactions input file and save the file unit number
        SUBROUTINE openinteractions_Input()

          OPEN(newunit=INTERACTION_FILEUNIT, file=INTERACTION_INPUT_FILE, &
               status='old')

        END SUBROUTINE openinteractions_Input



        ! Close the interactions input file
        SUBROUTINE closeinteractions_Input()

          CLOSE(unit=INTERACTION_FILEUNIT)

        END SUBROUTINE closeinteractions_Input



        ! Skip lines beginning with "!"
        SUBROUTINE skipcomments_Input(fileunit)
          ! Parameters
          INTEGER, INTENT(IN) :: fileunit
          ! Variables
          CHARACTER(1) :: first_character

          ! Loop until find a line that doesn't start with an !
          DO
            READ(fileunit,*) first_character
            IF(first_character /= "!") THEN
              ! Go back to start of line, ready to read data
              BACKSPACE(fileunit)
              EXIT
            ENDIF
          ENDDO

        END SUBROUTINE skipcomments_Input



        ! Read the general interactions settings
        SUBROUTINE readinteractions_Input(interaction_count, &
                                          gas_pressure, gas_temperature)
          ! Parameters
          INTEGER,          INTENT(OUT) :: interaction_count
          DOUBLE PRECISION, INTENT(OUT) :: gas_pressure, gas_temperature

          CALL skipcomments_Input(INTERACTION_FILEUNIT)
          READ(INTERACTION_FILEUNIT, *) interaction_count, &
                                        gas_pressure, gas_temperature

          ! Adjust units
          gas_temperature = gas_temperature + 273.15 ! K = deg C + 273.15

        END SUBROUTINE readinteractions_Input



        ! Read the settings for a single interaction
        SUBROUTINE readinteraction_Input(type, source_bunch, target_bunch, &
                                          interval, peak_energy, peak_cross)
          ! Parameters
          INTEGER,          INTENT(OUT) :: type, source_bunch, target_bunch
          DOUBLE PRECISION, INTENT(OUT) :: interval, peak_energy, peak_cross

          CALL skipcomments_Input(INTERACTION_FILEUNIT)
          READ(INTERACTION_FILEUNIT, *) type, source_bunch, target_bunch, &
                                        peak_energy, peak_cross, interval

        END SUBROUTINE readinteraction_Input



        ! Read in how many interaction bunches are to be created
        SUBROUTINE readinteractionbunchcount_Input(bunch_count)
          ! Parameters
          INTEGER, INTENT(OUT) :: bunch_count

          CALL skipcomments_Input(INTERACTION_FILEUNIT)
          READ(INTERACTION_FILEUNIT, *) bunch_count

        END SUBROUTINE readinteractionbunchcount_Input



        ! Read in the settings for a single interaction bunch
        SUBROUTINE readinteractionbunch_Input(current_in, energy_in, mass_in, &
                                              charge_in, phase_in, &
                                              npt_start, npt_max, npt_factor)
          ! Parameters
          DOUBLE PRECISION, INTENT(OUT) :: current_in, energy_in, mass_in
          DOUBLE PRECISION, INTENT(OUT) :: charge_in, phase_in
          INTEGER,          INTENT(OUT) :: npt_start, npt_max
          REAL,             INTENT(OUT) :: npt_factor
          ! Variables
          INTEGER :: bunch_id

          CALL skipcomments_Input(INTERACTION_FILEUNIT)
          READ(INTERACTION_FILEUNIT, *) bunch_id, current_in, energy_in, &
                                        mass_in, charge_in, phase_in, &
                                        npt_start, npt_max, npt_factor

        END SUBROUTINE readinteractionbunch_Input



      end module Inputclass

