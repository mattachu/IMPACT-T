!----------------------------------------------------------------
! (c) Copyright, 2017 by the Regents of the University of California.
! NumConstclass: numerical constants class in CONSTANTS module of
!                 DATA STRUCTURE layer.
!
! MODULE  : ... NumConstclass
! VERSION : ... 1.0
!> @author
!> Ji Qiang
!
! DESCRIPTION: 
!> This class defines the maximum size for numerical parameters
!> in the simulation.
! Comments: 
!----------------------------------------------------------------
      module NumConstclass
        implicit none
      
        !numerical parameters --------------------------
        !> phase dimension
        integer, parameter :: Pdim = 6 
        !> maximum \# of particles on a single processor.
        integer, parameter :: Nplcmax = 2000000
        !> @name
        
        !> @{
        !> maximum \# of grids in x,y,z dimension on a single processor.
        integer, parameter :: Nxlcmax = 1280
        integer, parameter :: Nylcmax = 1280
        integer, parameter :: Nzlcmax = 1280
        !> @}
        !> maximum \# of beam line elements
        integer, parameter :: Nblemtmax = 3000
        !> maximum \# of drift space
        integer, parameter :: Ndriftmax = 1400
        !> maximum \# of quadrupoles
        integer, parameter :: Nquadmax = 400
        integer, parameter :: Ncfmax = 100
        !> maximum \# of dipoles
        integer, parameter :: Ndipolemax = 100
        !> maximum # of RFQ cells
        integer, parameter :: Nrfqmax = 10000
        integer, parameter :: Nfldmpmax = 1000
        !> maximum # of rf gaps
        integer, parameter :: Nrfgapmax = 1000
        integer, parameter :: Ncclmax = 1000
        integer, parameter :: Nccdtlmax = 1000
        integer, parameter :: Ndtlmax = 1000
        integer, parameter :: Nscmax = 1000
        !> maximum \# of beam position monitors
        integer, parameter :: Nbpmmax = 2000
        !> maximum \# of magnetic solenoid
        integer, parameter :: Nslmax = 100
        integer, parameter :: Nsolmax = 1000
        !> maximum \# of magnetic solenoid with RF field
        integer, parameter :: Nslrfmax = 100
        !> maximum \# of field storage array size
        integer, parameter :: Maxiifile = 1000
        !> maximum \# of bunches/bins allowed
        integer, parameter :: Nbunchmax = 50
        !> maximum \# of field overlap allowed
        integer, parameter :: Maxoverlap = 1000
      end module NumConstclass
