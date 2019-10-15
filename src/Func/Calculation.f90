!------------------------------------------------------------------------------!
! (c) Copyright, 2018 by the Regents of the University of California.          !
! CalculationClass: Class to handle simple calculations used in multiple       !
!                    modules; part of the DATA STRUCTURES layer.               !
! Version:     for test, based on v2.0                                         !
! Created:     Matt Easton, Peking University, 2017-11-21                      !
! Modified:    Matt Easton, Peking University, 2018-06-08                      !
! Description: This class contains FUNCTIONs and procedures for carrying out   !
!                various calculations that are useful across modules.          !
!-------------------------------------------------------------------------------

MODULE CalculationClass

    ! Includes

    ! Declarations
    IMPLICIT NONE
    PRIVATE

    ! Parameters
    LOGICAL,          PARAMETER, PRIVATE :: debug_mode = .FALSE.  ! Include debugging output?
    LOGICAL,          PARAMETER, PRIVATE :: detailed_debug_mode = .FALSE.  ! Warning: SLOW!
    DOUBLE PRECISION, PARAMETER, PRIVATE :: almost_zero = 0.1d-10

    ! Variables

    ! Types

    ! Interfaces
    INTERFACE is_zero
        MODULE PROCEDURE is_zero_real
        MODULE PROCEDURE is_zero_dp
        MODULE PROCEDURE is_zero_int
        MODULE PROCEDURE is_zero_log
        MODULE PROCEDURE is_zero_char
    END INTERFACE

    INTERFACE remove_duplicates
        MODULE PROCEDURE remove_duplicates_int
    END INTERFACE

    ! Public procedures
    PUBLIC :: is_zero
    PUBLIC :: remove_duplicates

CONTAINS


    ! FUNCTION to return whether an integer is zero
    !  (trivial, included for completeness and compatibility)
    LOGICAL FUNCTION is_zero_int(compare)
        ! Parameters
        INTEGER, INTENT(IN) :: compare

        is_zero_int = (compare == 0)

    END FUNCTION is_zero_int



    ! FUNCTION to return whether a real number is zero
    !  (real numbers are often not exactly zero)
    LOGICAL FUNCTION is_zero_real(compare)
        ! Parameters
        REAL, INTENT(IN) :: compare

        is_zero_real = (compare < almost_zero)

    END FUNCTION is_zero_real



    ! FUNCTION to return whether a DOUBLE PRECISION number is zero
    !  (real numbers are often not exactly zero)
    LOGICAL FUNCTION is_zero_dp(compare)
        ! Parameters
        DOUBLE PRECISION, INTENT(IN) :: compare

        is_zero_dp = (compare < almost_zero)

    END FUNCTION is_zero_dp



    ! FUNCTION to return whether a LOGICAL variable is zero
    !  (assume .FALSE. means zero, mainly included for completeness)
    LOGICAL FUNCTION is_zero_log(compare)
        ! Parameters
        LOGICAL, INTENT(IN) :: compare

        is_zero_log = (.NOT.compare)

    END FUNCTION is_zero_log



    ! FUNCTION to return whether a character variable is zero
    !  (assume '0' means zero, mainly included for completeness)
    LOGICAL FUNCTION is_zero_char(compare)
        ! Parameters
        CHARACTER(LEN=*), INTENT(IN) :: compare

        is_zero_char = (TRIM(compare) == '0')

    END FUNCTION is_zero_char



    ! Remove duplicates from array_in and return (shorter) array_out
    SUBROUTINE remove_duplicates_int(array_in, array_out)
      ! Parameters
      INTEGER, ALLOCATABLE, INTENT(IN)  :: array_in(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: array_out(:)
      ! Variables
      INTEGER, ALLOCATABLE :: array_temp(:)
      INTEGER :: i, count

      ! Temporary array stores results but is incorrectly sized
      ALLOCATE(array_temp(SIZE(array_in)))
      array_temp = 0

      ! Loop through all values in array_in
      count = 1
      array_temp(1) = array_in(1)
      DO i = 2, SIZE(array_in)
        ! If the number already exists in array_temp then skip to next i
        IF(ANY(array_temp == array_in(i))) CYCLE
        ! No match found so add it to the output
        count = count + 1
        array_temp(count) = array_in(i)
      ENDDO

      ! Create output array
      ALLOCATE(array_out(count))
      array_out = array_temp(1:count)

    END SUBROUTINE remove_duplicates_int



END MODULE CalculationClass
