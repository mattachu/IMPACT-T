!------------------------------------------------------------------------------!
! (c) Copyright, 2018 by the Regents of the University of California.          !
! QuickSort: Module to sort the entries in an array;                           !
!              part of the FUNCTIONS layer.                                    !
! Version:     for test, based on v2.0                                         !
! Created:     Matt Easton, Peking University, 2017-12-22                      !
! Modified:    Matt Easton, Peking University, 2018-01-31                      !
! Description: This module is almost entirely taken from the FLIBS library.    !
!-------------------------------------------------------------------------------
!
! Notes from qsortarray_template from FLIBS v0.9:
!
! qsortarray_template.f90 --
!     QuickSort for arrays:
!
!     This template has to be preprocessed by defining the
!     two following macros :
!     - "_QSORTARRAY_TYPE" is the name of the derived type
!     - "_QSORTARRAY_MODULE" is the module defining the
!       derived type
!
!     Note:
!     The data that need to be sorted are stored in
!     an array. A comparison function is used to enable
!     sorting wrt different criteria.
!     The data have the generic type _QSORTARRAY_TYPE.
!
!     Note:
!     Algorithm translated from Kernighan and Pike,
!     The Practice of Programming.
!
!     $Id: qsortarray_template.f90,v 1.1 2008/04/17 14:39:08 relaxmike Exp $
!
! Routines and functions specific to dictionaries
!
! qsort_array --
!     Use the QuickSort algorithm to sort an array
! Arguments:
!     array      Array to be sorted
!     compare    Comparison function
!
MODULE QuickSort

    ! Declarations
    IMPLICIT NONE
    PRIVATE

    ! Types
    TYPE QuickSortArray
        INTEGER :: key
    END TYPE

    ! Interfaces
    INTERFACE qsort
        MODULE PROCEDURE qsort_int
        ! TODO: other types
    END INTERFACE qsort

    ! Public procedures
    PUBLIC :: qsort

CONTAINS

    ! Handler for integer arrays
    SUBROUTINE qsort_int(array)
        ! Parameters
        INTEGER, INTENT(INOUT) :: array(:)
        ! Variables
        TYPE(QuickSortArray), ALLOCATABLE :: internal_array(:)

        ! Allocate memory
        ALLOCATE(internal_array(SIZE(array)))

        ! Copy array
        internal_array%key = array

        ! Call the quicksort routine
        CALL qsort_array(internal_array)

        ! Return results
        array = internal_array%key

    END SUBROUTINE qsort_int

    ! Comparison of two items
    LOGICAL FUNCTION qsort_compare(a,b)
        ! Parameters
        TYPE(QuickSortArray), INTENT(IN) :: a, b

        qsort_compare = (a%key < b%key)

    END FUNCTION qsort_compare

    ! Wrapper
    SUBROUTINE qsort_array(array)
        ! Parameters
        TYPE(QuickSortArray), INTENT(INOUT) :: array(:)
        ! Variables
        TYPE(QuickSortArray), ALLOCATABLE   :: backup(:)
        INTEGER, ALLOCATABLE                :: order(:)
        INTEGER                             :: i

        ! Size of arrays
        ALLOCATE(backup(1:SIZE(array)))
        ALLOCATE(order(1:SIZE(array)))

        ! Initial ordering
        DO i = 1,SIZE(order)
            order(i) = i
        ENDDO

        ! Run main sort algorithm
        CALL qsort_sort(array, order, 1, SIZE(array))

        ! Copy newly sorted data to main array
        DO i = 1,SIZE(order)
            backup(i) = array(order(i))
        ENDDO
        array = backup

        ! Tidy up
        DEALLOCATE(backup)
        DEALLOCATE(order)

    END SUBROUTINE qsort_array

    ! qsort_sort --
    !     Sort the array according to the QuickSort algorithm
    ! Arguments:
    !     array      Array of data
    !     order      Array holding the order
    !     left       Start index
    !     right      End index
    !     compare    Comparison function
    !
    RECURSIVE SUBROUTINE qsort_sort(array, order, left, right)
        ! Parameters
        TYPE(QuickSortArray), INTENT(INOUT) :: array(:)
        INTEGER,              INTENT(INOUT) :: order(:)
        INTEGER, INTENT(IN)                 :: left
        INTEGER, INTENT(IN)                 :: right
        ! Variables
        INTEGER :: i, last

        IF (left >= right) RETURN
        CALL qsort_swap(order, left, qsort_rand(left,right))
        last = left
        DO i = left+1, right
            IF(qsort_compare(array(order(i)), array(order(left)))) THEN
                last = last + 1
                CALL qsort_swap(order, last, i)
            ENDIF
        ENDDO
        CALL qsort_swap(order, left, last)
        CALL qsort_sort(array, order, left, last-1)
        CALL qsort_sort(array, order, last+1, right)

    END SUBROUTINE qsort_sort

    ! qsort_swap --
    !     Swap two array elements
    ! Arguments:
    !     order      Array holding the order
    !     first      First index
    !     second     Second index
    !
    SUBROUTINE qsort_swap(order, first, second)
        ! Parameters
        INTEGER, INTENT(INOUT) :: order(:)
        INTEGER, INTENT(IN)    :: first, second
        ! Variables
        INTEGER                :: tmp

        ! Swap the values
        tmp           = order(first)
        order(first)  = order(second)
        order(second) = tmp

    END SUBROUTINE qsort_swap

    ! qsort_rand --
    !     Determine a random integer number
    ! Arguments:
    !     lower      Lowest value
    !     upper      Greatest value
    !
    INTEGER FUNCTION qsort_rand(lower, upper)
        ! Parameters
        INTEGER, INTENT(IN) :: lower, upper
        ! Variables
        REAL                :: r

        ! Get random number
        CALL random_number( r )
        ! Scale number to lower/upper bounds
        qsort_rand =  lower + NINT(r * (upper-lower))

    END FUNCTION qsort_rand

END MODULE QuickSort
