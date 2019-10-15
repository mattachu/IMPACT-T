!------------------------------------------------------------------------------!
! (c) Copyright, 2018 by the Regents of the University of California.          !
! InteractionClass: Class to handle interactions between particles as part of  !
!                    the Beam module of APPLICATION layer.                     !
! Version:     for test, based on v2.0                                         !
! Created:     Matt Easton, Peking University, 2017-11-09                      !
! Modified:    Matt Easton, Peking University, 2018-06-08                      !
! Description: This class defines interactions that can occur between          !
!              particles such as stripping, dissociation, collisions, etc.     !
!-------------------------------------------------------------------------------

MODULE InteractionClass

  ! Includes
  USE BeamBunchClass
  USE PhysConstClass
  USE CalculationClass
  USE InputClass
  USE QuickSort

  ! Declarations
  IMPLICIT NONE
  PRIVATE

  ! Parameters
  LOGICAL, PARAMETER, PRIVATE :: debug_mode = .FALSE.  ! Include debugging output?

  ! Enumeration of interaction types
  INTEGER, PARAMETER :: IT_DISSOCIATION = 1
  INTEGER, PARAMETER :: IT_CAPTURE = 2
  INTEGER, PARAMETER :: IT_CAPTURESPLIT = 3

  ! Types
  TYPE Interaction
    INTEGER :: type
    INTEGER :: source_bunch_id
    INTEGER :: target_bunch_id
    DOUBLE PRECISION :: peak_energy
    DOUBLE PRECISION :: peak_cross
    DOUBLE PRECISION :: interval
    DOUBLE PRECISION :: next
  END TYPE Interaction

  TYPE BeamBunchTemplate
    DOUBLE PRECISION :: current
    DOUBLE PRECISION :: energy
    DOUBLE PRECISION :: mass
    DOUBLE PRECISION :: charge
    DOUBLE PRECISION :: phase
    INTEGER          :: npt_start   ! Number of particles at start of simulation
    INTEGER          :: npt_max     ! Maximum number of particles (fixed)
    REAL             :: npt_factor  ! Maximum number of particles (multiple of source bunch)
  END TYPE BeamBunchTemplate

  ! Variables
  TYPE(Interaction), PRIVATE, ALLOCATABLE :: interactions(:)
  DOUBLE PRECISION, PRIVATE :: residual_gas_pressure, residual_gas_temperature

  ! Interfaces
  INTERFACE interact_all
    MODULE PROCEDURE interact_wrapper
  END INTERFACE

  ! Public procedures
  PUBLIC :: construct_interactions
  PUBLIC :: destruct_interactions
  PUBLIC :: interact_all
  PUBLIC :: interact_setnp0

CONTAINS

  ! Initialize the interaction bunches and their interactions
  !  - bunches is an array of all bunches
  !  - Nbunch is the number of bunches, which gets updated with the new
  !      interaction bunches added in
  !  - uses an array of templates to create new bunches for interaction products
  !  - uses an array of interaction definitions that is global to this module
  SUBROUTINE construct_interactions(bunches, n_bunches)
    ! Parameters
    TYPE(BeamBunch), INTENT(INOUT) :: bunches(:)
    INTEGER,         INTENT(INOUT) :: n_bunches
    ! Variables
    TYPE(BeamBunchTemplate), ALLOCATABLE :: bunch_templates(:)

    ! Read in the interaction parameters and new bunch templates
    CALL read_interaction_data(interactions, bunch_templates)

    ! Initialize the individual interactions
    CALL interact_initialize_all()

    ! Create interaction bunches for each input bunch
    CALL interact_createbunches(bunches, n_bunches, bunch_templates)

    ! Initialize random number generator
    CALL random_seed()

  END SUBROUTINE construct_interactions



  ! Wrapper for initialization
  SUBROUTINE interact_initialize_all()
    ! Variables
    INTEGER :: i

    ! Initialize the individual interactions
    DO i = 1, SIZE(interactions)
      CALL interact_initialize(interactions(i))
    ENDDO

  END SUBROUTINE interact_initialize_all



  ! Generic intialization routine
  SUBROUTINE interact_initialize(current_interaction)
    ! Parameters
    TYPE(Interaction), INTENT(INOUT) :: current_interaction
    ! Variables
    DOUBLE PRECISION :: z_start

    ! Check interaction type
    IF(current_interaction%type == IT_DISSOCIATION) THEN
      CALL dissociate_initialize(current_interaction)
    ELSEIF(current_interaction%type == IT_CAPTURE) THEN
      CALL capture_initialize(current_interaction)
    ELSEIF(current_interaction%type == IT_CAPTURESPLIT) THEN
      CALL capturesplit_initialize(current_interaction)
    ELSE
      STOP "ERROR: Invalid interaction type."
    ENDIF

    ! Set location for first interaction
    z_start = 0.0d0
    CALL interact_setnext(current_interaction, z_start)

  END SUBROUTINE interact_initialize



  ! Wrapper for reading in all relevant interaction data from input files
  SUBROUTINE read_interaction_data(interaction_in, templates_in)
    ! Parameters
    TYPE(Interaction), ALLOCATABLE :: interaction_in(:)
    TYPE(BeamBunchTemplate), ALLOCATABLE :: templates_in(:)

    ! Read in the interaction parameters and new bunch templates
    CALL openinteractions_Input()
    CALL read_interactions(interaction_in)
    CALL read_bunchtemplates(templates_in)
    CALL closeinteractions_Input()

  END SUBROUTINE read_interaction_data



  ! Read in interaction parameters for all interactions
  SUBROUTINE read_interactions(interactions_in)
    ! Parameters
    TYPE(Interaction), ALLOCATABLE, INTENT(OUT) :: interactions_in(:)
    ! Variables
    INTEGER :: i, interaction_count

    ! Load number of interactions and common parameters from input file
    CALL readinteractions_Input(interaction_count, &
                                residual_gas_pressure, residual_gas_temperature)

    ! Create array of interaction types
    ALLOCATE(interactions_in(interaction_count))

    ! Define each interaction
    DO i = 1, interaction_count
      CALL read_interaction(interactions_in(i))
    ENDDO

  END SUBROUTINE read_interactions



  ! Read in interaction parameters for a single interaction
  SUBROUTINE read_interaction(current_interaction)
    ! Parameters
    TYPE(Interaction), INTENT(OUT) :: current_interaction

    CALL readinteraction_Input(current_interaction%type, &
                               current_interaction%source_bunch_id, &
                               current_interaction%target_bunch_id, &
                               current_interaction%interval, &
                               current_interaction%peak_energy, &
                               current_interaction%peak_cross)

  END SUBROUTINE read_interaction



  ! Read in bunch definitions for all interaction bunches
  SUBROUTINE read_bunchtemplates(templates)
    ! Parameters
    TYPE(BeamBunchTemplate), ALLOCATABLE, INTENT(OUT) :: templates(:)
    ! Variables
    INTEGER :: i, template_count

    ! Load number of interaction bunches
    CALL readinteractionbunchcount_Input(template_count)
    ALLOCATE(templates(template_count))

    ! Define each interaction
    DO i = 1, template_count
      CALL read_bunchtemplate(templates(i))
    ENDDO

  END SUBROUTINE read_bunchtemplates



  ! Read in bunch definition for a single interaction bunch
  SUBROUTINE read_bunchtemplate(current_template)
    ! Parameters
    TYPE(BeamBunchTemplate), INTENT(OUT) :: current_template

    ! Set up templates for creating new bunches
    CALL readinteractionbunch_Input(current_template%current, &
                                    current_template%energy, &
                                    current_template%mass, &
                                    current_template%charge, &
                                    current_template%phase, &
                                    current_template%npt_start, &
                                    current_template%npt_max, &
                                    current_template%npt_factor)

    ! Convert units (input file has eV, internal energy is keV)
    current_template%energy = current_template%energy/1000.0d0

  END SUBROUTINE read_bunchtemplate



  ! Destructor
  SUBROUTINE destruct_interactions()

    IF(ALLOCATED(interactions)) DEALLOCATE(interactions)

  END SUBROUTINE destruct_interactions



  ! Create interaction bunches for each input bunch
  SUBROUTINE interact_createbunches(bunches, n_bunches, templates)
    ! Parameters
    TYPE(BeamBunch),         INTENT(INOUT) :: bunches(:)
    INTEGER,                 INTENT(INOUT) :: n_bunches
    TYPE(BeamBunchTemplate), INTENT(INOUT) :: templates(:)
    ! Variables
    INTEGER :: i, new_bunch_count, new_bunch_id, array_size, max_size

    ! Find the largest array size from existing bunches
    max_size = 0
    DO i = 1, n_bunches
      array_size = SIZE(bunches(i)%Pts1, 2)
      IF(array_size > max_size) max_size = array_size
    ENDDO

    ! Loop through all bunch templates
    new_bunch_count = n_bunches
    DO i = 1, SIZE(templates)
      ! Add one more bunch
      new_bunch_id = new_bunch_count + 1
      IF(debug_mode) PRINT *, "Creating interaction bunch: ", new_bunch_id
      ! Check whether or not to scale the array size from the existing bunches
      IF(templates(i)%npt_max == 0) &
        templates(i)%npt_max = NINT(max_size * templates(i)%npt_factor)
      ! Create the new bunch
      CALL interact_createbunch(bunches(new_bunch_id), templates(i))
      new_bunch_count = new_bunch_count + 1
    ENDDO

    ! Adjust the number of bunches
    n_bunches = new_bunch_count

  END SUBROUTINE interact_createbunches



  ! Create a new (empty) beam bunch object based on a template
  SUBROUTINE interact_createbunch(new_bunch, bunch_template)
    ! Parameters
    TYPE(BeamBunch),         INTENT(INOUT) :: new_bunch
    TYPE(BeamBunchTemplate), INTENT(IN)    :: bunch_template

    ! Create a new (empty) bunch
    CALL construct_BeamBunch(new_bunch, &
                             bunch_template%current, bunch_template%energy, &
                             bunch_template%mass, bunch_template%charge, &
                             bunch_template%npt_start, bunch_template%phase)

    ! Set the number of particles at the start
    CALL setnpt_BeamBunch(new_bunch, bunch_template%npt_start)

    ! Initialize the particle array
    ALLOCATE(new_bunch%Pts1(6, bunch_template%npt_max))
    new_bunch%Pts1 = 0.0d0

    IF(debug_mode) THEN
      PRINT *, " - starting current (A): ", bunch_template%current
      PRINT *, " - starting energy (keV): ", bunch_template%energy
      PRINT *, " - particle rest mass (MeV): ", bunch_template%mass/1.0d6
      PRINT *, " - particle charge (e): ", bunch_template%charge
      PRINT *, " - starting phase (rad): ", bunch_template%phase
      PRINT *, " - initial particle count: ", bunch_template%npt_start
      PRINT *, " - maximum particle count: ", bunch_template%npt_max
    ENDIF

  END SUBROUTINE interact_createbunch



  ! Set the values of maximum particle count for the np0 variables
  SUBROUTINE interact_setnp0(np0, nplocal0, bunches)
    ! Parameters
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: np0(:)
    INTEGER, ALLOCATABLE, INTENT(INOUT) :: nplocal0(:)
    TYPE(BeamBunch),      INTENT(IN)    :: bunches(:)
    ! Variables
    INTEGER, ALLOCATABLE :: interaction_bunch_ids(:)
    INTEGER :: i, current_bunch, ierr
    ! MPI
    include 'mpif.h'

    ! Get list of interaction bunch ids
    CALL interact_getinteractionbunchids(interaction_bunch_ids)

    ! Set np0 sizes for all interaction bunches based on particle array size
    DO i = 1, SIZE(interaction_bunch_ids)
      current_bunch = interaction_bunch_ids(i)
      nplocal0(current_bunch) = SIZE(bunches(current_bunch)%Pts1, 2)
      CALL MPI_ALLREDUCE(nplocal0(current_bunch), np0(current_bunch), &
                         1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    ENDDO

  END SUBROUTINE interact_setnp0



  ! Wrapper function to call all relevant interactions
  !  - receives a subset of bunches to act on
  SUBROUTINE interact_wrapper(bunches, current_location, dzz)
    ! Parameters
    TYPE(BeamBunch),  INTENT(INOUT) :: bunches(:)
    DOUBLE PRECISION, INTENT(IN)    :: current_location, dzz
    ! Variables
    INTEGER :: i

    ! Loop through all interactions
    DO i = 1, SIZE(interactions)
      IF(current_location + dzz >= interactions(i)%next) THEN
        ! Implement the interaction
        IF(debug_mode) PRINT *, "Interaction at z = ", current_location
        CALL interact(bunches, interactions(i), current_location, dzz)
        ! Set location for next interaction
        CALL interact_setnext(interactions(i), current_location)
        IF(debug_mode) PRINT *, " - next interaction: ", interactions(i)%next
      ENDIF
    ENDDO

  END SUBROUTINE interact_wrapper



  ! Generic interaction call
  SUBROUTINE interact(bunches, current_interaction, current_location, dzz)
    ! Parameters
    TYPE(BeamBunch),   INTENT(INOUT) :: bunches(:)
    TYPE(Interaction), INTENT(INOUT) :: current_interaction
    DOUBLE PRECISION,  INTENT(IN)    :: current_location, dzz
    ! Variables
    TYPE(BeamBunch)      :: source_bunch
    INTEGER              :: number_to_move
    INTEGER, ALLOCATABLE :: particles_to_move(:)
    DOUBLE PRECISION     :: probability

    ! Report interaction type
    IF(debug_mode) THEN
      IF(current_interaction%type == IT_DISSOCIATION) THEN
        PRINT *, " - dissociation interaction"
      ELSEIF(current_interaction%type == IT_CAPTURE) THEN
        PRINT *, " - electron capture interaction"
      ELSEIF(current_interaction%type == IT_CAPTURESPLIT) THEN
        PRINT *, " - electron capture and split interaction"
      ENDIF
    ENDIF

    ! Get source bunch for current interaction
    source_bunch = bunches(current_interaction%source_bunch_id)

    ! Calculate probability of interaction
    probability = interact_getprobability(source_bunch, current_interaction, &
                                          current_location, dzz)

    ! Translate probability to number of particles to interact
    number_to_move = INT(probability*source_bunch%Nptlocal)
    ALLOCATE(particles_to_move(number_to_move))
    CALL interact_randomselect(source_bunch%Nptlocal, particles_to_move)

    ! Check interaction type and move particles accordingly
    IF(current_interaction%type == IT_DISSOCIATION) THEN
      CALL dissociate_move(bunches, current_interaction, particles_to_move)
    ELSEIF(current_interaction%type == IT_CAPTURE) THEN
      CALL capture_move(bunches, current_interaction, particles_to_move)
    ELSEIF(current_interaction%type == IT_CAPTURESPLIT) THEN
      CALL capturesplit_move(bunches, current_interaction, particles_to_move)
    ELSE
      STOP "ERROR: Invalid interaction type."
    ENDIF

  END SUBROUTINE interact



  !-----------------------------------------------------------------------------
  ! Dissociation of composite ions into smaller ions



  ! Initialize dissociation interaction
  SUBROUTINE dissociate_initialize(current_interaction)
    ! Parameters
    TYPE(Interaction), INTENT(INOUT) :: current_interaction

    IF(debug_mode) THEN
      PRINT *, "Initializing dissociation interaction:"
      PRINT *, " - source bunch ID: ", current_interaction%source_bunch_id
      PRINT *, " - target bunch ID: ", current_interaction%target_bunch_id
      PRINT *, " - peak interaction energy (keV): ", current_interaction%peak_energy
      PRINT *, " - peak interaction cross-section (m^2): ", current_interaction%peak_cross
      PRINT *, " - interaction interval (m): ", current_interaction%interval
    ENDIF

  END SUBROUTINE dissociate_initialize



  ! Simple disocciation interaction: one source to two identical products
  SUBROUTINE dissociate_move(bunches, current_interaction, particle_ids)
    ! Parameters
    TYPE(BeamBunch), TARGET, INTENT(INOUT) :: bunches(:)
    TYPE(Interaction),       INTENT(IN)    :: current_interaction
    INTEGER,                 INTENT(IN)    :: particle_ids(:)
    ! Variables
    TYPE(BeamBunch), POINTER :: source_bunch, target_bunch

    ! Get source and target bunch
    source_bunch => bunches(current_interaction%source_bunch_id)
    target_bunch => bunches(current_interaction%target_bunch_id)

    ! Add new particles to target bunch (copy parameters from source bunch)
    CALL addto_BeamBunch(target_bunch, source_bunch, particle_ids)
    CALL addto_BeamBunch(target_bunch, source_bunch, particle_ids)
    ! Remove particles from source bunch
    CALL removefrom_BeamBunch(source_bunch, particle_ids)

  END SUBROUTINE dissociate_move



  !-----------------------------------------------------------------------------
  ! Electron capture by ions in the beam



  ! Initialize electron capture interaction
  SUBROUTINE capture_initialize(current_interaction)
    ! Parameters
    TYPE(Interaction), INTENT(INOUT) :: current_interaction

    IF(debug_mode) THEN
      PRINT *, "Initializing electron capture interaction:"
      PRINT *, " - source bunch ID: ", current_interaction%source_bunch_id
      PRINT *, " - target bunch ID: ", current_interaction%target_bunch_id
      PRINT *, " - peak interaction energy (keV): ", current_interaction%peak_energy
      PRINT *, " - peak interaction cross-section (m^2): ", current_interaction%peak_cross
      PRINT *, " - interaction interval (m): ", current_interaction%interval
    ENDIF

  END SUBROUTINE capture_initialize



  ! Electron capture interaction: one source to one product
  SUBROUTINE capture_move(bunches, current_interaction, particle_ids)
    ! Parameters
    TYPE(BeamBunch), TARGET, INTENT(INOUT) :: bunches(:)
    TYPE(Interaction),       INTENT(IN)    :: current_interaction
    INTEGER,                 INTENT(IN)    :: particle_ids(:)
    ! Variables
    TYPE(BeamBunch), POINTER :: source_bunch, target_bunch

    ! Get source and target bunch
    source_bunch => bunches(current_interaction%source_bunch_id)
    target_bunch => bunches(current_interaction%target_bunch_id)

    ! Add new particles to target bunch (copy parameters from source bunch)
    CALL addto_BeamBunch(target_bunch, source_bunch, particle_ids)
    ! Remove particles from source bunch
    CALL removefrom_BeamBunch(source_bunch, particle_ids)

  END SUBROUTINE capture_move



  !-----------------------------------------------------------------------------
  ! Electron capture by compound ions in the beam with subsequent split



  ! Initialize electron capture and split interaction
  SUBROUTINE capturesplit_initialize(current_interaction)
    ! Parameters
    TYPE(Interaction), INTENT(INOUT) :: current_interaction

    IF(debug_mode) THEN
      PRINT *, "Initializing electron capture and split interaction:"
      PRINT *, " - source bunch ID: ", current_interaction%source_bunch_id
      PRINT *, " - target bunch ID: ", current_interaction%target_bunch_id
      PRINT *, " - peak interaction energy (keV): ", current_interaction%peak_energy
      PRINT *, " - peak interaction cross-section (m^2): ", current_interaction%peak_cross
      PRINT *, " - interaction interval (m): ", current_interaction%interval
    ENDIF

  END SUBROUTINE capturesplit_initialize



  ! Electron capture and split interaction: one source to two identical products
  SUBROUTINE capturesplit_move(bunches, current_interaction, particle_ids)
    ! Parameters
    TYPE(BeamBunch), TARGET, INTENT(INOUT) :: bunches(:)
    TYPE(Interaction),       INTENT(IN)    :: current_interaction
    INTEGER,                 INTENT(IN)    :: particle_ids(:)
    ! Variables
    TYPE(BeamBunch), POINTER :: source_bunch, target_bunch

    ! Get source and target bunch
    source_bunch => bunches(current_interaction%source_bunch_id)
    target_bunch => bunches(current_interaction%target_bunch_id)

    ! Add new particles to target bunch (copy parameters from source bunch)
    CALL addto_BeamBunch(target_bunch, source_bunch, particle_ids)
    CALL addto_BeamBunch(target_bunch, source_bunch, particle_ids)
    ! Remove particles from source bunch
    CALL removefrom_BeamBunch(source_bunch, particle_ids)

  END SUBROUTINE capturesplit_move



  !-----------------------------------------------------------------------------
  ! General functions and procedures



  ! Select a certain number of random particles from a bunch
  !  - receives the total number of particles
  !  - returns an array of particle numbers selected
  !  - the number of particles selected is defined by
  !      the size of the selected_particles array
  SUBROUTINE interact_randomselect(n_total, selected_particles)
    ! Parameters
    INTEGER, INTENT(IN)  :: n_total
    INTEGER, INTENT(OUT) :: selected_particles(:)
    ! Variables
    INTEGER :: n_select, n_unshuffled, current_selection, i, swap
    INTEGER, ALLOCATABLE :: id(:)
    REAL :: random

    ! Get number of particles to select
    n_select = SIZE(selected_particles)
    IF(n_select > n_total) STOP "ERROR: Too many particles to be selected."
    IF(debug_mode) THEN
      PRINT *, " - selecting particles to move..."
      PRINT *, "   - number of particles to select from: ", n_total
      PRINT *, "   - number of particles to select: ", n_select
    ENDIF

    ! Make list of all particle ids
    ALLOCATE(id(n_total))
    DO i = 1, n_total
      id(i) = i
    ENDDO

    ! Start shuffling the ids with a Knuth shuffle
    n_unshuffled = n_total
    DO i = 1, n_select
      CALL random_number(random)
      current_selection = INT(random*n_unshuffled + 1)
      swap = id(n_unshuffled)
      id(n_unshuffled) = id(current_selection)
      id(current_selection) = swap
      n_unshuffled = n_unshuffled - 1
    ENDDO

    ! The last n_select ids have already been shuffled, are unique and random
    selected_particles = id(n_unshuffled+1:n_total)

    ! Sort the results
    CALL qsort(selected_particles)

  END SUBROUTINE interact_randomselect



  ! Return the bunch IDs of all interaction bunches
  SUBROUTINE interact_getinteractionbunchids(bunch_ids)
    ! Parameters
    INTEGER, ALLOCATABLE, INTENT(OUT) :: bunch_ids(:)
    ! Variables
    INTEGER, ALLOCATABLE :: id_list(:)
    INTEGER :: i

    ! Maximum size of list
    ALLOCATE(id_list(SIZE(interactions)))

    ! Check all interactions
    DO i = 1, SIZE(interactions)
      id_list(i) = interactions(i)%target_bunch_id
    ENDDO

    ! Sort results and remove duplicates
    CALL qsort(id_list)
    CALL remove_duplicates(id_list, bunch_ids)

  END SUBROUTINE interact_getinteractionbunchids



  ! Get energy of bunch
  !  - assume energy spread can be neglected for the purpose of calculating
  !    the interaction probability
  !  - this simplification means we can work from the reference particle
  DOUBLE PRECISION FUNCTION interact_getenergy(bunch)
    ! Parameters
    TYPE(BeamBunch), INTENT(IN) :: bunch

    ! Return value
    interact_getenergy = getreferenceenergy_BeamBunch(bunch)

  END FUNCTION interact_getenergy



  ! Get speed of bunch
  !  - assume energy spread can be neglected for the purpose of calculating
  !    the interaction probability
  !  - this simplification means we can work from the reference particle
  DOUBLE PRECISION FUNCTION interact_getspeed(bunch)
    ! Parameters
    TYPE(BeamBunch), INTENT(IN) :: bunch

    ! Return value
    interact_getspeed = getreferencespeed_BeamBunch(bunch)

  END FUNCTION interact_getspeed



  ! Get phase of bunch
  !  - assume phase spread can be neglected
  !  - this simplification means we can work from the reference particle
  DOUBLE PRECISION FUNCTION interact_getphase(bunch)
    ! Parameters
    TYPE(BeamBunch), INTENT(IN) :: bunch

    ! Return value
    interact_getphase = getreferencephase_BeamBunch(bunch)

  END FUNCTION interact_getphase



  ! Get particle current of bunch
  DOUBLE PRECISION FUNCTION interact_getparticlecurrent(bunch)
    ! Parameters
    TYPE(BeamBunch), INTENT(IN) :: bunch

    ! Return value
    interact_getparticlecurrent = getparticlecurrent_BeamBunch(bunch)

  END FUNCTION interact_getparticlecurrent



  ! Get the residual gas pressure at a certain location
  DOUBLE PRECISION FUNCTION interact_getresidualpressure(current_location)
    ! Parameters
    DOUBLE PRECISION, INTENT(IN) :: current_location ! currently unused

    ! Fixed value for now
    interact_getresidualpressure = residual_gas_pressure

  END FUNCTION interact_getresidualpressure



  ! Get the residual gas temperature at a certain location
  DOUBLE PRECISION FUNCTION interact_gettemperature(current_location)
    ! Parameters
    DOUBLE PRECISION, INTENT(IN) :: current_location ! currently unused

    ! Fixed value for now
    interact_gettemperature = residual_gas_temperature

  END FUNCTION interact_gettemperature



  ! Calculate a cross-section from an energy
  DOUBLE PRECISION FUNCTION interact_getcrosssection(energy, &
                                                     peak_energy, peak_cross)

    ! Parameters
    DOUBLE PRECISION, INTENT(IN) :: energy
    DOUBLE PRECISION, INTENT(IN) :: peak_energy
    DOUBLE PRECISION, INTENT(IN) :: peak_cross

    ! Simple linear/reciprocal curve
    IF(energy <= peak_energy) THEN
      interact_getcrosssection = peak_cross * energy / peak_energy
    ELSE
      interact_getcrosssection = peak_cross * peak_energy / energy
    ENDIF

  END FUNCTION interact_getcrosssection



  ! Calculate the probability of an interaction
  DOUBLE PRECISION FUNCTION interact_getprobability(source_bunch, &
                                                    current_interaction, &
                                                    current_location, dzz)
    ! Parameters
    TYPE(BeamBunch),   INTENT(IN) :: source_bunch
    TYPE(Interaction), INTENT(IN) :: current_interaction
    DOUBLE PRECISION,  INTENT(IN) :: current_location, dzz
    ! Variables
    DOUBLE PRECISION :: gas_pressure, gas_temperature, gas_density
    DOUBLE PRECISION :: bunch_energy, bunch_speed, particle_current
    DOUBLE PRECISION :: distance_interval, time_interval, cross_section

    ! Get values of required variables
    gas_pressure      = interact_getresidualpressure(current_location) ! Pa
    gas_temperature   = interact_gettemperature(current_location)      ! K
    gas_density       = gas_pressure/(k_B*gas_temperature)             ! m^-3
    bunch_energy      = interact_getenergy(source_bunch)               ! eV
    bunch_speed       = interact_getspeed(source_bunch)                ! m/s
    particle_current  = interact_getparticlecurrent(source_bunch)      ! s^-1
    distance_interval = current_interaction%interval                   ! m
    time_interval     = distance_interval / bunch_speed                ! s

    ! Calculate cross-section (m^2)
    cross_section   = interact_getcrosssection(bunch_energy, &
                                               current_interaction%peak_energy, &
                                               current_interaction%peak_cross)

    ! Probability of collision in this time-step
    interact_getprobability = gas_density * particle_current &
                            * cross_section &
                            * distance_interval * time_interval

    ! Check for unphysical probability
    interact_getprobability = MIN(interact_getprobability, 1.0d0)

    IF(debug_mode) THEN
      PRINT *, " - bunch reference energy (keV):     ", bunch_energy/1.0d3
      PRINT *, " - bunch reference speed (m/s):      ", bunch_speed
      PRINT *, " - particle current (1/s):           ", particle_current
      PRINT *, " - residual gas pressure (Pa):       ", gas_pressure
      PRINT *, " - residual gas temperature (deg C): ", gas_temperature - 273.15
      PRINT *, " - residual gas density (1/m^3):     ", gas_density
      PRINT *, " - interaction cross-section (m^2):  ", cross_section
      PRINT *, " - interaction interval (m):         ", distance_interval
      PRINT *, " - interaction time (s):             ", time_interval
      PRINT *, " - interaction probability:          ", interact_getprobability
    ENDIF

  END FUNCTION interact_getprobability



  ! Set the next location of an interaction
  SUBROUTINE interact_setnext(current_interaction, current_location)
    ! Parameters
    TYPE(Interaction), INTENT(INOUT) :: current_interaction
    DOUBLE PRECISION :: current_location

    current_interaction%next = current_location + current_interaction%interval

  END SUBROUTINE interact_setnext



END MODULE InteractionClass
