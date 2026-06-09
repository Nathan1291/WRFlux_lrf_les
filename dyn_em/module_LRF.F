
! AUTHOR: Yi-Ling Hwong (CCRC), with modifications by Tim Raupach, 2020.
!
! Functions for forcing vertical temperature or moisture tendencies to implement
! Herman & Kuang (2013) Linear Response Function (LRF).

module module_LRF
  
CONTAINS

  subroutine zero_tendencies(        &
       thforceten, qvforceten,       &
       threlaxten, qvrelaxten,       &
       urelaxten, vrelaxten,         &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       ips, ipe, jps, jpe, kps, kpe, &
       its, ite, jts, jte, kts, kte)
    ! Zero tendencies so other tendency functions can be additive.
    !
    ! thforceten: Potential temperature (theta) tendency field [K s-1].
    ! qvforceten: Water vapour tendency field [kg kg-1 s-1].
    ! threlaxten: Potential temperature (theta) relaxation tendency field [K s-1].
    ! qvrelaxten: Water vapour relaxation tendency field [kg kg-1 s-1].
    ! urelaxten: U wind relaxation tendency field [m s-2].
    ! vrelaxten: V wind relaxation tendency field [m s-2].
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.

    implicit none
    
    integer, intent(in) :: ids,ide,jds,jde,kds,kde, &
                           ims,ime,jms,jme,kms,kme, &  
                           ips,ipe,jps,jpe,kps,kpe, &
                           its,ite,jts,jte,kts,kte  

    ! Input/output variables.
    real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: &
         thforceten, qvforceten, threlaxten, qvrelaxten, urelaxten, vrelaxten

    ! Local variables.
    integer :: i,j,k

    ! Zero all values at mass points.
    do i = its, min(ite, ide-1)
       do k = kts, min(kte, kde-1)
          do j = jts, min(jte, jde-1)    
             thforceten(i,k,j) = 0
             qvforceten(i,k,j) = 0
             threlaxten(i,k,j) = 0
             qvrelaxten(i,k,j) = 0
             urelaxten(i,k,j) = 0
             vrelaxten(i,k,j) = 0
          end do
       end do
    end do
    
  end subroutine zero_tendencies
    
  subroutine force_LRF( &
       perturb_t, perturb_q, & 
       p_pert, TtendAmp, QtendAmp, &
       p, pi, thten, qvten, &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       ips, ipe, jps, jpe, kps, kpe, &
       its, ite, jts, jte, kts, kte)
    ! force_LRF: force vertical temperature and moisture tendencies, to
    ! implement Herman & Kuang (2013) Linear Response Function (LRF).
    !
    ! ---
    ! perturb_t: Perturb temperature tendencies?
    ! perturb_q: Perturb moisture tendencies?
    ! p_pert: The pressure around which to perturb [hPa].
    ! TtendAmp: Temperature (dT/dt) perturbation [K day-1].
    ! QtendAmp: Water vapour mixing ratio (dq/dt) perturbation [kg kg-1 day-1].
    ! p: Pressure field [Pa].
    ! pi: The Exner function (pi_phy in WRF).
    ! thten: Potential temperature (theta) tendency field [K s-1].
    ! qvten: Water vapour tendency field [kg kg-1 s-1].
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.
    
    implicit none

    ! Input variables.
    logical, intent(in)                                  :: perturb_t, perturb_q
    real, intent(in)                                     :: p_pert
    real, intent(in)                                     :: TtendAmp, QtendAmp
    real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: p, pi
    integer, intent(in)                                  :: ids,ide,jds,jde,kds,kde, &
                                                            ims,ime,jms,jme,kms,kme, &  
                                                            ips,ipe,jps,jpe,kps,kpe, &
                                                            its,ite,jts,jte,kts,kte  

    ! Input/output variables.
    real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: thten ! dT/dt (K s-1).
    real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: qvten ! dq/dt (kg kg-1 s-1).
    
    ! Local variables.
    integer :: i,j,k
    real    :: forcing

    ! If neither temperature or moisture are to be perturbed, return.
    if ((perturb_t .eqv. .false.) .and. (perturb_q .eqv. .false.)) return

    ! The forcing form is from the one used in Appendix A of Kuang 2012. Their
    ! forcing is around a vertical level and is the mean of a delta function and
    ! gamma function. We adapt the form to use only use the gamma function part
    ! and perturb around a pressure, not a level. This is because the pressure
    ! levels in the model may not exactly align with the pressure around which
    ! we want to perturb.

    ! Loop over unstaggered grid.
    do i = its, min(ite, ide-1)
       do k = kts, min(kte, kde-1)
          do j = jts, min(jte, jde-1)

             forcing = exp(-((p(i,k,j)-(p_pert*100))/7500.)**2.)

             ! Temperature forcing as a forcing on potential temperature (theta)
             ! tendency. Convert dT/dt to theta tendency by dividing by the
             ! exner function.
             if(perturb_t .eqv. .true.) then
                thten(i,k,j) = thten(i,k,j) + &
                     (((TtendAmp/86400) * forcing) / pi(i,k,j))
             endif
             
             ! Moisture forcing as a forcing on moisture tendency.
             if(perturb_q .eqv. .true.) then
                qvten(i,k,j) = qvten(i,k,j) + ((QtendAmp/86400) * forcing)
             endif
          enddo
       enddo
    enddo
  end subroutine force_LRF

subroutine relax_var_to_profile(    &
       ten, var, p, tau,              &
       target_profile, extend_relaxation, &
     relax_means, is_temp, pi,      &
       ids, ide, jds, jde, kds, kde,   &
       ims, ime, jms, jme, kms, kme,   &
       ips, ipe, jps, jpe, kps, kpe,   &
       its, ite, jts, jte, kts, kte)
    ! Relax given variable toward a profile read from file. If extend_relaxation
    ! is .false. then relax winds above 160 hPa to maximum relaxation at/above
    ! 100 hPa.  If 'extend_relaxation' is .true. then relax the whole
    ! profile. Targets should be given for the whole profile regardless of the
    ! value of extend_relaxation.
    !
    ! ---
    ! ten: The variable tendency field. 
    ! var: The variable itself. 
    ! p: Hydrostatic pressure [Pa].
    ! tau: Relaxation time, over which a full relaxation is made [s].
    ! target_profile: The target profile [m s-1].
    ! extend_relaxation: Extend relaxation to cover whole profile?
    ! relax_means: Relax the mean of each layer, instead of point-by-point?
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.
    ! 
    ! Optional Variables for temperature conversion
    ! is_temp: Check whether temp is being relaxed. 
    ! pi: The Exner function (pi_phy in WRF).

    implicit none

    real, dimension(ims:ime, kms:kme, jms:jme), intent(in)    :: var, p
    real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ten
    real, intent(in)                 :: tau
    real, dimension(kte), intent(in) :: target_profile
    logical, intent(in)              :: extend_relaxation
    logical, intent(in)              :: relax_means
    integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe, &
         its, ite, jts, jte, kts, kte
    logical, optional, intent(in)    :: is_temp
    real, optional, dimension(ims:ime,kms:kme,jms:jme), intent(in)    :: pi
    integer :: i,j,k
    real :: inv_tau
    real, dimension(kms:kme) :: means
    real :: relax_value

    if(relax_means .eqv. .true.) then
       call horizontal_means(v=var, out=means, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            ips=ips, ipe=ipe, jps=jps, jpe=jpe, kps=kps, kpe=kpe, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    endif

    ! Relax variable. 
    do j=jts,min(jte, jde-1)
       do k=kts,min(kte, kde-1)
          do i=its,min(ite, ide-1)
             ! Relax variable only above 160 hPa, or everywhere if
             ! extend_relaxation is true.
             if ((extend_relaxation .eqv. .true.) .or. (p(i,k,j) .lt. 16000)) then
                if ((extend_relaxation .eqv. .true.) .or. (p(i,k,j) .le. 10000)) then
                   inv_tau = 1/tau
                else
                   ! Between 160 hPa and 100 hPa, vary 1/tau linearly from 0 s-1
                   ! at 160 hPa to 1/tau s-1 at 100 hPa.
                   inv_tau = (1/tau)*(16000 - p(i,k,j))/6000
                endif

                if(relax_means .eqv. .true.) then
                   ! Apply relaxation to shift horizontal means to reference point.
                   relax_value = (inv_tau * -(means(k) - target_profile(k)))
                else
                   ! Apply relaxation point-to-point.
                   relax_value = (inv_tau * -(var(i,k,j) - target_profile(k)))
                endif

                if (present(is_temp)) then
                  if (is_temp) then
                     relax_value = relax_value / pi(i,k,j)
                  endif
               endif

                ten(i,k,j) = ten(i,k,j) + relax_value
             endif
          enddo
       enddo
    enddo
  end subroutine relax_var_to_profile

  subroutine relax_wind_to_profile(    &
       ten, wind, p, tau,              &
       target_wind, extend_relaxation, &
       relax_mean_wind,                &
       ids, ide, jds, jde, kds, kde,   &
       ims, ime, jms, jme, kms, kme,   &
       ips, ipe, jps, jpe, kps, kpe,   &
       its, ite, jts, jte, kts, kte)
    ! Relax winds toward a profiles read from file. If extend_relaxation
    ! is .false. then relax winds above 160 hPa to maximum relaxation at/above
    ! 100 hPa.  If 'extend_relaxation' is .true. then relax the whole
    ! profile. Targets should be given for the whole profile regardless of the
    ! value of extend_relaxation.
    !
    ! ---
    ! ten: The wind tendency field [m s-2].
    ! wind: The wind field (at mass points) [m s-1].
    ! p: Hydrostatic pressure [Pa].
    ! tau: Relaxation time, over which a full relaxation is made [s].
    ! target_wind: The target profile [m s-1].
    ! extend_relaxation: Extend relaxation to cover whole profile?
    ! relax_mean_wind: Relax the mean of each layer, instead of point-by-point?
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.

    implicit none
 
    real, dimension(ims:ime, kms:kme, jms:jme), intent(in)    :: wind, p
    real, dimension(ims:ime, kms:kme, jms:jme), intent(inout) :: ten
    real, intent(in)                 :: tau
    real, dimension(kte), intent(in) :: target_wind
    logical, intent(in)              :: extend_relaxation
    logical, intent(in)              :: relax_mean_wind
    integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe, &
         its, ite, jts, jte, kts, kte

    integer :: i,j,k
    real :: inv_tau
    real, dimension(kms:kme) :: means
    
    if(relax_mean_wind .eqv. .true.) then
       call horizontal_means(v=wind, out=means, &
            ids=ids, ide=ide, jds=jds, jde=jde, kds=kds, kde=kde, &
            ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
            ips=ips, ipe=ipe, jps=jps, jpe=jpe, kps=kps, kpe=kpe, &
            its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte)
    endif
    
    ! Relax winds. 
    do j=jts,min(jte, jde-1)
       do k=kts,min(kte, kde-1)
          do i=its,min(ite, ide-1)
             ! Relax winds only above 160 hPa, or everywhere if
             ! extend_relaxation is true.
             if ((extend_relaxation .eqv. .true.) .or. (p(i,k,j) .lt. 16000)) then
                if ((extend_relaxation .eqv. .true.) .or. (p(i,k,j) .le. 10000)) then
                   ! Above 100 hPa, or if extend_relaxation is true, use the
                   ! full inverse relaxation time, ie the specified tau value.
                   inv_tau = 1/tau
                else
                   ! Between 160 hPa and 100 hPa, vary 1/tau linearly from 0 s-1
                   ! at 160 hPa to 1/tau s-1 at 100 hPa.
                   inv_tau = (1/tau)*(16000 - p(i,k,j))/6000
                endif

                if(relax_mean_wind .eqv. .true.) then
                   ! Apply relaxation to shift horizontal means to reference point.
                   ten(i,k,j) = ten(i,k,j) + (inv_tau * -(means(k) - target_wind(k)))
                else
                   ! Apply relaxation point-to-point.
                   ten(i,k,j) = ten(i,k,j) + (inv_tau * -(wind(i,k,j) - target_wind(k)))
                endif
             endif
          enddo
       enddo
    enddo
  end subroutine relax_wind_to_profile
  
  subroutine horizontal_means(v, out, &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       ips, ipe, jps, jpe, kps, kpe, &
       its, ite, jts, jte, kts, kte)
    ! Calculate horizontal mean values per level.
    !
    ! v: The variable for which to find horizontal means.
    ! out: Output array, one value per level in v.
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.

    USE module_dm, ONLY: wrf_dm_sum_real
    
    implicit none
 
    real, dimension(ims:ime, kms:kme, jms:jme), intent(in) :: v
    real, dimension(kms:kme), intent(inout) :: out
    integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe, &
         its, ite, jts, jte, kts, kte
    
    integer :: i, j, k
    integer :: stag_x, stag_y
    integer :: i_from, i_to, j_from, j_to, k_from, k_to
    real :: sum, count

    do k = kts, min(kte, kde-1)
       sum = 0
       count = 0
       
       do j = jts, min(jte, jde-1)
          do i = its, min(ite, ide-1)
             sum = sum + v(i,k,j)
             count = count + 1
          enddo
       enddo
       
       ! Compute count and sum over whole layer in parallel.
       sum = wrf_dm_sum_real(sum)
       count = wrf_dm_sum_real(count)
       
       ! Check that all domain points have been counted.
       if(.not. count .eq. (ide-ids) * (jde-jds)) then
          call wrf_error_fatal("Wind nudging grid mean domain size mismatch.")
       endif

       out(k) = sum / count
    enddo
  end subroutine horizontal_means
  
  subroutine read_profile_file(num_levels, profile_file, vals)
    ! Read a profile from file. Expect one value per line.
    ! ---
    ! num_levels: The number of levels (values) to read.
    ! profile_file: The file to read from.
    ! vals: The read values (output).

    implicit none
    integer, intent(in) :: num_levels
    character(len=*), intent(in) :: profile_file
    real, dimension(num_levels), intent(inout) :: vals

    integer :: k
    
    open(unit=10, file=profile_file, form='formatted', status='old')
    rewind(10)

    k = 1
    do while(k .le. num_levels)
       read(10,*) vals(k)
       k = k + 1
    enddo

    close(unit=10, status='keep')
  end subroutine read_profile_file
  
  subroutine read_wind_target_profiles(num_levels, u_profile_file, &
       v_profile_file, u_val, v_val)
    ! Read U and V profiles from two files.
    ! ---
    ! num_levels: The number of levels (values) to read, one per line.
    ! u_profile_file: The file to read from for U.
    ! v_profile_file: The file to read from for V.
    ! u_val: The read values (output) for U.
    ! v_val: The read values (output) for V.
    
    implicit none
    integer, intent(in) :: num_levels
    character(len=*), intent(in) :: u_profile_file, v_profile_file
    real, dimension(num_levels), intent(inout) :: u_val, v_val
    write(6,*) 'Loading target profiles for wind relaxation.'

    call read_profile_file(num_levels, u_profile_file, u_val)
    call read_profile_file(num_levels, v_profile_file, v_val)
  end subroutine read_wind_target_profiles

  subroutine read_tq_target_profiles(num_levels, t_profile_file, &
       q_profile_file, t_val, q_val)
    ! Read U and V profiles from two files.
    ! ---
    ! num_levels: The number of levels (values) to read, one per line.
    ! t_profile_file: The file to read from for T.
    ! q_profile_file: The file to read from for q.
    ! t_val: The read values (output) for T.
    ! q_val: The read values (output) for q.
    
    implicit none
    integer, intent(in) :: num_levels
    character(len=*), intent(in) :: t_profile_file, q_profile_file
    real, dimension(num_levels), intent(inout) :: t_val, q_val
    write(6,*) 'Loading target profiles for T and q stratospheric relaxation.'

    call read_profile_file(num_levels, t_profile_file, t_val)
    call read_profile_file(num_levels, q_profile_file, q_val)
  end subroutine read_tq_target_profiles

  subroutine fix_radiation(          &
       p, pi, rthraten,              &
       ids, ide, jds, jde, kds, kde, &
       ims, ime, jms, jme, kms, kme, &
       ips, ipe, jps, jpe, kps, kpe, &
       its, ite, jts, jte, kts, kte)

    ! Set the radiation profile to a constant cooling profile as per Herman and
    ! Kuang (2013).
    !
    ! p: hydrostatic pressure at mass points (Pa).
    ! pi: the Exner function (absolute temperature/potential temperature).
    ! rthraten: theta tendency due to radiation; set by this function.
    ! ids, ide, jds, jde: kds, kde: Horizontal domain dimensions.
    ! ims, ime, jms, jme, kms, kme: Memory dimensions.
    ! ips, ipe, jps, jpe, kps, kpe: Patch dimensions.
    ! its, ite, jts, jte, kts, kte: Tile dimensions.
    
    implicit none
    
    real, dimension(ims:ime,kms:kme,jms:jme), intent(in) :: p, pi
    real, dimension(ims:ime,kms:kme,jms:jme), intent(inout) :: rthraten
    integer, intent(in) :: &
         ids, ide, jds, jde, kds, kde, &
         ims, ime, jms, jme, kms, kme, &
         ips, ipe, jps, jpe, kps, kpe, &
         its, ite, jts, jte, kts, kte
    
    integer :: i,j,k
    real :: t_tend
    
    do j = jts, min(jte, jde-1)
       do k = kts, min(kte, kde-1)
          do i = its, min(ite, ide-1)

             if (p(i,k,j) .ge. 20000) then             
                ! From surface to 200 hPa, constant radiative cooling value
                ! of Qrad = -1.5 K day-1 (converted to K s-1).
                t_tend = -1.5/86400
             else if (p(i,k,j) .le. 10000) then
                ! Above 100 hPa, radiative cooling is zero.
                t_tend = 0.0
             else
                ! Linear change between 200 hPa and 100 hPa.
                t_tend = -((1.5*p(i,k,j)/10000 - 1.5) / 86400)   
             endif
             
             ! Convert from temperature tendency to theta tendency.
             rthraten(i,k,j) = t_tend/pi(i,k,j) 
          enddo
       enddo
    enddo
  end subroutine fix_radiation
  
end module module_LRF
