module soilTemperatureTreatmentStreamMod

  ! Read the time-varying untreated (T0.00) soil temperature used by the
  ! deep-soil heater controller. The stream contains one gridcell field,
  ! T_SOIL_REFERENCE (K), evaluated at the configured control depth.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_strdata_mod
  use shr_mct_mod
  use mct_mod
  use spmdMod, only : mpicom, masterproc, comp_id
  use elm_varctl, only : iulog, inst_name, soil_heating_reference_file, &
       soil_heating_stream_year_first, soil_heating_stream_year_last, &
       soil_heating_model_year_align
  use decompMod, only : bounds_type, ldecomp, gsmap_lnd_gdc2glo
  use domainMod, only : ldomain

  implicit none
  private
  save

  public :: soil_temperature_treatment_stream_init
  public :: soil_temperature_treatment_stream_interp
  public :: soil_heating_reference_grc

  type(shr_strdata_type) :: sdat
  real(r8), allocatable :: soil_heating_reference_grc(:)

contains

  subroutine soil_temperature_treatment_stream_init(bounds)
    use elm_time_manager, only : get_calendar
    use ncdio_pio, only : pio_subsystem
    use shr_pio_mod, only : shr_pio_getiotype
    use ndepStreamMod, only : elm_domain_mct

    type(bounds_type), intent(in) :: bounds
    type(mct_ggrid) :: dom_elm

    allocate(soil_heating_reference_grc(bounds%begg:bounds%endg))
    soil_heating_reference_grc(:) = 0._r8

    call elm_domain_mct(bounds, dom_elm)
    call shr_strdata_create(sdat, name='elm_soil_heating_reference', &
         pio_subsystem=pio_subsystem, &
         pio_iotype=shr_pio_getiotype(inst_name), &
         mpicom=mpicom, compid=comp_id, &
         gsmap=gsmap_lnd_gdc2glo, ggrid=dom_elm, &
         nxg=ldomain%ni, nyg=ldomain%nj, &
         yearFirst=soil_heating_stream_year_first, &
         yearLast=soil_heating_stream_year_last, &
         yearAlign=soil_heating_model_year_align, &
         offset=0, &
         domFilePath='', domFileName=trim(soil_heating_reference_file), &
         domTvarName='time', domXvarName='lon', domYvarName='lat', &
         domAreaName='area', domMaskName='mask', &
         filePath='', filename=(/trim(soil_heating_reference_file)/), &
         fldListFile='T_SOIL_REFERENCE', &
         fldListModel='T_SOIL_REFERENCE', &
         fillalgo='none', mapalgo='bilinear', &
         calendar=get_calendar(), taxmode='extend')

    if (masterproc) then
       write(iulog,*) 'Deep-soil heating reference stream: ', &
            trim(soil_heating_reference_file)
       call shr_strdata_print(sdat, 'ELM deep-soil heating reference data')
    end if
  end subroutine soil_temperature_treatment_stream_init

  subroutine soil_temperature_treatment_stream_interp(bounds)
    use elm_time_manager, only : get_curr_date

    type(bounds_type), intent(in) :: bounds
    integer :: year, mon, day, sec, mcdate
    integer :: g, ig

    call get_curr_date(year, mon, day, sec)
    mcdate = year * 10000 + mon * 100 + day
    call shr_strdata_advance(sdat, mcdate, sec, mpicom, &
         'soil_temperature_treatment')

    ig = 0
    do g = bounds%begg, bounds%endg
       ig = ig + 1
       soil_heating_reference_grc(g) = sdat%avs(1)%rAttr(1, ig)
    end do
  end subroutine soil_temperature_treatment_stream_interp

end module soilTemperatureTreatmentStreamMod
