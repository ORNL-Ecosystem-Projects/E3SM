module MicrobeMethaneReactionMod

  ! Side-effect-free reaction kernel for the revised microbial methane backend.
  ! All carbon-bearing solid/organic inputs and tendencies use g C m-3 soil
  ! and all dissolved-gas inputs and tendencies use mol gas m-3. Rates use
  ! the corresponding units per second. The kernel never mutates its inputs.

  use shr_kind_mod, only : r8 => shr_kind_r8
  use elm_varcon, only : catomw, secspday
  use MicrobeMethaneParamsMod, only : MicrobeMethaneParamsType

  implicit none
  private
  save

  real(r8), parameter :: mmol_to_mol = 1.e-3_r8
  real(r8), parameter :: q10_temperature_interval = 10._r8

  type, public :: microbe_methane_reaction_state_type
     real(r8) :: dom_c = 0._r8
     real(r8) :: acetate_c = 0._r8
     real(r8) :: acetate_methanogen_c = 0._r8
     real(r8) :: h2_methanogen_c = 0._r8
     real(r8) :: aerobic_methanotroph_c = 0._r8
     real(r8) :: anaerobic_methanotroph_c = 0._r8
     real(r8) :: conc_ch4 = 0._r8
     real(r8) :: conc_o2 = 0._r8
     real(r8) :: conc_co2 = 0._r8
     real(r8) :: conc_h2 = 0._r8
  end type microbe_methane_reaction_state_type

  type, public :: microbe_methane_reaction_environment_type
     real(r8) :: soil_temperature = 273.15_r8
     real(r8) :: soil_ph = 7._r8
     ! Hydrology supplies these bounded scalars when Step 3 connects the
     ! partition/repartition logic. They keep that policy outside chemistry.
     real(r8) :: dom_fermentation_scalar = 0._r8
     real(r8) :: aerobic_acetate_oxidation_scalar = 0._r8
     ! Potential ELM decomposition, root-respiration, and nitrification O2
     ! demand. The reaction limiter competes this demand with its microbial
     ! aerobic reactions, following the shared-stress structure in CH4Mod.
     real(r8) :: elm_aerobic_o2_demand = 0._r8
  end type microbe_methane_reaction_environment_type

  type, public :: microbe_methane_reaction_rates_type
     ! Carbon-basis process extents in mol C m-3 s-1.
     real(r8) :: dom_to_acetate_c = 0._r8
     real(r8) :: acetogenesis_c = 0._r8
     real(r8) :: acetoclastic_methanogenesis_c = 0._r8
     real(r8) :: hydrogenotrophic_methanogenesis_c = 0._r8
     real(r8) :: aerobic_acetate_oxidation_c = 0._r8
     real(r8) :: aerobic_methane_oxidation_c = 0._r8
     ! Diagnostic copy after CH4 competition and before shared-O2 scaling.
     real(r8) :: aerobic_methane_oxidation_pre_o2_c = 0._r8
     real(r8) :: anaerobic_methane_oxidation_c = 0._r8
     real(r8) :: acetate_methanogen_mortality_c = 0._r8
     real(r8) :: h2_methanogen_mortality_c = 0._r8
     real(r8) :: aerobic_methanotroph_mortality_c = 0._r8
     real(r8) :: anaerobic_methanotroph_mortality_c = 0._r8
     ! Non-microbial O2 consumption accepted by the shared limiter
     ! [mol O2 m-3 s-1] and the dimensionless common O2 stress.
     real(r8) :: elm_aerobic_o2_consumption = 0._r8
     real(r8) :: oxygen_stress = 1._r8
     real(r8) :: effective_soil_ph = 0._r8
     real(r8) :: ph_response = 0._r8
  end type microbe_methane_reaction_rates_type

  type, public :: microbe_methane_reaction_tendencies_type
     real(r8) :: dom_c = 0._r8
     real(r8) :: acetate_c = 0._r8
     real(r8) :: acetate_methanogen_c = 0._r8
     real(r8) :: h2_methanogen_c = 0._r8
     real(r8) :: aerobic_methanotroph_c = 0._r8
     real(r8) :: anaerobic_methanotroph_c = 0._r8
     real(r8) :: conc_ch4 = 0._r8
     real(r8) :: conc_o2 = 0._r8
     real(r8) :: conc_co2 = 0._r8
     real(r8) :: conc_h2 = 0._r8
  end type microbe_methane_reaction_tendencies_type

  public :: computeMicrobeMethaneReactionTendencies
  public :: computeMicrobeMethanePotentialRates
  public :: limitMicrobeMethaneReactionRates
  public :: assembleMicrobeMethaneReactionTendencies
  public :: microbeMethanePHResponse
  public :: microbeMethaneQ10Response
  public :: microbeMethaneCarbonResidual

contains

  pure subroutine computeMicrobeMethaneReactionTendencies(state, environment, parameters, dt, rates, tendencies)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    type(microbe_methane_reaction_environment_type), intent(in) :: environment
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    real(r8), intent(in) :: dt
    type(microbe_methane_reaction_rates_type), intent(out) :: rates
    type(microbe_methane_reaction_tendencies_type), intent(out) :: tendencies

    call computeMicrobeMethanePotentialRates(state, environment, parameters, rates)
    call limitMicrobeMethaneReactionRates(state, environment, parameters, dt, rates)
    call assembleMicrobeMethaneReactionTendencies(parameters, rates, tendencies)
  end subroutine computeMicrobeMethaneReactionTendencies

  pure subroutine computeMicrobeMethanePotentialRates(state, environment, parameters, rates)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    type(microbe_methane_reaction_environment_type), intent(in) :: environment
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    type(microbe_methane_reaction_rates_type), intent(out) :: rates
    real(r8) :: dom_mmol_c, acetate_mmol_c
    real(r8) :: ch4_mmol, o2_mmol, co2_mmol, h2_mmol
    real(r8) :: acetate_methanogen_mol_c, h2_methanogen_mol_c
    real(r8) :: aerobic_methanotroph_mol_c, anaerobic_methanotroph_mol_c
    real(r8) :: acetate_feedback, co2_inhibition, o2_inhibition
    real(r8) :: dom_scalar, aerobic_acetate_scalar

    rates = microbe_methane_reaction_rates_type()
    dom_mmol_c = gCToMmolC(state%dom_c)
    acetate_mmol_c = gCToMmolC(state%acetate_c)
    ch4_mmol = 1000._r8 * max(0._r8, state%conc_ch4)
    o2_mmol = 1000._r8 * max(0._r8, state%conc_o2)
    co2_mmol = 1000._r8 * max(0._r8, state%conc_co2)
    h2_mmol = 1000._r8 * max(0._r8, state%conc_h2)
    acetate_methanogen_mol_c = gCToMolC(state%acetate_methanogen_c)
    h2_methanogen_mol_c = gCToMolC(state%h2_methanogen_c)
    aerobic_methanotroph_mol_c = gCToMolC(state%aerobic_methanotroph_c)
    anaerobic_methanotroph_mol_c = gCToMolC(state%anaerobic_methanotroph_c)

    rates%effective_soil_ph = effectiveSoilPH(environment%soil_ph, acetate_mmol_c, parameters)
    rates%ph_response = microbeMethanePHResponse(rates%effective_soil_ph, parameters%ph_min, &
         parameters%ph_opt, parameters%ph_max)
    dom_scalar = clampUnitInterval(environment%dom_fermentation_scalar)
    aerobic_acetate_scalar = clampUnitInterval(environment%aerobic_acetate_oxidation_scalar)
    acetate_feedback = parameters%acetate_feedback_half_saturation / &
         (acetate_mmol_c + parameters%acetate_feedback_half_saturation)
    co2_inhibition = max(0._r8, 1._r8 - min(1._r8, &
         co2_mmol / parameters%h2_methanogenesis_co2_inhibition_scale))
    o2_inhibition = max(0._r8, 1._r8 - min(1._r8, &
         o2_mmol / parameters%aom_o2_inhibition_scale))

    ! The legacy 2/3 factor is the acetate-C share of the DOM conversion.
    rates%dom_to_acetate_c = mmolCPerDayToMolCPerSecond((2._r8 / 3._r8) * &
         parameters%acetate_prod_max * monod(dom_mmol_c, parameters%k_acetate) * &
         microbeMethaneQ10Response(parameters%dom_to_acetate_q10, environment%soil_temperature, &
         parameters%reaction_t_ref) * rates%ph_response * acetate_feedback * &
         dom_scalar)

    rates%acetogenesis_c = mmolCPerDayToMolCPerSecond(parameters%acetogenesis_max * &
         monod(h2_mmol, parameters%k_acetogenesis_h2) * &
         monod(co2_mmol, parameters%k_acetogenesis_co2) * &
         microbeMethaneQ10Response(parameters%acetogenesis_q10, environment%soil_temperature, &
         parameters%reaction_t_ref) * rates%ph_response)

    ! This extent is CH4-C production. Applying the biomass yield to this
    ! extent makes the named growth rate the maximum specific growth rate;
    ! the legacy extra factor of four in biomass growth is not retained.
    rates%hydrogenotrophic_methanogenesis_c = &
         parameters%h2_methanogen_growth_rate / parameters%h2_methanogen_yield / secspday * &
         h2_methanogen_mol_c * monod(h2_mmol, parameters%k_h2_methanogenesis_h2) * &
         monod(co2_mmol, parameters%k_h2_methanogenesis_co2) * &
         microbeMethaneQ10Response(parameters%h2_methanogenesis_q10, environment%soil_temperature, &
         parameters%reaction_t_ref) * rates%ph_response * co2_inhibition

    rates%acetoclastic_methanogenesis_c = &
         parameters%acetate_methanogen_growth_rate / parameters%acetate_methanogen_yield / secspday * &
         acetate_methanogen_mol_c * monod(acetate_mmol_c, &
         parameters%k_acetoclastic_methanogenesis_acetate) * &
         microbeMethaneQ10Response(parameters%acetoclastic_methanogenesis_q10, &
         environment%soil_temperature, parameters%reaction_t_ref) * rates%ph_response

    rates%aerobic_acetate_oxidation_c = parameters%aerobic_acetate_oxidation_rate / secspday * &
         gCToMolC(state%acetate_c) * monod(o2_mmol, parameters%k_acetate_prod_o2) * &
         aerobic_acetate_scalar

    rates%aerobic_methane_oxidation_c = &
         parameters%aerobic_methanotroph_growth_rate / parameters%aerobic_methanotroph_yield / secspday * &
         aerobic_methanotroph_mol_c * monod(ch4_mmol, parameters%k_aerobic_oxidation_ch4) * &
         monod(o2_mmol, parameters%k_aerobic_oxidation_o2) * &
         microbeMethaneQ10Response(parameters%aerobic_oxidation_q10, environment%soil_temperature, &
         parameters%reaction_t_ref) * rates%ph_response

    rates%anaerobic_methane_oxidation_c = &
         parameters%anaerobic_methanotroph_growth_rate / parameters%anaerobic_methanotroph_yield / secspday * &
         anaerobic_methanotroph_mol_c * monod(ch4_mmol, parameters%k_anaerobic_oxidation_ch4) * &
         microbeMethaneQ10Response(parameters%anaerobic_oxidation_q10, environment%soil_temperature, &
         parameters%aom_t_ref) * rates%ph_response * o2_inhibition

    rates%acetate_methanogen_mortality_c = parameters%acetate_methanogen_death_rate / secspday * &
         acetate_methanogen_mol_c * rates%ph_response
    rates%h2_methanogen_mortality_c = parameters%h2_methanogen_death_rate / secspday * &
         h2_methanogen_mol_c * rates%ph_response
    rates%aerobic_methanotroph_mortality_c = parameters%aerobic_methanotroph_death_rate / secspday * &
         aerobic_methanotroph_mol_c * rates%ph_response
    rates%anaerobic_methanotroph_mortality_c = &
         parameters%anaerobic_methanotroph_death_rate / secspday * &
         anaerobic_methanotroph_mol_c * rates%ph_response
  end subroutine computeMicrobeMethanePotentialRates

  pure subroutine limitMicrobeMethaneReactionRates(state, environment, parameters, dt, rates)
    type(microbe_methane_reaction_state_type), intent(in) :: state
    type(microbe_methane_reaction_environment_type), intent(in) :: environment
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    real(r8), intent(in) :: dt
    type(microbe_methane_reaction_rates_type), intent(inout) :: rates
    real(r8) :: available_h2_rate, available_co2_rate, available_acetate_rate
    real(r8) :: available_ch4_rate, available_o2_rate
    real(r8) :: h2_demand, co2_demand, acetate_demand, ch4_demand, o2_demand
    real(r8) :: scale

    if (dt <= 0._r8) then
       rates = microbe_methane_reaction_rates_type()
       return
    end if

    ! Stage 1: DOM conversion. Its 1.5 mol-C donor requirement yields one
    ! mol-C acetate, 0.5 mol CO2, and 1/6 mol H2 per mol-C process extent.
    rates%dom_to_acetate_c = min(rates%dom_to_acetate_c, &
         gCToMolC(state%dom_c) / (1.5_r8 * dt))

    ! Stage 2: acetogenesis and hydrogenotrophic methanogenesis compete for
    ! the H2 and CO2 present after same-step DOM conversion.
    available_h2_rate = max(0._r8, state%conc_h2) / dt + rates%dom_to_acetate_c / 6._r8
    available_co2_rate = max(0._r8, state%conc_co2) / dt + 0.5_r8 * rates%dom_to_acetate_c
    h2_demand = 2._r8 * rates%acetogenesis_c + &
         4._r8 * rates%hydrogenotrophic_methanogenesis_c
    co2_demand = rates%acetogenesis_c + &
         (1._r8 + parameters%h2_methanogen_yield) * rates%hydrogenotrophic_methanogenesis_c
    scale = min(supplyScale(available_h2_rate, h2_demand), &
         supplyScale(available_co2_rate, co2_demand))
    rates%acetogenesis_c = rates%acetogenesis_c * scale
    rates%hydrogenotrophic_methanogenesis_c = rates%hydrogenotrophic_methanogenesis_c * scale

    ! Stage 3: both acetate consumers see acetate produced in stages 1-2.
    available_acetate_rate = gCToMolC(state%acetate_c) / dt + &
         rates%dom_to_acetate_c + rates%acetogenesis_c
    acetate_demand = rates%acetoclastic_methanogenesis_c + rates%aerobic_acetate_oxidation_c
    scale = supplyScale(available_acetate_rate, acetate_demand)
    rates%acetoclastic_methanogenesis_c = rates%acetoclastic_methanogenesis_c * scale
    rates%aerobic_acetate_oxidation_c = rates%aerobic_acetate_oxidation_c * scale

    ! Stage 4: aerobic and anaerobic oxidation compete for initial plus
    ! same-step methanogenic CH4 production.
    available_ch4_rate = max(0._r8, state%conc_ch4) / dt + &
         parameters%acetoclastic_methanogenesis_ch4_yield * &
         (1._r8 - parameters%acetate_methanogen_yield) * &
         rates%acetoclastic_methanogenesis_c + rates%hydrogenotrophic_methanogenesis_c
    ch4_demand = rates%aerobic_methane_oxidation_c + rates%anaerobic_methane_oxidation_c
    scale = supplyScale(available_ch4_rate, ch4_demand)
    rates%aerobic_methane_oxidation_c = rates%aerobic_methane_oxidation_c * scale
    rates%anaerobic_methane_oxidation_c = rates%anaerobic_methane_oxidation_c * scale
    rates%aerobic_methane_oxidation_pre_o2_c = rates%aerobic_methane_oxidation_c

    ! Stage 5: ELM decomposition, root respiration, nitrification, aerobic
    ! acetate oxidation, and aerobic methane oxidation compete for O2. ELM's
    ! carbon and nitrogen fluxes were resolved earlier using the preceding
    ! timestep's shared stress, as in the legacy CH4Mod coupling.
    available_o2_rate = max(0._r8, state%conc_o2) / dt
    o2_demand = max(0._r8, environment%elm_aerobic_o2_demand) + &
         parameters%aerobic_decomp_o2_c_ratio * rates%aerobic_acetate_oxidation_c + &
         parameters%aerobic_oxidation_o2_ch4_ratio * rates%aerobic_methane_oxidation_c
    scale = supplyScale(available_o2_rate, o2_demand)
    rates%oxygen_stress = scale
    rates%elm_aerobic_o2_consumption = &
         max(0._r8, environment%elm_aerobic_o2_demand) * scale
    rates%aerobic_acetate_oxidation_c = rates%aerobic_acetate_oxidation_c * scale
    rates%aerobic_methane_oxidation_c = rates%aerobic_methane_oxidation_c * scale

    ! Mortality cannot consume biomass produced in this call or the functional
    ! biomass seed. Reserving the floor here avoids a post-update clamp that
    ! would create carbon. Mortality carbon is returned to DOM when tendencies
    ! are assembled.
    rates%acetate_methanogen_mortality_c = min(rates%acetate_methanogen_mortality_c, &
         gCToMolC(max(0._r8, state%acetate_methanogen_c - parameters%mfg_biomass_min)) / dt)
    rates%h2_methanogen_mortality_c = min(rates%h2_methanogen_mortality_c, &
         gCToMolC(max(0._r8, state%h2_methanogen_c - parameters%mfg_biomass_min)) / dt)
    rates%aerobic_methanotroph_mortality_c = min(rates%aerobic_methanotroph_mortality_c, &
         gCToMolC(max(0._r8, state%aerobic_methanotroph_c - parameters%mfg_biomass_min)) / dt)
    rates%anaerobic_methanotroph_mortality_c = min(rates%anaerobic_methanotroph_mortality_c, &
         gCToMolC(max(0._r8, state%anaerobic_methanotroph_c - parameters%mfg_biomass_min)) / dt)
  end subroutine limitMicrobeMethaneReactionRates

  pure subroutine assembleMicrobeMethaneReactionTendencies(parameters, rates, tendencies)
    type(MicrobeMethaneParamsType), intent(in) :: parameters
    type(microbe_methane_reaction_rates_type), intent(in) :: rates
    type(microbe_methane_reaction_tendencies_type), intent(out) :: tendencies
    real(r8) :: acetate_methanogen_growth, h2_methanogen_growth
    real(r8) :: aerobic_methanotroph_growth, anaerobic_methanotroph_growth
    real(r8) :: acetate_nonbiomass, acetate_ch4, acetate_co2
    real(r8) :: aerobic_oxidation_co2, anaerobic_oxidation_co2

    tendencies = microbe_methane_reaction_tendencies_type()
    acetate_methanogen_growth = parameters%acetate_methanogen_yield * &
         rates%acetoclastic_methanogenesis_c
    ! All guild yields are biomass-C per carbon-substrate process extent.
    h2_methanogen_growth = parameters%h2_methanogen_yield * &
         rates%hydrogenotrophic_methanogenesis_c
    aerobic_methanotroph_growth = parameters%aerobic_methanotroph_yield * &
         rates%aerobic_methane_oxidation_c
    anaerobic_methanotroph_growth = parameters%anaerobic_methanotroph_yield * &
         rates%anaerobic_methane_oxidation_c
    acetate_nonbiomass = (1._r8 - parameters%acetate_methanogen_yield) * &
         rates%acetoclastic_methanogenesis_c
    acetate_ch4 = parameters%acetoclastic_methanogenesis_ch4_yield * acetate_nonbiomass
    acetate_co2 = (1._r8 - parameters%acetoclastic_methanogenesis_ch4_yield) * acetate_nonbiomass
    aerobic_oxidation_co2 = (1._r8 - parameters%aerobic_methanotroph_yield) * &
         rates%aerobic_methane_oxidation_c
    anaerobic_oxidation_co2 = (1._r8 - parameters%anaerobic_methanotroph_yield) * &
         rates%anaerobic_methane_oxidation_c

    tendencies%dom_c = catomw * (-1.5_r8 * rates%dom_to_acetate_c + &
         rates%acetate_methanogen_mortality_c + rates%h2_methanogen_mortality_c + &
         rates%aerobic_methanotroph_mortality_c + rates%anaerobic_methanotroph_mortality_c)
    tendencies%acetate_c = catomw * (rates%dom_to_acetate_c + rates%acetogenesis_c - &
         rates%acetoclastic_methanogenesis_c - rates%aerobic_acetate_oxidation_c)
    tendencies%acetate_methanogen_c = catomw * (acetate_methanogen_growth - &
         rates%acetate_methanogen_mortality_c)
    tendencies%h2_methanogen_c = catomw * (h2_methanogen_growth - &
         rates%h2_methanogen_mortality_c)
    tendencies%aerobic_methanotroph_c = catomw * (aerobic_methanotroph_growth - &
         rates%aerobic_methanotroph_mortality_c)
    tendencies%anaerobic_methanotroph_c = catomw * (anaerobic_methanotroph_growth - &
         rates%anaerobic_methanotroph_mortality_c)

    tendencies%conc_ch4 = acetate_ch4 + rates%hydrogenotrophic_methanogenesis_c - &
         rates%aerobic_methane_oxidation_c - rates%anaerobic_methane_oxidation_c
    tendencies%conc_o2 = -rates%elm_aerobic_o2_consumption - &
         parameters%aerobic_decomp_o2_c_ratio * &
         rates%aerobic_acetate_oxidation_c - parameters%aerobic_oxidation_o2_ch4_ratio * &
         rates%aerobic_methane_oxidation_c
    tendencies%conc_co2 = 0.5_r8 * rates%dom_to_acetate_c - rates%acetogenesis_c - &
         (1._r8 + parameters%h2_methanogen_yield) * &
         rates%hydrogenotrophic_methanogenesis_c + acetate_co2 + &
         rates%aerobic_acetate_oxidation_c + aerobic_oxidation_co2 + anaerobic_oxidation_co2
    tendencies%conc_h2 = rates%dom_to_acetate_c / 6._r8 - 2._r8 * rates%acetogenesis_c - &
         4._r8 * rates%hydrogenotrophic_methanogenesis_c
  end subroutine assembleMicrobeMethaneReactionTendencies

  pure real(r8) function microbeMethanePHResponse(ph, ph_min, ph_opt, ph_max) result(response)
    real(r8), intent(in) :: ph, ph_min, ph_opt, ph_max
    real(r8) :: numerator, denominator

    if (ph <= ph_min .or. ph >= ph_max) then
       response = 0._r8
       return
    end if
    numerator = (ph - ph_min) * (ph - ph_max)
    denominator = numerator - (ph - ph_opt) * (ph - ph_opt)
    if (abs(denominator) <= tiny(denominator)) then
       response = 0._r8
    else
       response = clampUnitInterval(numerator / denominator)
    end if
  end function microbeMethanePHResponse

  pure real(r8) function microbeMethaneQ10Response(q10, temperature, reference_temperature) result(response)
    real(r8), intent(in) :: q10, temperature, reference_temperature
    response = q10 ** ((temperature - reference_temperature) / q10_temperature_interval)
  end function microbeMethaneQ10Response

  pure real(r8) function microbeMethaneCarbonResidual(tendencies) result(residual)
    type(microbe_methane_reaction_tendencies_type), intent(in) :: tendencies
    residual = (tendencies%dom_c + tendencies%acetate_c + &
         tendencies%acetate_methanogen_c + tendencies%h2_methanogen_c + &
         tendencies%aerobic_methanotroph_c + tendencies%anaerobic_methanotroph_c) / catomw + &
         tendencies%conc_ch4 + tendencies%conc_co2
  end function microbeMethaneCarbonResidual

  pure real(r8) function effectiveSoilPH(soil_ph, acetate_mmol_c, parameters) result(ph)
    real(r8), intent(in) :: soil_ph, acetate_mmol_c
    type(MicrobeMethaneParamsType), intent(in) :: parameters

    if (soil_ph > parameters%acetate_ph_trigger .and. acetate_mmol_c > 0._r8) then
       ph = -log10(10._r8 ** (-soil_ph) + parameters%acidification_coefficient * acetate_mmol_c)
    else
       ph = soil_ph
    end if
  end function effectiveSoilPH

  pure real(r8) function monod(substrate, half_saturation) result(response)
    real(r8), intent(in) :: substrate, half_saturation
    real(r8) :: nonnegative_substrate

    nonnegative_substrate = max(0._r8, substrate)
    response = nonnegative_substrate / (nonnegative_substrate + half_saturation)
  end function monod

  pure real(r8) function clampUnitInterval(value) result(clamped)
    real(r8), intent(in) :: value
    clamped = min(1._r8, max(0._r8, value))
  end function clampUnitInterval

  pure real(r8) function supplyScale(available_rate, demand_rate) result(scale)
    real(r8), intent(in) :: available_rate, demand_rate
    if (demand_rate <= 0._r8) then
       scale = 1._r8
    else
       scale = clampUnitInterval(max(0._r8, available_rate) / demand_rate)
    end if
  end function supplyScale

  pure real(r8) function gCToMolC(carbon_mass) result(carbon_moles)
    real(r8), intent(in) :: carbon_mass
    carbon_moles = max(0._r8, carbon_mass) / catomw
  end function gCToMolC

  pure real(r8) function gCToMmolC(carbon_mass) result(carbon_mmoles)
    real(r8), intent(in) :: carbon_mass
    carbon_mmoles = 1000._r8 * gCToMolC(carbon_mass)
  end function gCToMmolC

  pure real(r8) function mmolCPerDayToMolCPerSecond(rate) result(converted_rate)
    real(r8), intent(in) :: rate
    converted_rate = max(0._r8, rate) * mmol_to_mol / secspday
  end function mmolCPerDayToMolCPerSecond

end module MicrobeMethaneReactionMod
