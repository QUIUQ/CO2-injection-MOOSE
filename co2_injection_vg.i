[Mesh]
    type = GeneratedMesh
    dim = 2
    ny = 25
    nx = 50
    ymax = 100
    xmin = 0.1
    xmax = 1200
    bias_x = 1.05
    coord_type = RZ
    rz_coord_axis = Y
[]

[GlobalParams]
    PorousFlowDictator = dictator
    gravity = '0 -9.81 0'
    temperature_unit = Celsius
[]

[Variables]
    [pgas]
        #scaling = 1e-7
    []
    [zi]
        initial_condition = 0
    []
[]

[ICs]
    [pp_initial]
        type = FunctionIC
        variable = pgas
        function = pp_ic
    []
[]

[Functions]
    [pp_ic]
        type = SolutionFunction
        solution = soln
    []
    [injection_rate]
        type = ParsedFunction
        symbol_values = injection_area
        symbol_names = area
        expression = '-5/area'
    []
[]


[UserObjects]
    [soln]
        type = SolutionUserObject
        mesh = initial_condition.e
        system_variables = porepressure
    []
    [dictator]
        type = PorousFlowDictator
        porous_flow_vars = 'pgas zi'
        number_fluid_phases = 2
        number_fluid_components = 2
    []
    [pc]
        type = PorousFlowCapillaryPressureVG
        alpha = 1e-5
        m = 0.5
        sat_lr = 0.3
    []
    [pc0]
        type = PorousFlowCapillaryPressureConst
        pc = 0
    []
    [fs]
        type = PorousFlowBrineCO2
        brine_fp = brine
        co2_fp = co2
        capillary_pressure = pc
    []
[]

[AuxVariables]
    [temperature]
        initial_condition = 50
    []
    [xnacl]
        initial_condition = 0.1
    []
    [pressure_gas]
        order = CONSTANT
        family = MONOMIAL
    []
    [pressure_water]
        order = CONSTANT
        family = MONOMIAL
    []
    [saturation_gas]
        order = CONSTANT
        family = MONOMIAL
    []
    [saturation_water]
        order = CONSTANT
        family = MONOMIAL
    []
    [density_water]
        order = CONSTANT
        family = MONOMIAL
    []
    [density_gas]
        order = CONSTANT
        family = MONOMIAL
    []
    [x0_water]
        order = CONSTANT
        family = MONOMIAL
    []
    [x0_gas]
        order = CONSTANT
        family = MONOMIAL
    []
    [x1_water]
        order = CONSTANT
        family = MONOMIAL
    []
    [x1_gas]
        order = CONSTANT
        family = MONOMIAL
    []
    [cap_pre]
        order = CONSTANT
        family = MONOMIAL
    []
[]

[Kernels]
    [mass0]
        type = PorousFlowMassTimeDerivative
        variable = pgas
        fluid_component = 0
    []
    [flux0]
        type = PorousFlowAdvectiveFlux
        variable = pgas
        fluid_component = 0
    []
    [mass1]
        type = PorousFlowMassTimeDerivative
        variable = zi
        fluid_component = 1
    []
    [flux1]
        type = PorousFlowAdvectiveFlux
        variable = zi
        fluid_component = 1
    []
[]

[AuxKernels]
    [pressure_water]
        type = PorousFlowPropertyAux
        variable = pressure_water
        property = pressure
        phase = 0
        execute_on = timestep_end
    []
    [pressure_gas]
        type = PorousFlowPropertyAux
        variable = pressure_gas
        property = pressure
        phase = 1
        execute_on = timestep_end
    []
    [saturation_water]
        type = PorousFlowPropertyAux
        variable = saturation_water
        property = saturation
        phase = 0
        execute_on = timestep_end
    []
    [saturation_gas]
        type = PorousFlowPropertyAux
        variable = saturation_gas
        property = saturation
        phase = 1
        execute_on = timestep_end
    []
    [density_water]
        type = PorousFlowPropertyAux
        variable = density_water
        property = density
        phase = 0
        execute_on = timestep_end
    []
    [density_gas]
        type = PorousFlowPropertyAux
        variable = density_gas
        property = density
        phase = 1
        execute_on = timestep_end
    []
    [x1_water]
        type = PorousFlowPropertyAux
        variable = x1_water
        property = mass_fraction
        phase = 0
        fluid_component = 1
        execute_on = timestep_end
    []
    [x1_gas]
        type = PorousFlowPropertyAux
        variable = x1_gas
        property = mass_fraction
        phase = 1
        fluid_component = 1
        execute_on = timestep_end
    []
    [x0_water]
        type = PorousFlowPropertyAux
        variable = x0_water
        property = mass_fraction
        phase = 0
        fluid_component = 0
        execute_on = timestep_end
    []
    [x0_gas]
        type = PorousFlowPropertyAux
        variable = x0_gas
        property = mass_fraction
        phase = 1
        fluid_component = 0
        execute_on = timestep_end
    []
    [capillary_pressure]
        type = PorousFlowPropertyAux
        variable = cap_pre
        property = capillary_pressure
        execute_on = timestep_end
    []
[]

[BCs]
    [methane_injection]
        type = PorousFlowSink
        boundary = left
        variable = zi
        flux_function = injection_rate
        fluid_phase = 1
    []
    [outer_pressure_fixed]
        type = FunctionDirichletBC
        variable = pgas
        boundary = right
        function = pp_ic
    []
[]

[FluidProperties]
    [water]
        type = Water97FluidProperties
    []
    [watertab]
        type = TabulatedBicubicFluidProperties
        fp = water
        out_of_bounds_behavior = ignore

    []
    [brine]
        type = BrineFluidProperties
        water_fp = watertab
    []
    [co2sw]
        type = CO2FluidProperties
    []
    [co2]
        type = TabulatedBicubicFluidProperties
        fp = co2sw
        out_of_bounds_behavior = set_to_closest_bound
    []

[]

[Materials]
    [temperature]
        type = PorousFlowTemperature
        temperature = temperature
    []
    [brineco2]
        type = PorousFlowFluidState
        gas_porepressure = pgas
        z = zi
        temperature = temperature
        xnacl = xnacl
        capillary_pressure = pc
        fluid_state = fs
    []
    [porosity]
        type = PorousFlowPorosityConst
        porosity = 0.1
    []
    [permeability]
        type = PorousFlowPermeabilityConst
        permeability = '1.5833e-13 0 0  0 1.5833e-13 0  0 0 1.5833e-13' #160.43 mD
    []
    [relperm_water]
        type = PorousFlowRelativePermeabilityVG
        m = 0.45946
        phase = 0
        s_res = 0.300
        sum_s_res = 0.3
      []
      [relperm_gas]
        type = PorousFlowRelativePermeabilityBC
        phase = 1
        s_res = 0.0
        sum_s_res = 0.3
        lambda = 2
        nw_phase = true
      []
[]

[Preconditioning]
  active = preferred_but_might_not_be_installed
  [basic]
    type = SMP
    full = true
  []
  [preferred_but_might_not_be_installed]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
[]

  [Executioner]
    type = Transient
    solve_type = NEWTON
    end_time = 31536000
    nl_abs_tol = 1e-7
    nl_rel_tol = 1e-5
    dtmax = 1e5
    nl_max_its = 15
    automatic_scaling = true
    #compute_scaling_once = true
    [TimeStepper]
      type = IterationAdaptiveDT
      dt = 1e1
      growth_factor = 1.25
    []
  []

  [Postprocessors]
    [mass_ph0]
      type = PorousFlowFluidMass
      fluid_component = 0
      execute_on = 'initial timestep_end'
    []
    [mass_ph1]
      type = PorousFlowFluidMass
      fluid_component = 1
      execute_on = 'initial timestep_end'
    []
    [injection_area]
      type = AreaPostprocessor
      boundary = left
      execute_on = initial
    []
  []

  [Outputs]
    execute_on = 'initial timestep_end'
    exodus = true
    perf_graph = true
    [csv]
        type = CSV
        sync_times = '2.592e6 8.64e6 31536000'
        sync_only = true
    []
  []
