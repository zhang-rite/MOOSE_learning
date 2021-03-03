[Mesh]
  [cartesian]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10000 4000 1000 0.2 1000 4000 10000'
    ix = '5 10 50 2 50 10 5'
    dy = '10000 3000 600 1400'
    iy = '5 10 60 10'
    dz = '10000 4000 1000 0.2 1000 4000 10000'
    iz = '5 10 50 2 50 10 5'
  []
  [shift_down]
    type = TransformGenerator
    transform = TRANSLATE
    vector_value = '-15000.1 -15000 -15000.1'
    input = cartesian
  []
  [well_body]
    type = SubdomainBoundingBoxGenerator
    input = shift_down
    block_id = 1
    bottom_left = '-0.1 -1900 -0.1'
    top_right = '0.1 0 0.1'
  []
  [delete_well_body]
    type = BlockDeletionGenerator
    input = well_body
    block_id = 1
  []
  [aquifer]
    type = SubdomainBoundingBoxGenerator
    block_id = 2
    bottom_left = '-16000 -1900 -16000'
    top_right = '16000 -1700 16000'
    input = delete_well_body
  []
  [injection_area]
    type = ParsedGenerateSideset
    combinatorial_geometry = 'x>-0.15 & x<0.15 & z>-0.15 & z<0.15'
    included_subdomain_ids = '2'
    new_sideset_name = injection_area
    input = aquifer
  []
  [rename]
    type = RenameBlockGenerator
    old_block_id = '0 2'
    new_block_name = 'caps aquifer'
    input = injection_area
  []
[]

[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  PorousFlowDictator = 'dictator'
  gravity = '0 -9.81 0'
  biot_coefficient = '1.0'
[]

[Variables]
  [pwater]
  []
  [disp_x]
    scaling = '1e-7'
  []
  [disp_y]
    scaling = '1e-7'
  []
  [disp_z]
    scaling = '1e-7'
  []
[]

[ICs]
  [pwater]
    type = FunctionIC
    function = ppic
    variable = pwater
  []
[]

[Functions]
  [perm_md_fcn]
    type = PiecewiseMultilinear
    data_file = k.csv
  []
  [ppic]
    type = ParsedFunction
    value = '101325-1000*9.81*(0-y)'
  []
  [weight_fcn]
    type = ParsedFunction
    vars = 'g rho0 biot'
    vals = '9.81 1000 1.0'
    value = 'biot*(101325-rho0*g*y)+2386.0*g*y'
  []
[]

[Kernels]
  # [./mass_water_dot]
  # type = PorousFlowMassTimeDerivative
  # fluid_component = 0
  # # use_displaced_mesh = false
  # variable = pwater
  # [../]
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
    gravity = '0 -9.81 0'
  []
  [grad_stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    component = 1
    use_displaced_mesh = false
  []
  [grad_stress_z]
    type = StressDivergenceTensors
    variable = disp_z
    component = 2
    use_displaced_mesh = false
  []
  [gravity_y]
    type = Gravity
    variable = disp_y
    value = -9.81
    use_displaced_mesh = false
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    biot_coefficient = 1.0
    use_displaced_mesh = false
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    biot_coefficient = 1.0
    use_displaced_mesh = false
    component = 1
  []
  [poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    biot_coefficient = 1.0
    use_displaced_mesh = false
    component = 2
  []
[]

[BCs]
  [top]
    type = FunctionDirichletBC
    variable = pwater
    boundary = 'top'
    use_displaced_mesh = false
    function = 101325-y*9810
  []
  [xmax_drained]
    type = FunctionDirichletBC
    variable = pwater
    boundary = 'right left back front'
    use_displaced_mesh = false
    function = 101325-y*9810
  []
  [fixed_outer_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right left'
    use_displaced_mesh = false
    value = 0
  []
  [fixed_outer_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'back front'
    use_displaced_mesh = false
    value = 0
  []
  [fixed_well]
    type = DirichletBC
    variable = disp_x
    boundary = 'injection_area'
    use_displaced_mesh = false
    value = 0
  []
  [fixed_well2]
    type = DirichletBC
    variable = disp_z
    boundary = 'injection_area'
    use_displaced_mesh = false
    value = 0
  []
  [fixed_bottom_z]
    type = DirichletBC
    variable = disp_y
    boundary = 'bottom'
    use_displaced_mesh = false
    value = 0
  []
  [top_z]
    type = Pressure
    variable = disp_y
    use_displaced_mesh = false
    boundary = 'top'
    component = 1
    function = -2386.0*9.81*y
  []
[]

[AuxVariables]
  [stress_xx]
    order = FIRST
    family = MONOMIAL
  []
  [stress_zz]
    order = FIRST
    family = MONOMIAL
  []
  [stress_yy]
    order = FIRST
    family = MONOMIAL
  []
  [strain_yy]
    family = MONOMIAL
    order = FIRST
  []
  [strain_yy_0]
    family = MONOMIAL
    order = FIRST
  []
  [strain_yy_2]
    family = MONOMIAL
    order = FIRST
  []
  [strain_xx]
    family = MONOMIAL
    order = FIRST
  []
  [strain_zz]
    family = MONOMIAL
    order = FIRST
  []
  [perm_md]
    family = MONOMIAL
    order = FIRST
  []
  [perm]
    family = MONOMIAL
    order = FIRST
  []
[]

[AuxKernels]
  [perm_md]
    type = FunctionAux
    function = perm_md_fcn
    variable = perm_md
    execute_on = 'initial'
  []
  [perm]
    type = ParsedAux
    variable = perm
    args = 'perm_md'
    function = 'perm_md' # '9.869233e-13*perm_md'
    execute_on = 'initial'
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
  [stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  []
  [strain_yy]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy
    index_i = 1
    index_j = 1
  []
  [strain_yy_0]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_yy_0
    index_i = 1
    index_j = 1
    execute_on = 'initial'
  []
  [strain_yy_2]
    type = ParsedAux
    variable = strain_yy_2
    args = 'strain_yy_0 strain_yy'
    function = 'strain_yy-strain_yy_0'
  []
  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 2
    index_j = 2
  []
  [strain_xx]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_xx
    index_i = 0
    index_j = 0
  []
[]

[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
[]

[Modules]
  [FluidProperties]
    [water]
      type = SimpleFluidProperties
      bulk_modulus = 2.1e9
      density0 = 1000
      viscosity = 8.9e-4
      cv = 4149.0
      cp = 4149.0
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
    [co2]
      type = SimpleFluidProperties
      bulk_modulus = 2.1e8
      density0 = 600
      viscosity = 1e-5
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = '45'
  []
  [massfrac]
    type = PorousFlowMassFraction
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [porosity_reservoir]
    type = PorousFlowPorosityConst
    porosity = '0.2'
  []
  [permeability_aquifer]
    type = PorousFlowPermeabilityConst
    block = 'aquifer'
    permeability = '2E-12 0 0   0 2E-12 0   0 0 2E-12'
  []
  [permeability_caps]
    type = PorousFlowPermeabilityConst
    block = 'caps'
    permeability = '1E-19 0 0   0 1E-19 0   0 0 1E-19'
  []
  [relperm]
    type = PorousFlowRelativePermeabilityCorey
    n = 1
    phase = 0
  []
  [ppss]
    type = PorousFlow1PhaseP
    porepressure = 'pwater'
    capillary_pressure = pc
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 6.0E9
    poissons_ratio = 0.2
  []
  [strain]
    # # use_displaced_mesh = false
    type = ComputeSmallStrain # ComputeAxisymmetricRZSmallStrain
    eigenstrain_names = 'ini_stress'
  []
  [ini_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'weight_fcn 0 0  0 weight_fcn 0  0 0 weight_fcn' 
    eigenstrain_name = ini_stress
  []
  [stress]
    type = ComputeLinearElasticStress
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    outputs = 'exodus'
  []
  [density]
    type = GenericConstantMaterial
    prop_names = 'density'
    prop_values = '2386.0' # = (1-0.1)*2540 + 0.1*999.526
  []
[]

[Preconditioning]
  [andy]
    type = SMP
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -pc_factor_shift_type  -snes_atol -snes_rtol'
    petsc_options_value = 'gmres      lu       mumps        NONZERO                5.7e1     1E-8'
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
  []
[]

[Executioner]
  type = Steady
  solve_type = Newton
  l_max_its = 200
  nl_max_its = 500
[]

[Outputs]
  exodus = true
  perf_graph = true
  [console]
    type = Console
    output_linear = true
    output_nonlinear = true
  []
[]
