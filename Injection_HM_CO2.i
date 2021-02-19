# Two phase, temperature-dependent, with mechanics, radial with fine mesh, 
# constant injection of cold co2 into a overburden-reservoir-underburden containing mostly water
# species=0 is water
# species=1 is co2
# phase=0 is liquid, and since massfrac_ph0_sp0 = 1, this is all water
# phase=1 is gas, and since massfrac_ph1_sp0 = 0, this is all co2
#
# The mesh used below has very high resolution, so the simulation takes a long time to complete.
# Some suggested meshes of different resolution:
# nx=50, bias_x=1.2
# nx=100, bias_x=1.1
# nx=200, bias_x=1.05
# nx=400, bias_x=1.02
# nx=1000, bias_x=1.01
# nx=2000, bias_x=1.003
# [Mesh]
#   type = GeneratedMesh
#   dim = 2
#   nx = 100
#   bias_x = 1.05
#   xmin = 0
#   xmax = 5000
#   ny = 100
#   ymin = 0
#   ymax = 5000
# []

# [Mesh]
#  type = MeshGeneratorMesh
#  [./cartesian]
#    type = CartesianMeshGenerator
#    dim = 2
#    dx = '0.2 1000 4500 6000'
#    ix = '2 20 20 6'
#    dy = '3000 1500 500'
#    iy = '30 30 10'
#    # subdomain_id = '1 2'
#  [../]
# [] 

# [Mesh]
#  type = MeshGeneratorMesh
#  [./cartesian]
#    type = CartesianMeshGenerator
#    dim = 2
#    dx = '0.2 1000 4500 10000'
#    ix = '2 50 45 10'
#    dy = '3000 1500 500'
#    iy = '30 100 10'
#    # subdomain_id = '1 2'
#  [../]
# [] 



[Mesh]
  type = FileMesh
  file = Initialization_HM_CO2_exodus.e
[]



# [MeshModifiers]
#   [shift_down]
#     type = Transform
#     transform = TRANSLATE
# 	  vector_value = '0.1 -5000 0'
#   []
#   [aquifer]
#     type = SubdomainBoundingBox
#     block_id = 1
#     bottom_left = '0.1 -1500 0' 
#     top_right = '30000.1 -1000 0'
# 	depends_on = shift_down
#   []
#   [injection_area]
#     type = ParsedAddSideset
#     combinatorial_geometry = x<0.101
#     included_subdomain_ids = '1'
#     new_sideset_name = injection_area
#     depends_on = 'aquifer'
#   []
#   [rename]
#     type = RenameBlock
#     old_block_id = '0 1'
#     new_block_name = 'caps aquifer'
#     depends_on = 'injection_area'
#   []
# []


[Problem]
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_r disp_z'
  PorousFlowDictator = dictator
  gravity = '0 -9.81 0'
  biot_coefficient = 1.0
[]

[UserObjects]
  [./ini_soln]
    type = SolutionUserObject
    mesh = Initialization_HM_CO2_exodus.e
    timestep = 'LATEST'
    execute_on = 'INITIAL'
    # system_variables = u
  [../]
[]

[Variables]
  [./pwater]
    # initial_condition = 18.3e6
    initial_from_file_var = pwater
  [../]
  [./sgas]
    initial_condition = 0.0
    # initial_from_file_var = sgas
    
  [../]
  # [./temp]
  #   initial_condition = 358
  # [../]
  [./disp_r]
  [../]
  [./disp_z]
  [../]  
[]

[ICs]
  [pwater]
    type = FunctionIC
    function = ppic
    variable = pwater
  []
[]
[AuxVariables]
  [./stress_ini1]
    order = FIRST
    family = LAGRANGE
  [../]
  [./stress_ini2]
    order = FIRST
    family = LAGRANGE
  [../]
  [./stress_ini3]
    order = FIRST
    family = LAGRANGE
  [../]
  [./stress_ini0]
    initial_condition = 0
  [../]  

  [./rate]
  [../]
  # [./disp_z]
  # [../]
  [./massfrac_ph0_sp0]
    initial_condition = 1 # all H20 in phase=0
  [../]
  [./massfrac_ph1_sp0]
    initial_condition = 0 # no H2O in phase=1
  [../]
  [./pgas]
    family = MONOMIAL
    order = FIRST
  [../]
  [./pgas0]
    family = MONOMIAL
    order = FIRST
  [../]  
  [./pgas2]
    family = MONOMIAL
    order = FIRST
  [../]  
  [./swater]
    family = MONOMIAL
    order = FIRST
  [../]
  [./stress_rr]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_tt]
    order = FIRST
    family = MONOMIAL
  [../]
  [./stress_zz]
    order = FIRST
    family = MONOMIAL
  [../]
  [strain_zz]
    family = MONOMIAL
    order = FIRST
  []
  [strain_zz_0]
    family = MONOMIAL
    order = FIRST
  []    
  [strain_zz_2]
    family = MONOMIAL
    order = FIRST
    # execute_on = initial
  []    
  [strain_rr]
    family = MONOMIAL
    order = FIRST
  [] 
  [strain_tt]
    family = MONOMIAL
    order = FIRST
  [] 
[]

[Kernels]
  [./weight]
    type = BodyForce
    variable = disp_z
    value = -2.6e4 # this is density*gravity
  [../]
  [./mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  [../]
  [./flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  [../]
  [./mass_co2_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    use_displaced_mesh = false
    variable = sgas
  [../]
  [./flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    use_displaced_mesh = false
    variable = sgas
  [../]
  # [./energy_dot]
  #   type = PorousFlowEnergyTimeDerivative
  #   use_displaced_mesh = false
  #   variable = temp
  # [../]
  # [./advection]
  #   type = PorousFlowHeatAdvection
  #   use_displaced_mesh = false
  #   variable = temp
  # [../]
  # [./conduction]
  #   type = PorousFlowExponentialDecay
  #   use_displaced_mesh = false
  #   variable = temp
  #   reference = 358
  #   rate = rate
  # [../]
  [./grad_stress_r]
    type = StressDivergenceRZTensors
    # temperature = temp
    # thermal_eigenstrain_name = thermal_contribution
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  [../]
  [./grad_stress_z]
    type = StressDivergenceRZTensors
    variable = disp_z
    component = 1
    use_displaced_mesh = false
  []

  [./poro_r]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_r
    use_displaced_mesh = false
    component = 0
  [../]
  [./poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    use_displaced_mesh = false
    component = 1
  [../]    
[]

[AuxKernels]
  [./initial_stress1]
    type = SolutionAux
    solution = ini_soln
    from_variable = stress_rr
    execute_on = initial
    variable = stress_ini1 
  [../]
  [./initial_stress2]
    type = SolutionAux
    solution = ini_soln
    from_variable = stress_zz
    execute_on = initial
    variable = stress_ini2
  [../]  
  [./initial_stress3]
    type = SolutionAux
    solution = ini_soln
    from_variable = stress_tt
    execute_on = initial
    variable = stress_ini3
  [../]    
  [./rate]
    type = FunctionAux
    variable = rate
    execute_on = timestep_begin
    function = decay_rate
  [../]
  [./pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas
  [../]
  [./swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  [../]
  [./stress_rr]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_rr
    index_i = 0
    index_j = 0
  [../]
  [./stress_tt]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_tt
    index_i = 2
    index_j = 2
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 1
    index_j = 1
  [../]
  [strain_zz]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz
    index_i = 1
    index_j = 1
  []
  [strain_zz_0]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_zz_0
    index_i = 1
    index_j = 1
    execute_on = initial
  []  
  [strain_zz_2]
    type = ParsedAux
    variable = strain_zz_2
    args = 'strain_zz_0 strain_zz'
    function = 'strain_zz-strain_zz_0'#'9.869233e-13*perm_md'
    # execute_on = initial
  []    

  [pgas0]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas0
    execute_on = initial
  []  
  [pgas2]
    type = ParsedAux
    variable = pgas2
    args = 'pgas0 pgas'
    function = 'pgas-pgas0'#'9.869233e-13*perm_md'
    # execute_on = initial
  []   

  [strain_tt]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_tt
    index_i = 2
    index_j = 2
  []  
  [strain_rr]
    type = RankTwoAux
    rank_two_tensor = elastic_strain
    variable = strain_rr
    index_i = 0
    index_j = 0
  [] 
[]

[Functions]
  [ppic]
    type = ParsedFunction
    value = '1000*9.81*(0-y)'
  []
  [./rhogh]
    type = ParsedFunction
    value = '0*1000*9.81*(0-y)' # initial stress that should result from the weight force
  [../]    
  [./decay_rate]
# Eqn(26) of the first paper of LaForce et al.
# Ka * (rho C)_a = 10056886.914
# h = 11
    type = ParsedFunction
    value = 'sqrt(10056886.914/t)/11.0'
  [../]
  # [./insitu_pp]
  #   type = ParsedFunction
  #   value = if(t<1e0,18.3e6,18.3e6+t) #'18.3e6+t' #'min(t/100.0,1)*(-2.294001475)'  if(t<1e3,1e2,1e4)
  # [../]  
[]

[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater sgas disp_r'
    number_fluid_phases = 2
    number_fluid_components = 2
  [../]
  [./pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  [../]
[]

[Modules]
  [./FluidProperties]
    [./water]
      type = SimpleFluidProperties
      bulk_modulus = 2.1e9
      density0 = 970.0
      viscosity = 8.9e-4
      cv = 4149.0
      cp = 4149.0
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    [../]
    [./co2]
      type = SimpleFluidProperties
      bulk_modulus = 2.1e8
      density0 = 600
      viscosity = 1e-5
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    [../]
  [../]
[]

[Materials]
  # [./temperature]
  #   type = PorousFlowTemperature
  #   temperature = temp
  # [../]
  [./temperature]
    type = PorousFlowTemperature
    temperature = '45'
  [../]  
  [./ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc
  [../]
  # [./brineco2]
  #   type = PorousFlowFluidState
  #   gas_porepressure = 'pgas'
  #   z = 'zi'
  #   temperature_unit = Celsius
  #   xnacl = 'xnacl'
  #   capillary_pressure = pc
  #   fluid_state = fs
  # [../]  
  [./massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  [../]
  [./water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  [../]
  [./gas]
    type = PorousFlowSingleComponentFluid
    fp = co2
    phase = 1
  [../]
  [./porosity_reservoir]
    type = PorousFlowPorosityConst
    porosity = 0.2
  [../]
  # [./permeability_reservoir]
  #   type = PorousFlowPermeabilityConst
  #   permeability = '2e-12 0 0  0 2e-12 0  0 0 2e-12'
  # [../]

  [./permeability_aquifer]
    type = PorousFlowPermeabilityConst
    block = aquifer
    permeability = '2E-12 0 0   0 2E-12 0   0 0 2E-12'
  [../]
  [./permeability_caps]
    type = PorousFlowPermeabilityConst
    block = caps
    permeability = '1E-19 0 0   0 1E-19 0   0 0 1E-19'
  [../]

  [./relperm_liquid]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    phase = 0
    s_res = 0.200
    sum_s_res = 0.405
  [../]
  [./relperm_gas]
    type = PorousFlowRelativePermeabilityBC
    phase = 1
    s_res = 0.205
    sum_s_res = 0.405
    nw_phase = true
    lambda = 2
  [../]
  # [./thermal_conductivity_reservoir]
  #   type = PorousFlowThermalConductivityIdeal
  #   dry_thermal_conductivity = '0 0 0  0 1.320 0  0 0 0'
  #   wet_thermal_conductivity = '0 0 0  0 3.083 0  0 0 0'
  # [../]
  # [./internal_energy_reservoir]
  #   type = PorousFlowMatrixInternalEnergy
  #   specific_heat_capacity = 1100
  #   density = 2350.0
  # [../]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 6.0E9
    poissons_ratio = 0.2
  [../]
  [./strain]
    type = ComputeAxisymmetricRZSmallStrain
    eigenstrain_names = 'ini_stress'
  [../]
  [./ini_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = '0 0 0  0 0 0  0 0 0'
    initial_stress_aux = 'stress_ini1 stress_ini0 stress_ini0  stress_ini0 stress_ini2 stress_ini0  stress_ini0 stress_ini0 stress_ini3'
    eigenstrain_name = ini_stress
  [../]
  # [./thermal_contribution]
  #   type = ComputeThermalExpansionEigenstrain
  #   temperature = temp
  #   stress_free_temperature = 358
  #   thermal_expansion_coeff = 5E-6
  #   eigenstrain_name = thermal_contribution
  # [../]
  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    outputs = exodus
  [../]
  [./vol_strain]
    type = PorousFlowVolumetricStrain
    consistent_with_displaced_mesh = false
  [../]
[]

[BCs]
  # [./outer_pressure_fixed]
  #   type = DirichletBC
  #   boundary = right
  #   value = 18.3e6
  #   variable = pwater
  # [../]
  [top]
    type = DirichletBC
    variable = pwater
    value = 0
    boundary = 'top'
  []       
  [xmax_drained]
    type = FunctionDirichletBC
    variable = pwater
    boundary = 'right'
    function = ppic #pressure_BC 
  [] 
  [./outer_saturation_fixed]
    type = DirichletBC
    boundary = right
    value = 0.0
    variable = sgas
  [../]
  # [./outer_temp_fixed]
  #   type = DirichletBC
  #   boundary = right
  #   value = 358
  #   variable = temp
  # [../]
  [./fixed_outer_r]
    type = PresetBC
    variable = disp_r
    value = 0
    boundary = right
  [../]
  [./fixed_bottom_z]
    type = PresetBC
    variable = disp_z
    value = 0
    boundary = bottom
  [../]  
  [./co2_injection]
    type = PorousFlowSink
    boundary = 'injection_area' # left # 
    variable = sgas
    use_mobility = false
    use_relperm = false
    fluid_phase = 1
    flux_function = 0 #'min(t/100.0,1)*(-1.294001475)' # 5.0E5 T/year = 15.855 kg/s, over area of 2Pi*0.1*11
  [../]

  # [./constant_injection_porepressure]
  #   type = FunctionDirichletBC
  #   function = insitu_pp
  #   # type = DirichletBC
  #   variable = pwater
  #   # value = 19.e6
  #   boundary =  left #injection_area
  # [../]

  # [./cold_co2]
  #   type = PresetBC
  #   boundary = left
  #   variable = temp
  #   value = 294
  # [../]
  [./cavity_pressure_x]
    type = Pressure
    boundary = left
    variable = disp_r
    component = 0
    postprocessor = p_bh # note, this lags
    use_displaced_mesh = false
  [../]

[]

[Postprocessors]
  [./p_bh]
    type = PointValue
    variable = pwater
    point = '0.1 0 0'
    execute_on = timestep_begin
    use_displaced_mesh = false
  [../]
[]

# [VectorPostprocessors]
#   [./ptsuss]
#     type = LineValueSampler
#     use_displaced_mesh = false
#     start_point = '0.1 0 0'
#     end_point = '5000 0 0'
#     sort_by = x
#     num_points = 50000
#     outputs = csv
#     variable = 'pwater sgas disp_r stress_rr stress_tt'
#   [../]
# []

[Preconditioning]
  active = 'mumps'
  [./smp]
    type = SMP
    full = true
    #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E2       1E-5        500'
  [../]
  [./mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-9       1E-4       50'
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 1.5768e10
  automatic_scaling = True
  compute_scaling_once = False #True #

  #dtmax = 1e6
  [./TimeStepper]
    type = IterationAdaptiveDT
    dt = 1
    growth_factor = 2
  [../]
[]

[Outputs]
  print_linear_residuals = false
  # execute_on = 'initial timestep_end'
  # sync_times = '3600 86400 2.592E6 1.5768E8'
  perf_graph = true
  exodus = true
  # [./csv]
  #   type = CSV
  #   sync_only = true
  # [../]
[]
