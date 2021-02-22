# 2D RZ hydro-mechanical model--Initialization with gravity equilibrium of stress strain and pore pressure
# based on example "https://github.com/idaholab/falcon/blob/devel/tests/THM_injection/thm_steady.i"

[Mesh]
 type = MeshGeneratorMesh
 [./cartesian]
   type = CartesianMeshGenerator
   dim = 2
   dx = '0.2 15000'
   ix = '2 150'
   dy = '15000'
   iy = '150'
 [../]
[] 

[MeshModifiers]
  [shift_down]
    type = Transform
    transform = TRANSLATE
	  vector_value = '0.1 -15000 0'
  []
  [aquifer]
    type = SubdomainBoundingBox
    block_id = 1
    bottom_left = '0.1 -1500 0' 
    top_right = '30000.1 -1000 0'
	depends_on = shift_down
  []
[]

############################################################
[Problem]
  coord_type = RZ
[]

[GlobalParams]
  displacements = 'disp_r disp_z'
  PorousFlowDictator = dictator
  gravity = '0 -9.81 0'
  biot_coefficient = 1.0
[]
############################################################
[Variables]
  [./pwater]
  
  [../]
  [./disp_r]
  [../]
  [./disp_z]
  [../]  
[]
############################################################
[ICs]
  [pwater]
    type = FunctionIC
    function = ppic
    variable = pwater
  []
[]
############################################################

[Functions]
  [ppic]
    type = ParsedFunction
    value = '101325-1000*9.81*(0-y)'
  []

  [./weight_fcn]
    type = ParsedFunction
    vars = 'g rho0 biot'
    vals = '9.81 1000 1.0'
    value = 'biot*(101325-rho0*g*y)+2386.0*g*y'
  [../]
[]
############################################################

[Kernels]
  # [./mass_water_dot]
  #   type = PorousFlowMassTimeDerivative
  #   fluid_component = 0
  #   # use_displaced_mesh = false
  #   variable = pwater
  # [../]
  [./flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    # use_displaced_mesh = false
    variable = pwater
    gravity = '0 -9.81 0'
  [../]
  [./grad_stress_r]
    type = StressDivergenceRZTensors
    variable = disp_r
    # use_displaced_mesh = false
    component = 0
  [../]
  [./grad_stress_z]
    type = StressDivergenceRZTensors
    variable = disp_z
    component = 1
    # use_displaced_mesh = false
  []
  [./gravity_z]
    type = Gravity
    variable = disp_z
    value = -9.81
  [../]  
  [./poro_r]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_r
    biot_coefficient = 1.0
    # use_displaced_mesh = false
    component = 0
  [../]
  [./poro_z]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_z
    biot_coefficient = 1.0
    # use_displaced_mesh = false
    component = 1
  [../]  
[]
############################################################

[BCs]
  [top]
    type = FunctionDirichletBC
    variable = pwater
    boundary = 'top'
    # use_displaced_mesh = false
    function = '101325-z*9810'
  []    
     
  [xmax_drained]
    type = FunctionDirichletBC
    variable = pwater
    boundary = 'right left'
    function = '101325-z*9810'  
    # use_displaced_mesh = false 
  []  

  [./fixed_outer_r]
    type = PresetBC
    variable = disp_r
    value = 0
    boundary = right
    # use_displaced_mesh = false
  [../]

  [./fixed_outer_l]
    type = PresetBC
    variable = disp_r
    value = 0
    boundary = 'left' 
    # use_displaced_mesh = false
  [../]

  [./fixed_bottom_z]
    type = PresetBC
    variable = disp_z
    value = 0
    boundary = bottom
    # use_displaced_mesh = false
  [../]  

  [./top_z]
    type = Pressure
    variable = disp_z
    component = 1
    function = '-2386.0*9.81*z'
    use_displaced_mesh = false
    boundary = 'top'
  [../]  
[]

############################################################
[AuxVariables]
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
############################################################

[AuxKernels]
  [./stress_rr]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_rr
    index_i = 0
    index_j = 0
  [../]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 1
    index_j = 1
  [../]
  [./stress_tt]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_tt
    index_i = 2
    index_j = 2
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
    function = 'strain_zz-strain_zz_0' 
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

############################################################
[UserObjects]
  [./dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater disp_r disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
  [./pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  [../]
[]
############################################################
[Modules]
  [./FluidProperties]
    [./water]
      type = SimpleFluidProperties
      bulk_modulus = 2.1e9
      density0 = 1000
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
############################################################
[Materials]
  [./temperature]
    type = PorousFlowTemperature
    temperature = '45'
  [../]  
  [./massfrac]
    type = PorousFlowMassFraction
    # at_nodes = true
  [../]

  [./water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
    # at_nodes = true
  [../]

  [./porosity_reservoir]
    type = PorousFlowPorosityConst
    # at_nodes = true
    porosity = 0.2
  [../]

  [./permeability_aquifer]
    type = PorousFlowPermeabilityConst
    # block = aquifer
    permeability = '2E-12 0 0   0 2E-12 0   0 0 2E-12'
  [../]

  [./relperm]
    type = PorousFlowRelativePermeabilityCorey
    # at_nodes = true
    n = 1
    phase = 0
  [../]

  [./ppss]
    type = PorousFlow1PhaseP
    # at_nodes = true
    porepressure = pwater
    capillary_pressure = pc
  [../]

  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 6.0E9
    poissons_ratio = 0.2
  [../]
  [./strain]
    type = ComputeAxisymmetricRZSmallStrain
    eigenstrain_names = 'ini_stress'
    # # use_displaced_mesh = false
  [../]
  [./ini_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'weight_fcn 0 0  0 weight_fcn 0  0 0 weight_fcn' # '-12.8E6 0 0  0 -12.8E6 0  0 0 -12.8E6'
    eigenstrain_name = ini_stress
    # # use_displaced_mesh = false
  [../]

  [./stress]
    type = ComputeLinearElasticStress
  [../]
  [./eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
    # at_nodes = true
    outputs = exodus
  [../]
  [./density]
    type = GenericConstantMaterial
    prop_names = density
    prop_values = 2386.0 # = (1-0.1)*2540 + 0.1*999.526
  [../]  
[]
############################################################

[Preconditioning]
  active = superlu #'superlu' #
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
    petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-9       1E-5       50'
  [../]
  [./superlu]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package'
    petsc_options_value = 'gmres lu superlu_dist'
  [../]  
[]
############################################################

[Executioner]
  type = Steady
  solve_type = 'NEWTON' # default = PJFNK | NEWTON
  automatic_scaling = True
  compute_scaling_once = False #True #
  l_max_its  = 50
  l_tol      = 1e-4
  nl_max_its = 500
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
[]

############################################################
[Outputs]
  exodus         = true
  [./console]
    type = Console
    output_linear = true
    output_nonlinear = true
  [../]
[]
