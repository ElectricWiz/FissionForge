module reactor_simulation
  use geometry_module
  use one_group_diffusion_solver
  use output_module
  
  implicit none
  
  private
  public :: calculate_reactor
  
  real, parameter :: NU = 2.6
  real, parameter :: ETA = 2.2
  real, parameter :: SIGMA_F_MACRO_U235 = 0.6832 
  real, parameter :: SIGMA_A_MACRO_U235 = 0.8052
  real, parameter :: SIGMA_A_MACRO_H20 = 0.0197
  real, parameter :: SIGMA_TR = 3.318

  real, parameter :: SIGMA_A = SIGMA_A_MACRO_U235 + SIGMA_A_MACRO_H20
  real, parameter :: DIFF_COEFF = 1/(3*SIGMA_TR)
  real, parameter :: DIFF_AREA = DIFF_AREA/SIGMA_A

  real, parameter :: K_INFTY = ETA*(SIGMA_A_MACRO_U235/SIGMA_A)
  real, parameter :: BUCKLING_SQUARED = (K_INFTY - 1)/DIFF_AREA

contains

  subroutine calculate_reactor(reactor_geom, flux, k_eff, power)
    type(geometry), intent(in) :: reactor_geom
    real, allocatable, intent(out) :: flux(:), power(:)
    real, intent(out) :: k_eff
    real :: atomic_density, sigma_f, sigma_a
    
    call solve_diffusion(reactor_geom, NU, sigma_f, sigma_a, flux, k_eff)
    call calculate_power(flux, sigma_f, power)
    call output_results(flux, power, k_eff)
    
  end subroutine calculate_reactor
  
  subroutine calculate_atomic_density(density, atomic_mass, avogadro, atomic_density)
    real, intent(in) :: density, atomic_mass, avogadro
    real, intent(out) :: atomic_density
    
    atomic_density = (density * avogadro) / atomic_mass
  end subroutine calculate_atomic_density
  
  ! Other subroutines for solving diffusion equation, calculating power, output etc 
end module reactor_simulation