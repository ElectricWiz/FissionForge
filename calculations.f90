module reactor_simulation
  use geometry_module
  use one_group_diffusion_solver
  use output_module
  
  implicit none
  
  private
  public :: calculate_reactor
  
  real, parameter :: NU = 2.6
  real, parameter :: ETA = 2.2
  real, parameter :: SIGMA_F_MACRO = 0.6832 
  real, parameter :: SIGMA_A_MACRO = 0.8052
  real, parameter :: U235_DENSITY = 19.1 ! g/cm^3
  real, parameter :: U235_ATOMIC_MASS = 235.0439299 ! g/mol
  real, parameter :: AVOGADRO = 6.022E23 ! atoms/mol
  real, parameter :: SIGMA_H20_MACRO = 0.0197

contains

  subroutine calculate_reactor(reactor_geom, flux, k_eff, power)
    type(geometry), intent(in) :: reactor_geom
    real, allocatable, intent(out) :: flux(:), power(:)
    real, intent(out) :: k_eff
    real :: atomic_density, sigma_f, sigma_a
    
    call calculate_atomic_density(U235_DENSITY, U235_ATOMIC_MASS, AVOGADRO, atomic_density)
    call calculate_macroscopic_cross_sections(atomic_density, SIGMA_F_MICRO, SIGMA_A_MICRO, sigma_f, sigma_a)
    
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
  
  subroutine calculate_k_infty(eta, SIGMA_A_MACRO, SIGMA_H20_MACRO, k_infty)
    real, intent(in) :: eta, SIGMA_A_MACRO, SIGMA_H20_MACRO
    real, intent(out) :: k_infty

    k_infty = eta*(SIGMA_A_MACRO/(SIGMA_H20_MACRO+SIGMA_A_MACRO))
  end subroutine
end module reactor_simulation