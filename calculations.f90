module reactor_simulation
    use geometry_module
    use one_group_diffusion_solver
    use output_module
    
    implicit none
    
    private
    public :: calculate_reactor
    
    real, parameter :: NU = 2.6
    real, parameter :: SIGMA_F_MICRO = 1.4
    real, parameter :: SIGMA_A_MICRO = 1.65
    real, parameter :: U235_DENSITY = 19.1 ! g/cm^3
    real, parameter :: U235_ATOMIC_MASS = 235.0439299 ! g/mol
    real, parameter :: AVOGADRO = 6.022E23 ! atoms/mol
    real :: ATOMIC_DENSITY
    real :: SIGMA_F, SIGMA_A
    
  contains
  
    subroutine calculate_reactor(reactor_geom, flux, k_eff, power)
      type(geometry), intent(in) :: reactor_geom
      real, allocatable, intent(out) :: flux(:), power(:)
      real, intent(out) :: k_eff
      
      call calculate_atomic_density()
      call calculate_macroscopic_cross_sections()
      
      call solve_diffusion(reactor_geom, NU, SIGMA_F, SIGMA_A, flux, k_eff)
      call calculate_power(flux, SIGMA_F, power)
      call output_results(flux, power, k_eff)
      
    end subroutine calculate_reactor
    
    subroutine calculate_atomic_density()
      ! Calculate atomic density from material density, atomic mass, and Avogadro's number
      ATOMIC_DENSITY = (U235_DENSITY * AVOGADRO) / U235_ATOMIC_MASS
    end subroutine calculate_atomic_density
    
    subroutine calculate_macroscopic_cross_sections()
      ! Calculate macroscopic cross-sections from microscopic and atomic density
      SIGMA_F = ATOMIC_DENSITY * SIGMA_F_MICRO
      SIGMA_A = ATOMIC_DENSITY * SIGMA_A_MICRO
    end subroutine calculate_macroscopic_cross_sections
    
    ! Other subroutines for solving diffusion equation, calculating power, output etc 
    
  end module reactor_simulation