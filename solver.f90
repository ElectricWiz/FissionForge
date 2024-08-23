module helmholtz_solver
  use calculations

  implicit none
  private
  public :: initialize, solve

  integer :: dim  ! Problem dimension (1, 2, or 3)
  integer :: nx, ny, nz  ! Grid sizes
  real, parameter :: k = 1.0  ! Wavenumber
  real, allocatable :: u(:,:,:), f(:,:,:)  ! Solution and source term
  real, allocatable :: mesh(:,:,:)  ! Mesh grid

contains

  subroutine gauss_seidel_iteration()
    integer :: i, j, l
    real :: h_x, h_y, h_z
    
    select case (dim)
      case (1)
        do i = 2, nx-1
          h_x = mesh(i+1,1,1) - mesh(i-1,1,1)

          mesh(i,1,1) = -mesh(i-1,1,1) + (2 - BUCKLING_SQUARED**h_x^2)*mesh(i,1,1) - mesh(i+1,1,1)
        end do
      case (2)
        do j = 2, ny-1
          do i = 2, nx-1
            h_x = mesh(i+1,j,1) - mesh(i-1,j,1)
            h_y = mesh(i,j+1,1) - mesh(i,j-1,1)

            mesh(i,j,1) = mesh(i+1,j,1) + mesh(i-1,j,1) + mesh(i,j+1,1) + mesh(i,j-1,1) - (4 - BUCKLING_SQUARED*SQRT(h_x**2+h_y**2))*mesh(i,j,1)
          end do
        end do
      case (3)
        do l = 2, nz-1
          do j = 2, ny-1
            do i = 2, nx-1
              h_x = mesh(i+1,j,l) - mesh(i-1,j,l)
              h_y = mesh(i,j+1,l) - mesh(i,j-1,l)
              h_z = mesh(i,j,l+1) - mesh(i,j,l-1)
              
              mesh(i,j,l) = (6-BUCKLING_SQUARED*SQRT(h_x**2 + h_y**2 + h_z**2))*mesh(i,j,l) - (mesh(i+1,j,k)+mesh(i-1,j,k)+mesh(i,j-1,k)+mesh(i,j+1,k)+mesh(i,j,k-1)+mesh(i,j,k+1))
            end do
          end do
        end do
    end select
  end subroutine gauss_seidel_iteration

  subroutine solve()
    integer :: iter, max_iter
    real :: tolerance, error
    
    max_iter = 100000
    tolerance = 1.0e-8
    
    do iter = 1, max_iter
      call gauss_seidel_iteration()
      ! Calculate error
      ! Check for convergence
      if (error < tolerance) exit
    end do
  end subroutine solve

end module helmholtz_solver
