module 3d_boundary
    implicit none
    private
    public :: partition_sphere, partition_cuboid

    contains 

    subroutine partition_sphere(mesh, radius, fill_value, interior_mask)
        real, intent(inout) :: mesh(:,:,:)  ! 3D mesh
        integer, intent(in) :: radius       ! Radius in mesh units
        real, intent(in) :: fill_value      ! Value to fill inside the sphere
        logical, intent(out) :: interior_mask(:,:,:)
        
        integer :: Nx, Ny, Nz, i, j, k, center_x, center_y, center_z
        real :: distance, dx, dy, dz, mesh_diagonal
        
        Nx = size(mesh, 1)
        Ny = size(mesh, 2)
        Nz = size(mesh, 3)
        
        ! Calculate center indices
        center_x = (Nx + 1) / 2
        center_y = (Ny + 1) / 2
        center_z = (Nz + 1) / 2
        
        ! Calculate mesh spacing (assuming uniform spacing)
        dx = 2.0 / (Nx - 1)
        dy = 2.0 / (Ny - 1)
        dz = 2.0 / (Nz - 1)
        
        ! Calculate diagonal of a mesh cell
        mesh_diagonal = sqrt(dx**2 + dy**2 + dz**2)
        
        ! Determine points inside the sphere and update mesh
        do k = 1, Nz
            do j = 1, Ny
                do i = 1, Nx
                    distance = sqrt(real((i - center_x)**2 + (j - center_y)**2 + (k - center_z)**2))
                    if ((distance <= real(radius) + 0.5*mesh_diagonal) .and. &
                        (distance <= real(radius) - 0.5*mesh_diagonal)) then
                        mesh(i,j,k) = fill_value
                        interior_mask(i,j,k) = .true.
                    else
                        interior_mask(i,j,k) = .false.
                    end if
                end do
            end do
        end do
    end subroutine partition_sphere
end module 3d_boundary