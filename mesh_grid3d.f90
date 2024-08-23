module meshgrid_3d
    implicit none
    private
    public :: create_meshgrid

contains

    subroutine create_meshgrid(nx, ny, nz, meshgrid)
        integer, intent(in) :: nx, ny, nz
        real, allocatable, intent(out) :: meshgrid(:,:,:,:)  ! 4D array to hold x, y, z coordinates
        integer :: i, j, k
        real :: dx, dy, dz
        
        allocate(meshgrid(nx, ny, nz, 3))
        meshgrid = 0.0  ! Initialize all elements to zero
        
        dx = 2.0 / real(nx - 1)
        dy = 2.0 / real(ny - 1)
        dz = 2.0 / real(nz - 1)
        
        ! Generate x coordinates
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    meshgrid(i,j,k,1) = real(i-1) * dx - 1.0
                end do
            end do
        end do
        
        ! Generate y coordinates
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    meshgrid(i,j,k,2) = real(j-1) * dy - 1.0
                end do
            end do
        end do
        
        ! Generate z coordinates
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    meshgrid(i,j,k,3) = real(k-1) * dz - 1.0
                end do
            end do
        end do
    end subroutine create_meshgrid

end module meshgrid_3d
