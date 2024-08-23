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
        
        dx = 1.0 / real(nx - 1)
        dy = 1.0 / real(ny - 1)
        dz = 1.0 / real(nz - 1)
        
        ! Generate x coordinates
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    meshgrid(i,j,k,1) = real(i-1) * dx
                end do
            end do
        end do
        
        ! Generate y coordinates
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    meshgrid(i,j,k,2) = real(j-1) * dy
                end do
            end do
        end do
        
        ! Generate z coordinates
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    meshgrid(i,j,k,3) = real(k-1) * dz
                end do
            end do
        end do
    end subroutine create_meshgrid

end module meshgrid_3d

program test_meshgrid
    use meshgrid_3d
    implicit none
    
    integer :: nx, ny, nz
    real, allocatable :: result(:,:,:,:)
    integer :: i, j, k

    ! Get input from user
    print *, "Enter nx, ny, nz:"
    read *, nx, ny, nz

    ! Create 3D mesh
    call create_meshgrid(nx, ny, nz, result)

    ! Print a few values to verify
    print *, "Sample values:"
    print *, "At (1,1,1): ", result(1,1,1,:)
    print *, "At (nx,ny,nz): ", result(nx,ny,nz,:)

    ! Optional: Print full grid (be cautious with large grids)
    if (nx*ny*nz <= 1000) then  ! Only print if total points <= 1000
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    print *, "At (", i, ",", j, ",", k, "): ", result(i,j,k,:)
                end do
            end do
        end do
    else
        print *, "Grid too large to print all values."
    end if

    deallocate(result)

end program test_meshgrid