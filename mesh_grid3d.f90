module meshgrid
    implicit none
    private
    public :: create_meshgrid

contains

    subroutine create_meshgrid(x, y, z, xx, yy, zz)
        real, intent(in) :: x(:)
        real, intent(in), optional :: y(:), z(:)
        real, allocatable, intent(out) :: xx(:,:,:), yy(:,:,:), zz(:,:,:)
        integer :: nx, ny, nz, i, j, k

        nx = size(x)
        ny = 1
        nz = 1

        if (present(y)) ny = size(y)
        if (present(z)) nz = size(z)

        allocate(xx(nx,ny,nz))
        if (ny > 1) allocate(yy(nx,ny,nz))
        if (nz > 1) allocate(zz(nx,ny,nz))

        ! 1D case
        if (ny == 1 .and. nz == 1) then
            xx(:,1,1) = x
            return
        end if

        ! 2D case
        if (nz == 1) then
            do j = 1, ny
                xx(:,j,1) = x
            end do
            do i = 1, nx
                yy(i,:,1) = y
            end do
            return
        end if

        ! 3D case
        do k = 1, nz
            do j = 1, ny
                xx(:,j,k) = x
            end do
        end do
        do k = 1, nz
            do i = 1, nx
                yy(i,:,k) = y
            end do
        end do
        do j = 1, ny
            do i = 1, nx
                zz(i,j,:) = z
            end do
        end do
    end subroutine create_meshgrid

end module meshgrid