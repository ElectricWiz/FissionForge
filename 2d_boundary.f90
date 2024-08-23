module geometry_partitioner
    implicit none
    private
    public :: partition_circle, partition_square, partition_rectangle
  
    contains
    
    subroutine partition_circle(mesh, radius, fill_value, interior_mask)
        real, intent(inout) :: mesh(:,:)  ! 2D mesh
        integer, intent(in) :: radius    ! Radius in mesh units
        real, intent(in) :: fill_value   ! Value to fill inside the circle
        logical, intent(out) :: interior_mask(:,:)
        
        integer :: Nx, Ny, i, j, center_x, center_y
        real :: distance
        
        Nx = size(mesh, 1)
        Ny = size(mesh, 2)
        
        ! Check if Nx = Ny
        if (Nx /= Ny) then
            print *, "Error: Mesh must be square (Nx = Ny)"
            return
        end if
        
        ! Calculate center indices
        center_x = (Nx + 1) / 2
        center_y = (Ny + 1) / 2
        
        ! Determine points inside the circle and update mesh
        do j = 1, Ny
            do i = 1, Nx
                distance = sqrt(real((i - center_x)**2 + (j - center_y)**2))
                if ((distance <= real(radius)+1.5/(SQRT(Nx**2+Ny**2)-1)) .and. (distance <= real(radius)-1.5/(SQRT(Nx**2+Ny**2)-1) )) then
                mesh(i,j) = fill_value
                interior_mask(i,j) = .true.
                else
                interior_mask(i,j) = .false.
                end if
            end do
        end do
    end subroutine partition_circle
 
    subroutine partition_rectangle(mesh, width, height, boundary_value, boundary_mask, status)
        real, intent(inout) :: mesh(:,:)  ! 2D mesh
        integer, intent(in) :: width, height  ! Rectangle dimensions
        real, intent(in) :: boundary_value  ! Value to fill in the boundary layer
        logical, intent(out) :: boundary_mask(:,:)  ! True for boundary points
        integer, intent(out) :: status  ! Status code: 0 for success, non-zero for error
        
        integer :: Nx, Ny, i, j
        integer :: left, right, top, bottom
        integer :: center_x, center_y
        
        ! Initialize status to success
        status = 0
        
        Nx = size(mesh, 1)
        Ny = size(mesh, 2)
        
        ! Calculate center of the mesh
        center_x = (Nx + 1) / 2
        center_y = (Ny + 1) / 2
        
        ! Error checking
        if (width <= 2 .or. height <= 2) then
            print *, "Error: dimensions must be greater than 2"
            status = 1
            return
        end if
        
        ! Calculate rectangle boundaries
        left = center_x - width / 2
        right = center_x + (width - 1) / 2
        top = center_y - height / 2
        bottom = center_y + (height - 1) / 2
        
        ! Check if rectangle fits within the mesh
        if (left < 1 .or. right > Nx .or. top < 1 .or. bottom > Ny) then
            print *, "Error: rectangle extends outside the mesh"
            status = 2
            return
        end if
        
        ! Initialize boundary_mask
        boundary_mask = .false.
        
        ! Fill boundary layer
        do j = top, bottom
            do i = left, right
                if (i == left .or. i == right .or. j == top .or. j == bottom) then
                    mesh(i,j) = boundary_value
                    boundary_mask(i,j) = .true.
                end if
            end do
        end do
    end subroutine partition_rectangle
    
end module geometry_partitioner