program reactor_simulation
    use reactor_modules
    implicit none

    integer :: geometry_type, calculation_type
    logical :: continue_simulation

    integer :: nx
    integer :: ny
    integer :: nz
    
    real :: radius
    real :: fill_value
    call initialize()

    do
        call display_menu()
        read(*, *) dimensions

        select case (dimensions)
        case (1)  
            print *,"How many points for the X axis?"
            read (*,*) nx
    
            ny = 1
            nz = 1 
            call create_meshgrid(nx, ny, nz, meshgrid)
        case (2)
            print *,"How many points for the X axis?"
            read (*,*) nx
            print *,"How many points for the Y axis?"
            read (*,*) ny
    
            nz = 1 
            call create_meshgrid(nx, ny, nz, meshgrid)
            select case (geometry_type)
            case (1)
                print *,"Radius Size"
                read (*,*) radius
                print *,"value at boundaries (DIRICHLET CONDITIONS)"
                read (*,*) fill_value

                call partition_circle(meshgrid, radius, fill_value, interior_mask)
            case (2)
                print *,"Width length"
                read (*,*) width
                print *,"Height length"
                read (*,*) height
                print *,"value at boundaries (DIRICHLET CONDITIONS)"
                read (*,*) fill_value

                
                call partition_rectangle(meshgrid, width, height, fill_value_value, boundary_mask, status)
        case (3)
                call calculate_spherical_mesh(sizeMesh)

            case default
                print *, "Invalid geometry type"
                cycle
            end select

            ! Solve one-group diffusion equation
            ! call solve_one_group_diffusion()

            ! Perform parametric studies
            ! call perform_parametric_studies()

            ! Optimize reactor design
            ! call optimize_reactor_design()

            ! Generate output and plots
            ! call generate_output()
            ! call plot_results()

            ! Save results to file
            ! call save_results()

        case (2) ! Load previous calculation
            ! call load_previous_calculation()

        case (3) ! Exit program
            exit

        case default
            print *, "Invalid choice. Please try again."
        end select

        ! Ask user if they want to continue
        print *, "Do you want to perform another calculation? (Y/N)"
        read(*, *) continue_simulation
        if (.not. continue_simulation) exit
    end do

    ! Clean up and finalize
    call finalize()

end program reactor_simulation