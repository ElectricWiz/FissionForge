program reactor_simulation
    use reactor_modules
    implicit none

    integer :: geometry_type, calculation_type
    logical :: continue_simulation

    integer :: sizeMesh
    
    print *,"How many points for the total mesh?"
    read (*,*) sizeMesh

    call initialize()

    do
        call display_menu()
        read(*, *) calculation_type

        select case (calculation_type)
        case (1) ! New calculation
            call get_reactor_geometry(geometry_type)

            ! Perform calculation based on geometry
            select case (geometry_type)
            case (1) 
                 call calculate_slab_mesh(sizeMesh)
            case (2)  
                 call calculate_cylindrical_mesh(sizeMesh)
            case (3)
                call calculate_spherical_mesh(sizeMesh)
            case (4) 
                call calculate_reflected_mesh(sizeMesh)
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