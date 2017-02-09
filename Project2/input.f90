module input
    implicit none    
    integer(kind=4) :: Nmax, Nmin, j
    character(len=2) :: solver_type

contains
    subroutine read_inputs()
    !{ -----------------------------------------------------------------
    !  - Reads inputs from namelist file
    !  ------------------------------------------------------------------
    !}
        implicit none
        Nmax = 1
        Nmin = 1
        solver_type = "HO"
    end subroutine read_inputs

end module input
