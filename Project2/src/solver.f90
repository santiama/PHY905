program solver 
    use input
    use unit_tests

    character infile*40
    ! UNIT TESTS / 
    call init_tests()
    call max_test()
    call gram_test()
    ! / 

    call getarg(1,infile)
    call read_inputs(infile)
    call init()
    call solve()
    call savedata()
    call dealloc()

end program solver
