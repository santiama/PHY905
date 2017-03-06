program solver 
    use input
    character infile*40

    call getarg(1,infile)
    call read_inputs(infile)
    call init()
    call solve()
    call savedata()
    call dealloc()

end program solver
