subroutine initialize_mapp(x0_mc,p0_mc)

    use global

    real*8, intent(out) :: x0_mc(nstate,nb),p0_mc(nstate,nb)

    real*8 :: gran,sigma

    integer :: i,j,k

    sigma = sqrt(0.5)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Initialize positions 
    do j = 1,nstate
        do i=1,nb
            x0_mc(j,i) = gran()*sigma
            p0_mc(j,i) = gran()*sigma
        enddo
    enddo   

end subroutine initialize_mapp