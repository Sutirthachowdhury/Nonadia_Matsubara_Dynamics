subroutine potential(r,hel)

    use global,only : nstate,ndof,nb,k1,c,delta
   
    implicit none
  
    real*8, intent(in) :: r(ndof)
    real*8, intent(out) :: hel(nstate,nstate)
  
    integer :: i
  
    hel(1,1) = 0.0d0; hel(2,2) = 0.0d0 
  
       !*** v11
    do i = 1 , ndof
       hel(1,1) = hel(1,1) + k1*r(i) + c
       !*** v22
       hel(2,2) = hel(2,2) - k1*r(i) 
    end do
  
    hel(1,2) = delta
    hel(2,1) = hel(1,2)

    return
  
  end subroutine potential
  