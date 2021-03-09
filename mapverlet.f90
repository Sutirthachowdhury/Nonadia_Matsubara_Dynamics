subroutine mapverlet(hel,x0,p0)
  use global
  implicit none

  real*8 xrhs(nstate),prhs(nstate)
  real*8 force(nstate)
  real*8 xsumdum,psumdum,ddt

  real*8  hel(nstate,nstate)
  real*8  x0(nstate),p0(nstate)

  integer i,n,m,nlit

  nlit=10
  ddt=0.5*dt/nlit 
  
  do i=1,nlit

     do n=1,nstate
        xrhs(n)=hel(n,n)*p0(n)
        prhs(n)=-hel(n,n)*x0(n)
     end do
  
     do n=1,nstate
        xsumdum=0.
        psumdum=0.
        do m=1,nstate
           if(m.ne.n) then
              xsumdum=xsumdum+hel(n,m)*p0(m)
              psumdum=psumdum+hel(n,m)*x0(m)
           end if
        end do

        xrhs(n)=xrhs(n)+xsumdum
        prhs(n)=prhs(n)-psumdum
     end do

!     Generate current second derivatives

     do n=1,nstate
        force(n)=hel(n,n)*prhs(n)
     end do


     do n=1,nstate
        xsumdum=0.
        do m=1,nstate
           if(m.ne.n) then
              xsumdum=xsumdum+hel(n,m)*prhs(m)
           end if
        end do
        force(n)=force(n)+xsumdum
     end do

! other mapping force from gaussian overlap
!     do n=1,nstate
!        xrhs(n,i)=xrhs(n,i)+p0(n,i)*real(nb/beta)
!        prhs(n,i)=prhs(n,i)-x0(n,i)*real(nb/beta)
!        force(n,i)=force(n,i)-x0(n,i)*(nb/beta)**2
!     end do

!     Advance Ver step

     do n=1,nstate
        x0(n)=x0(n)+xrhs(n)*ddt+0.5*force(n)*ddt*ddt
        p0(n)=p0(n)+0.5*prhs(n)*ddt
     end do

!     Compute new first derivatives

     do n=1,nstate
        prhs(n)=-hel(n,n)*x0(n)
     end do
  
     do n=1,nstate
        psumdum=0.
        do m=1,nstate
           if(m.ne.n) then
              psumdum=psumdum+hel(n,m)*x0(m)
           end if
        end do
        prhs(n)=prhs(n)-psumdum
     end do

! other mapping force from gaussian   
!     do n=1,nstate
!        prhs(n,i)=prhs(n,i)-x0(n,i)*real(nb/beta)
!        force(n,i)=force(n,i)-x0(n,i)*(nb/beta)**2
!     end do

!    Advance let step

     do n=1,nstate
        p0(n)=p0(n)+0.5*prhs(n)*ddt
     end do

   ! end loop of the mapping propagation
  end do

  return
end subroutine mapverlet

