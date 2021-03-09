subroutine monte_carlo(Q_mats,x0_mc,p0_mc)

    use global
    implicit none

    real*8, allocatable :: Q_new(:,:),x0_new(:,:),p0_new(:,:)
    !real*8, intent(in) :: P_mats(ndof,lb_m:ub_m)
    real*8, intent(inout) :: Q_mats(ndof,lb_m:ub_m),x0_mc(nstate,nb),p0_mc(nstate,nb)

    integer :: tot_mode,imode,imove,ibead
    integer :: init_nmc
    integer :: i,j,iup,idn,k,init_mc,l

    real*8 :: iacc_x,iatt_x,iacc_s,iacc_s1,iatt_s,iatt_s1,rnum,rand
    real*8 :: wt,wt1,wtt
    real*8 :: energy,energy_old,del_e
    complex*16 :: bigphase, bigphase1   

    complex*16 :: gamma
    complex*16 :: gamma1

   
    allocate(Q_new(ndof,lb_m:ub_m),x0_new(nstate,nb),p0_new(nstate,nb))

     !=========MONTE CARLO INTEGRATION RUN==========   
 
    !initialize all zeroes

      iatt_x=0.d0
      iacc_x=0.d0
      iatt_s=0.d0
      iatt_s1=0.d0
      iacc_s=0.d0
      iacc_s1=0.d0

  do init_mc=1,nmc

    !-----randomly pick a mode to move-----                                 
    call random_number(rand)
    tot_mode = int(rand*mats_mode)
    imode = int(tot_mode-1)   ! here "-1" refers the lowest mats modes, in this case as M=3 so lowest mat mode is "-1"
    !--------------------------------
    
    !randomly pick a mapping bead to move                                                                                                                                            
    call random_number(rand)
    ibead = int(rand*nb)+1

   !---- cyclic condition-----
    iup=ibead+1
    if(iup.gt.nb) iup=1
    idn=ibead-1
    if(idn.lt.1) idn=nb
  !---------------------------
    
    !move eletronic state or nuclear coord
    !either 1 or 2 ,pick 1 for matsubara modes, 2 for mapping DOF
    call random_number(rand)
    imove=int(rand*3)+1

    !check  boundary cond.
    if(imove.lt.1) imove=1
    if(imove.gt.3) imove=3

    !=======================================================
    !=======================================================
    !move individual modes, generating trial moves
    
    if(imove==1) then
    
      Q_new = Q_mats
    
      do j=1,ndof
        call random_number(rand)
        rnum = rand
        Q_new(j,imode) = Q_mats(j,imode) + (rnum-0.5d0)*step(1)  ! here need to check with gran as well
      
      end do

    call sampling_u0(Q_new,wtt)
    call prefactor(Q_new,x0_mc,p0_mc,gamma)

    energy = wtt - (1.0/(beta))*log(abs(gamma))

    iatt_x = iatt_x + 1.0

    elseif(imove==2) then

      x0_new = x0_mc

      do j = 1,nstate
        call random_number(rand)
        rnum = rand
        x0_new(j,ibead) = x0_mc(j,ibead)+(rnum-0.5d0)*step(j+1)
      enddo

      call sampling_elec(x0_new,p0_mc,wtt)
      call prefactor(Q_mats,x0_new,p0_mc,gamma)

      energy = wtt - (1.0/(beta))*log(abs(gamma))
      iatt_s=iatt_s+1.0
    
    else

      p0_new = p0_mc
                  
      do j=1,nstate
         call random_number(rand)
         rnum=rand
         p0_new(j,ibead) = p0_mc(j,ibead)+(rnum-0.5d0)*step(j+1)
      enddo

      call sampling_elec(x0_mc,p0_new,wtt)
      call prefactor(Q_mats,x0_mc,p0_new,gamma)

      energy = wtt - (1.0/(beta))*log(abs(gamma))
      iatt_s1=iatt_s1+1.0

    end if

    
    !checking against original weight

    if (imove==1) then
     
      call sampling_u0(Q_mats,wt1)
      call prefactor(Q_mats,x0_mc,p0_mc,gamma1)

      energy_old = wt1 - (1.0/(beta))*log(abs(gamma1))


    elseif (imove==2) then
     
      call sampling_elec(x0_mc,p0_mc,wt1)
      call prefactor(Q_mats,x0_mc,p0_mc,gamma1)

      energy_old = wt1 - (1.0/(beta))*log(abs(gamma1))

    else
   
      call sampling_elec(x0_mc,p0_mc,wt1)
      call prefactor(Q_mats,x0_mc,p0_mc,gamma1)


      energy_old = wt1 - (1.0/(beta))*log(abs(gamma1))
    

    end if

  ! calculating the energy difference
  del_e = (energy-energy_old)

  !accepting/rejecting step
  wt = dexp(-beta*del_e)

   call random_number(rand)
   rnum=rand

   if(rnum.lt.wt) then
      if(imove==1) then
         Q_mats(:,imode) = Q_new(:,imode)
         iacc_x=iacc_x+1.d0
         !write(*,*) int_mc,imode,"nuc"
      else if(imove==2) then
         x0_mc(:,ibead)=x0_new(:,ibead)
         iacc_s=iacc_s+1.d0
         !write(*,*) int_mc,ibead,"accept map position"
      else
         p0_mc(:,ibead)=p0_new(:,ibead)
         iacc_s1=iacc_s1+1.d0
         !write(*,*) int_mc,ibead,"accept map momentum"
      end if
        
   end if
        
  enddo 

  open(20,file='final_config.dat')
  write(20,*) Q_mats
  write(20,*) x0_mc
  write(20,*) p0_mc
  close(20)

end subroutine monte_carlo
!============================================
subroutine sampling_u0(Q_new,U0)

  use global
  implicit none
 
  real*8, intent(in) :: Q_new(ndof,lb_m:ub_m)
  real*8, intent(out) :: U0
 
  real*8 :: q(ndof,nb),pot
  integer :: l
 
   !------------- transform the modes onto bead representation ----------------
  call back_transform_matrix(Q_new,q)
 
  !--- compute the potential in bead representation------
 
 pot = 0.
 
 do l = 1,nb
  pot = pot + 0.5*k1*q(1,l)**2   !state-indepedent part is hard-coded here, one can do it more generally
 enddo 
 
 pot = pot/real(nb)
 
 U0 = pot
 
 !------------------------------------------------------
 
 end subroutine sampling_u0
!============================================================
!               ELECTRONIC COORD. SAMPLING                                                                                                 
!=====================================================================                                                                     
 subroutine sampling_elec(x,p,wtt)

  use global,only : nstate,nb,ndof,beta,betap

  implicit none

  real*8, intent(in) :: x(nstate,nb),p(nstate,nb)
  real*8, intent(out) :: wtt

  real*8 ::  hel(nstate,nstate)

  integer :: i,j,k
  real*8 :: sums
  wtt = 1.0d0

  !gaussian term from e-ic variables and prefactor term                                                                                    
  sums=0.0 
   
  do i = 1,nb
    do j=1,nstate
      !gaussian part                                                                                                                        
      sums=sums + (x(j,i)*x(j,i) + p(j,i)*p(j,i)) 
    end do
  end do

  wtt=(1.0/(beta))*sums
 
end subroutine sampling_elec
!==============================================================
subroutine prefactor(Q_mats,x0_mc,p0_mc,theta)

  use global
  implicit none

  real*8, intent(in) :: Q_mats(ndof,lb_m:ub_m),x0_mc(nstate,nb),p0_mc(nstate,nb)
  complex*16, intent(out) :: theta

  
  integer :: i,iup,idn,j,k,l
  real*8 :: identity(nstate,nstate)
  real*8 :: coeff1(nstate,nstate)
  real*8 :: inv1(nstate,nstate)
  real*8 :: hel(nstate,nstate)
  complex*16 :: weight(nstate,nstate)
  complex*16 :: cmat(nstate,nstate)
  complex*16 :: full_mat(nstate,nstate)
  real*8 :: q(ndof,nb)

  !definition of initial gamma matrix--------
  cmat(1,1) = cmplx(1.,0.)
  cmat(1,2) = cmplx(0.,0.)
  cmat(2,1) = cmat(1,2)
  cmat(2,2) = cmat(1,1)
  !---------------------------------------

  full_mat(1,1) = cmplx(1.,0.)
  full_mat(1,2) = cmplx(0.,0.)
  full_mat(2,1) = full_mat(1,2)
  full_mat(2,2) = full_mat(1,1)

  !---- transform back to bead coordinate
  call back_transform_matrix(Q_mats,q)

  do i=1,nb

    iup=i+1
    if(iup.gt.nb) iup=1
    idn=i-1
    if(idn.lt.1) idn=nb
 
    call potential(q(:,i),hel)
    call coeff(hel,coeff1)
    inv1=coeff1

     !initialization of weight
   weight=cmplx(0.,0.)
 
   do j=1,nstate
     do k=1,nstate
        !for prefactor
        weight(j,k) = ((x0_mc(j,i)+eye*p0_mc(j,i))*(x0_mc(k,i)-eye*p0_mc(k,i)))&
             -(0.5d0*del(j,k))
     enddo
  enddo

  full_mat = MATMUL(weight,inv1)

!full prefactor
  
 cmat = MATMUL(cmat,full_mat)

enddo

!----- trace----------
theta=cmplx(0.,0.)

do l = 1,nstate
   theta = theta + cmat(l,l)
enddo
!------------------------

end subroutine prefactor
!==========================================================================