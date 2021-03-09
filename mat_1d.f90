program mat_1d

    use global
    implicit none

    integer :: istep,iup,idn,isample,ieq     
    integer,allocatable :: seed(:)
    integer seed_dimension,time


    integer i,j,k,l,m,n,ncount

    !======= premitive bead variables=============

    !real*8, allocatable :: xe_mc(:,:),xe(:,:) ! nuclear bead positions 
    !real*8, allocatable :: ve(:,:) ! nuclear beads momentum
    !=======================================


    real*8, allocatable:: Q_mats(:,:),P_mats(:,:),f_mats(:,:) !"Matsubara" positions, momentum and force
    real*8, allocatable:: x0_mc(:,:),p0_mc(:,:) ! mapping positons and momentum (in bead rep.)
    real*8, allocatable:: Q_mats_dyn(:,:),P_mats_dyn(:,:),f_mats_dyn(:,:),x0(:,:),p0(:,:) ! this variables will use in dynamics
    real*8, allocatable :: matsfreq(:)
    real*8, allocatable :: corr_xx(:) ! this is RR auto-correlation

    real*8 :: q(ndof,nb)

    real*8 :: dp,gran,rand,theta 
    real*8 :: corr_xx0 ! initial value of RR correlation function

    complex*16 :: pref ! this is the weighting term from non-adiabatic part
    complex*16 :: phase ! this is the real part of e^{i*beta*theta}

    real*8 :: partition ! weighting the cosine term

    ! getting the inputs
    call parameter_input
 
    !-------- random seed-------------------
    CALL RANDOM_SEED(size=seed_dimension)
    ALLOCATE (seed(seed_dimension))
    do i=1,seed_dimension
      seed(i) =time()+3*i-1 + ourseed
    end do
    CALL RANDOM_SEED(PUT=seed)
    !-----------------------------------------

    !allocate(xe_mc(ndof,nb),xe(ndof,nb),ve(ndof,nb))
    allocate(Q_mats(ndof,lb_m:ub_m),P_mats(ndof,lb_m:ub_m),f_mats(ndof,lb_m:ub_m))
    allocate(Q_mats_dyn(ndof,lb_m:ub_m),P_mats_dyn(ndof,lb_m:ub_m),f_mats_dyn(ndof,lb_m:ub_m))
    allocate(x0(nstate,nb),p0(nstate,nb))
    allocate(matsfreq(lb_m:ub_m))
    allocate(x0_mc(nstate,nb),p0_mc(nstate,nb))
    
    allocate(corr_xx(nstep))

    !--- initialization of correlation function 
    corr_xx = 0.0
    corr_xx0 = 0.0
    !----------------------------

    !------- initializing the partition function----

    partition = 0.0

    !---- initialize the matsubara pos,momentum, transform matrix coeff and mapping variables----
    call initialize_mat_pos(Q_mats) ! for matsubara positions
    call initialize_mapp(x0_mc,p0_mc) ! for mapping DOF
    !call momentum_sample(P_mats) ! for initializing momentum from Maxwell-Boltzmann distribution
    !call init_pot() ! initialize potential variables (coeff for infinity N limit tranformation)
    !===========================================================

       
    !----- run thermalization trajectory---
    do j = 1,neq
      !call momentum_sample(P_mats)
      call monte_carlo(Q_mats,x0_mc,p0_mc)
  
    enddo 
    !----------------------------------------


    do j = 1,ntraj


    
      call monte_carlo(Q_mats,x0_mc,p0_mc)
       
       !store sample configuration
       Q_mats_dyn = Q_mats
       x0=x0_mc
       p0=p0_mc
  
      call momentum_sample(P_mats)
      !------------- calculation of prefactor part----------

      call estimator(Q_mats_dyn,x0,p0,pref)

      call bigphase_sub(Q_mats_dyn,P_mats,phase)
      !---------- calculating partiton function (accumulation of real(phase))-----


      partition = partition + real((phase)*(pref/abs(pref)))
      !partition = partition + phase*(real(pref)/abs(real(pref)))

      !----------- C_xx(t=0)----------------------
      
     corr_xx0 = corr_xx0 &
        + real((phase)*(pref/abs(pref)))*Q_mats_dyn(1,0)*Q_mats_dyn(1,0)

      !call back_transform_matrix(Q_mats_dyn,q)

      !   corr_xx0 = corr_xx0 &
      !     + phase*(real(pref)/abs(real(pref)))*Q_mats_dyn(1,0)*Q_mats_dyn(1,0)
      !================dynamics step ==========================================
     
      do istep = 1,nstep

        call evolve(Q_mats_dyn,P_mats,x0,p0)


        corr_xx(istep) = corr_xx(istep) &
          + real((phase)*(pref/abs(pref)))*Q_mats(1,0)*Q_mats_dyn(1,0)

        !corr_xx(istep) = corr_xx(istep) &
        !  + phase*(real(pref)/abs(real(pref)))*Q_mats(1,0)*Q_mats_dyn(1,0)

      enddo 


    enddo 

    corr_xx = corr_xx/real(partition)
    corr_xx0 = corr_xx0/real(partition)

    write(3001,*) corr_xx0

    do ncount=1,nstep
      write(200,222) ncount*dt,corr_xx(ncount)
    end do

222      format(60(e13.6,2x))
    stop
end program mat_1d