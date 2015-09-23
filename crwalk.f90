program crwalk
  use accuracy
  use inout
  use omp_lib

  implicit none

  !Declarations
  integer :: jj, kk, dim, steps, it
  integer, allocatable :: laplacian(:,:)
  real(dp), allocatable, target :: pp(:,:)
  real(dp), pointer :: tmpvec(:)
  real(dp), allocatable :: pp2(:,:), current(:,:), voltage(:)
  real(dp) :: tmp, dt, tper, extin, extout

  !Read Input
  call read_laplacian(steps, dim, extin, extout, laplacian, pp, dt, it)

  tmp = 0.0_dp
  !Time for one entire iteration of the code
  tper = dt * steps
  !Walk
  do jj = 1, it
     do kk = 1, steps
        tmpvec => pp(:,kk)
        pp(:,kk+1) = tmpvec - dt *  matmul(laplacian, tmpvec)
        !Stop-Condition(optional)
        !if(pp(dim,kk+1) > 1.0_dp /(2.0_dp * dim)) then
        !   write(*,"(A,F8.3)") "Probability reached at t = ", dt * kk
        !   stop
        !end if
        !Adding external particles
        pp(1,kk+1) = pp(1,kk+1) + extin
        if( pp(dim,kk+1) >= extout ) then
           pp(dim,kk+1) = pp(dim,kk+1) - extout
        end if
     end do
     
     !Calculate Currents(for serial circuits)
     !allocate(current(dim-1, steps))
     !do kk = 1, steps
     !   do ii = 1, dim-1
     !      current(ii,kk) = - laplacian(ii+1,ii) * pp(ii,kk)&
     !           & + laplacian(ii,ii+1) * pp(ii+1,kk)
     !   end do
     !end do
     
     !Write currents into file
     !call write_current(current,steps,jj,tper,dt)
     !Write Voltage into file
     call write_voltage(voltage,steps,pp,dim,jj,tper,dt)
     !Write probabilities into files
     call write_prob(dim,pp,pp2,steps,jj,tper,dt)
     !Reset
     pp(:,1) = pp(:,steps)
     !deallocate(current)
     !Write final probability vector and voltage on screen
     write(*,*) "Probability vector at the end:"
     write(*,*) pp(:,steps)
     write(*,"(A)",advance="no") "Voltage at the end: "
     write(*,"(F8.4)") pp(1,steps) - pp(dim,steps)
  end do
end program crwalk
