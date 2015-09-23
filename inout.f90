!> Module containing in/output-related routines
module inout
  use accuracy
  use omp_lib
  implicit none
  
contains
  !> Reads the Laplacian and other Input
  !!
  !! \details This routine reads all the important input from the 
  !! laplacian.inp - file.
  !!
  !! \param steps  Timesteps for the simulation.
  !! \param dim  The system's dimension.
  !! \param extin  Particle Input per timestep.
  !! \param extout  Particle Output per timestep.
  !! \param laplacian  Laplacian matrix of the system.
  !! \param pp  Contains initial probability vector on output.
  !! \param dt  Time-discretization (length of one timestep).
  !! \param it  Tells the code how many times to run.
  !!
  !! \note The laplacian.example.inp - file contains an example which shows
  !! how the input-file has to look.
  subroutine read_laplacian(steps, dim, extin, extout, laplacian, pp, dt, it)
    integer, intent(out) :: steps, dim, it
    real(dp), intent(out) :: extin, extout
    integer, allocatable, intent(out) :: laplacian(:,:)
    real(dp), allocatable, intent(out) :: pp(:,:)
    real(dp), intent(out) :: dt
    integer :: ii
    
    open(10, file="laplacian.inp", status="old",&
         &form="formatted", action="read")
    read(10,*) steps
    read(10,*) it
    read(10,*) dim
    allocate(laplacian(dim,dim))
    do ii = 1, dim
       read(10,*) laplacian(ii,:)
    end do
    allocate(pp(dim,steps+1))
    pp = 0
    read(10,*) dt
    read(10,*) pp(:,1)
    read(10,*) extin, extout

    close(10, status="keep")    
  end subroutine read_laplacian

  
  !> Write currents into file
  !!
  !! \details Calculates the currents from the populations and saves
  !! them into a file.
  !!
  !! \param current  Previously allocated array which will contain the
  !! calculated currents.
  !! \param steps  Number of timesteps which have been used.
  !! \param jj  Integer which specifies which iteration of the code is
  !! running.
  !! \param tper  Simulation time of one iteration.
  !! \param dt  Discretization timestep.
  subroutine write_current(current,steps,jj,tper,dt)
    real(dp), intent(in) :: current(:,:), tper, dt
    integer, intent(in) :: steps, jj
    integer :: ii
    !Open File
    if(jj == 1) then
       open(11, file="current.dat", status="replace",&
            &form="formatted", action="write")
       else if(jj > 1) then
          open(11, file="current.dat", status="old",&
               &form="formatted", action="write",&
               &position="append")
       else
          write(*,"(A)") "Bad iteration number!"
          stop
       end if
       !Write currents
       do ii = 1, steps
          write(11,*) (jj-1) * tper + (ii - 1) * dt, current(:,ii)
       end do
       close(11, status="keep")
     end subroutine write_current
     

     
     !> Calculate and write voltage.
     !!
     !! \details  Calculates the voltages from the previously calculated
     !! populations and saves them into a file.
     !!
     !! \param voltage  Array which will contain the voltage.
     !! \param steps  Timesteps of the simulation.
     !! \param pp  Array which contains the populations.
     !! \param dim  The system's dimension.
     !! \param jj  Number which specifies which iteration of the code is
     !! running.
     !! \param tper  Simulation time of one iteration.
     !! \param dt  Time-discretization (length of one timestep).
     subroutine write_voltage(voltage,steps,pp,dim,jj,tper,dt)
       real(dp), allocatable, intent(out) :: voltage(:)
       real(dp), intent(in) :: pp(:,:), tper, dt
       integer, intent(in) :: steps, dim, jj
       integer :: ii
       
       allocate(voltage(steps))
       !Open File
       if(jj == 1) then
          open(11, file="volt.dat", status="replace",&
               &form="formatted", action="write")
       else if(jj > 1) then
          open(11, file="volt.dat", status="old",&
               &form="formatted", action="write",&
               &position="append")
       else
          write(*,"(A)") "Bad iteration number!"
          stop
       end if
       !Write voltage
       do ii = 1, steps
          voltage(ii) = pp(1,ii) - pp(dim,ii)
          write(11,*) (jj-1) * tper + (ii - 1) * dt, voltage(ii)
       end do
       close(11, status="keep")
       deallocate(voltage)
     end subroutine write_voltage
     
     !> Writes the populations/probabilities into a file.
     !!
     !! \details The previously calculated probabilities/populations
     !! are written in 2 file. One file contains them as rows, the second
     !! one as columns.
     !!
     !! \param dim  The system's dimension.
     !! \param pp  The previously calculated populations.
     !! \param pp2  Array which will contain the populations as columns
     !! on exit.
     !! \param steps  Timesteps of the simulation.
     !! \param jj  Integer which specifies which iteration is running.
     !! \param tper  Simulation time of one iteration.
     !! \param dt  Time-discretization (length of a timestep).
     subroutine write_prob(dim,pp,pp2,steps,jj,tper,dt)
       integer, intent(in) :: dim, steps, jj
       real(dp), intent(in) :: pp(:,:), tper, dt
       real(dp), allocatable, intent(out) :: pp2(:,:)
       integer :: ii
       !Rows
       !Open Files
       if(jj == 1) then
          open(11, file="histo.dat", status="replace",&
               &form="formatted", action="write")
          open(12, file="histo2.dat", status="replace",&
               &form="formatted", action="write")
       else if(jj > 1) then
          open(11, file="histo.dat", status="old",&
               &form="formatted", action="write",&
               &position="append")
          open(12, file="histo2.dat", status="old",&
               &form="formatted", action="write",&
               &position="append")
       else
          write(*,"(A)") "Bad iteration number!"
          stop
       end if
       !$OMP PARALLEL SECTIONS NUM_THREADS(2)
       !$OMP SECTION
       !Write as Rows
       do ii = 1, dim
          write(11,*) ii, pp(ii,:)
       end do
       close(11, status="keep")
       !$OMP SECTION
       !Write as Columns
       allocate(pp2(steps+1,dim))
       pp2 = transpose(pp)
       do ii = 1, steps+1
          write(12,*) (jj-1) * tper + (ii - 1) * dt, pp2(ii,:)
       end do
       close(12, status="keep")
       deallocate(pp2)
       !$OMP END PARALLEL SECTIONS
     end subroutine write_prob
     
   end module inout
