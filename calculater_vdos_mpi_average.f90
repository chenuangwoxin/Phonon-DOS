!yangyu's workstation
!mpiifort -o vdos_mpi calculater_vdos_mpi.f90 /home/yangyu/Data/zhongwei/fftw/lib/libfftw3.a
!Todai's workstation
!mpiifort -o calculater_vdos_mpi calculater_vdos_mpi.f90 /work/po0043/o00043/software/fftw/lib/libfftw3.a
!-----to calculate the vibrational density of states from a parallel running code


program main
use mpi
include '/public1/home/jingchao/fftw/1/include/fftw3.f'

integer i,j,k,m,n
integer i_sample,n_sample,i_step,atom_id
integer Ln,num,num_cell,i_v,n_omega,n_calculation                      !settings
integer n_interval,nbins,nlest,recvcount
double precision,allocatable:: atom_number_node(:),atom_number(:),q_mode(:,:),averaged_p(:,:)
!double precision,allocatable:: modal_real(:,:),modal_imag(:,:),modal_ve_real(:,:),modal_ve_imag(:,:)
double precision,allocatable:: ve_x(:),ve_y(:),ve_z(:),vdos_omega(:,:),vdos_output(:,:),vdos_all(:,:),mass(:)
double precision:: para_units,ratio
double complex,allocatable:: vdos0(:),id_velocity(:),md_velocity(:)
double complex::  temp1,temp2,q
real, parameter :: pi=3.1415926,m2kg=1.66054e-27,Atom=1.0e-10,pstos=1.0e-12,ev2J=1.60218e-19
character tmp 
integer	status(MPI_STATUS_SIZE)
integer(kind=MPI_OFFSET_KIND)	offset
integer  fh
integer,allocatable:: sendcount(:),displs(:)
character*16 infile1,infile2,infile3
character*7 read_num
character(len=10) :: file_id
character(len=50) :: file_name

integer ierr,my_id,num_procs,root_proc,filetype

integer :: plan

!----------------------------------------------------------------
 
integer :: tx
character( Len = 85 ) :: cStr

call CPU_TIME(time0)
    
root_proc=0
    
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)
 
!-----------------------all nodes setting------------------------

open(1,file='parameters.dat')
read(1,*)
read(1,*) num,dt,Ln
read(1,*)
read(1,*) n_sample,ratio
close(1)

para_units=28.0855*(m2kg/ev2J)*(Atom/pstos)**2
para_units=para_units/2.0 *dt           !coefficient in 1/2mv^2 and fdt in fourier

allocate(mass(num))

mass=28.0855
!open(1,file='water.EMD')
!do i=1,15
!read(1,*)
!enddo
!do i=1,num
!    read(1,*)j,k,m
!    if(m==1) then
!    mass(j)=15.9994
!    else
!    mass(j)=1.00794
 !   endif
!enddo
!----------------------parallel message sending------------------

allocate(sendcount(num_procs))
allocate(displs(num_procs)) 

nbins=num/num_procs
nlest=mod(num,num_procs)

do i=1,num_procs	

	if(i>nlest) then		
		sendcount(i)=nbins
	else
		sendcount(i)=nbins+1
	endif

enddo
 
displs(1)=0
do i=2,num_procs	 
	 displs(i)=displs(i-1)+sendcount(i-1)
	 displs(i)=displs(i)
enddo

do i=1,num_procs
	
	if(my_id==(i-1)) then
		allocate(atom_number_node(sendcount(i)))
		atom_number_node(:)=0.0
		recvcount=sendcount(i)
	endif
enddo


!-----------------------------------------------------------------

if (my_id==root_proc) then
	
	allocate(atom_number(num))

	do i=1,num
		atom_number(i)=i
	enddo

	write(*,*) '------------------------------------------------------------'
	write(*,*) '-----------Atomic-Vibrational Density of States-------------'
	write(*,*) '---------------------By Zhongwei Zhang----------------------'
	write(*,*) '---------------Email: zhongwei@tongji.edu.cn----------------'
	write(*,*) '----------------------Date: Feb. 1 2021---------------------'
	write(*,*) '------------------------------------------------------------'
	write(*,*) ' '
	write(*,*) "-----------Successfully prepare the data!"
	write(*,*) "------------------------------------------------------------"
	write(*,*) ' '

	write(*,*) "--------------------------------------------------" 
	write(*,*) 'Number of atoms:',num
    	write(*,*) 'Number of steps:',Ln
    	write(*,*) 'Timesteps:',dt
	write(*,*) "*******************"
	write(*,*) 'Number of step:',int(Ln*ratio/2)
	write(*,*) "*******************"
	write(*,*) "--------------------------------------------------"

	!-------------------------------------------------------------------------------------------------
   
	write(*,*) "------------------The calculation is begin:"

endif


!------------------------------------------------------------------
if (my_id==root_proc) then
write(*,*) "-----------------Atomic_velocity reading is begin:"
endif

allocate(ve_x(num*Ln))
allocate(ve_y(num*Ln))
allocate(ve_z(num*Ln))


call MPI_File_open(MPI_COMM_WORLD,'vx.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
call MPI_File_read_all(fh,ve_x,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
call MPI_File_close(fh, ierror)

call MPI_File_open(MPI_COMM_WORLD,'vy.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
call MPI_File_read_all(fh,ve_y,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
call MPI_File_close(fh, ierror)

call MPI_File_open(MPI_COMM_WORLD,'vz.bin',MPI_MODE_RDONLY,MPI_INFO_NULL,fh,ierr)
call MPI_File_read_all(fh,ve_z,1*num*Ln, MPI_DOUBLE_PRECISION, status,ierror)
call MPI_File_close(fh, ierror)

!------------------------------------------------------------------

if (my_id==root_proc) then
	write(*,*) '------------------------------------------------------------'
	write(*,*) 'Vdos mpirun calculating..'
	write(*,*) '------------------------------------------------------------'
endif

! Distribute the matrix A into different processes and store in the local matrix A1 

call MPI_SCATTERV(atom_number(:),sendcount,displs,MPI_DOUBLE_PRECISION,atom_number_node(:),recvcount,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)


!--------------------------calculation module---------------------------------------

allocate(id_velocity(Ln))
id_velocity=cmplx(0.0,0.0)
n_interval=int(Ln*ratio)
n_omega=int(n_interval/2)
allocate(vdos_omega(3,n_omega))
vdos_omega=0.0
allocate(md_velocity(n_interval))
allocate(vdos0(n_interval))
md_velocity=cmplx(0.0,0.0)
vdos0=cmplx(0.0,0.0)


do i=1,sendcount(my_id+1)
 
        if(my_id==root_proc) then
            write(*,'(A21,2I7)') ' In calculating step',i,sendcount(my_id+1)
        endif

!--------------------atom numer and velocity--------------------------
	
		atom_id=int(atom_number_node(i))

!--------------------fftw to vdos--------------------------
 
		do i_step=1, Ln
			id_velocity(i_step)=cmplx(ve_x(atom_id+(i_step-1)*num),0.0)
		enddo	 

        do i_sample=1,n_sample

            n1=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+1
            n2=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+n_interval
         
            md_velocity=id_velocity(n1:n2)
            md_velocity(:)=md_velocity(:)-sum(md_velocity(:))/float(n_interval)
            
            call dfftw_plan_dft_1d(plan,n_interval,md_velocity(:),vdos0,FFTW_FORWARD,FFTW_ESTIMATE)
            call dfftw_execute_dft(plan,md_velocity(:),vdos0)
            call dfftw_destroy_plan(plan)

            vdos_omega(1,:)=vdos_omega(1,:)+mass(atom_id)*(real(vdos0(:n_omega))**2+imag(vdos0(:n_omega))**2)/float(n_interval)
	     
        enddo 
        
		do i_step=1, Ln
			id_velocity(i_step)=cmplx(ve_y(atom_id+(i_step-1)*num),0.0)
		enddo	

        do i_sample=1,n_sample

            n1=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+1
            n2=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+n_interval

            md_velocity=id_velocity(n1:n2)
            md_velocity(:)=md_velocity(:)-sum(md_velocity(:))/float(n_interval)
            
            call dfftw_plan_dft_1d(plan,n_interval,md_velocity(:),vdos0,FFTW_FORWARD,FFTW_ESTIMATE)
            call dfftw_execute_dft(plan,md_velocity(:),vdos0)
            call dfftw_destroy_plan(plan)

            vdos_omega(2,:)=vdos_omega(2,:)+mass(atom_id)*(real(vdos0(:n_omega))**2+imag(vdos0(:n_omega))**2)/float(n_interval)

        enddo
 
		do i_step=1, Ln 
			id_velocity(i_step)=cmplx(ve_z(atom_id+(i_step-1)*num),0.0)
		enddo	
        id_velocity(:)=id_velocity(:)-sum(id_velocity(:))/float(Ln)
 
        do i_sample=1,n_sample

            n1=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+1
            n2=(i_sample-1)*int((Ln-n_interval)/(n_sample-1))+n_interval

            md_velocity=id_velocity(n1:n2)
            md_velocity(:)=md_velocity(:)-sum(md_velocity(:))/float(n_interval)
            
            call dfftw_plan_dft_1d(plan,n_interval,md_velocity(:),vdos0,FFTW_FORWARD,FFTW_ESTIMATE)
            call dfftw_execute_dft(plan,md_velocity(:),vdos0)
            call dfftw_destroy_plan(plan)

            vdos_omega(3,:)=vdos_omega(3,:)+mass(atom_id)*(real(vdos0(:n_omega))**2+imag(vdos0(:n_omega))**2)/float(n_interval)

        enddo
	 
enddo 

vdos_omega=vdos_omega/float(n_sample)
 
allocate(vdos_all(n_omega,3))

allocate(vdos_output(n_omega,num_procs))

call MPI_GATHER(vdos_omega(1,:),n_omega,MPI_DOUBLE_PRECISION,vdos_output,n_omega,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)

call MPI_Barrier(MPI_COMM_WORLD,IERROR)

forall(i=1:n_omega)vdos_all(i,1)=sum(vdos_output(i,:),dim=1)


call MPI_GATHER(vdos_omega(2,:),n_omega,MPI_DOUBLE_PRECISION,vdos_output,n_omega,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)

call MPI_Barrier(MPI_COMM_WORLD,IERROR)

forall(i=1:n_omega)vdos_all(i,2)=sum(vdos_output(i,:),dim=1) 

call MPI_GATHER(vdos_omega(3,:),n_omega,MPI_DOUBLE_PRECISION,vdos_output,n_omega,MPI_DOUBLE_PRECISION,&
                    root_proc,MPI_COMM_WORLD,ierr)

call MPI_Barrier(MPI_COMM_WORLD,IERROR)

forall(i=1:n_omega)vdos_all(i,3)=sum(vdos_output(i,:),dim=1)

if(my_id==root_proc) then

open(1,file='vdos_mpi.dat')
vdos_all=vdos_all/float(num)
do i=1,n_omega
	write(1,'(5e12.5)') 1.0/dt/float(n_interval)*float(i),vdos_all(i,1:3),sum(vdos_all(i,1:3))
enddo
close(1)
!------------------------------------------------------------------------
write(*,*)
write(*,*) "Vdos calculation is successfully done!"
write(*,*)
!------------------------------------------------------------------------
  
write(*,*)
write(*,*) "Work is well done!"

endif

call MPI_FINALIZE(ierr)

call CPU_TIME(time1)


if(my_id==root_proc) then

!------------------------------------------------------------------------
write(*,*)
write(*,*) "The total used CPU time is:",time1-time0,'second.' 
write(*,*)
!------------------------------------------------------------------------

endif

end program
