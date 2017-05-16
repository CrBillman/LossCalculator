program transform 
! DESCRIPTION: After running the automated wall relaxation and calculation of changes in energy with strain for each minimum in each TLS, this post-processing code calculates the slope of the asymmetry vs strain curve for each TLS.
! INPUT: 1) pre_proc.dat: file generated from pre_processing code to determine 'good' TLSs (non-duplicates below asymmetry threshold)
!        2) 1eng.dat: energy vs strain curve calculated for first minimum of TLS (created by cycle_postcouple_script)
!        3) 2eng.dat: energy vs strain curve calculated for second minimum of TLS
! OUTPUT: 1) result_cpl: asymmetry vs strain curve for each TLS and for each strain direction.
!                        Column 1: strain magnitude, Column 2-7: Energy for xx, yy, zz, xy, xz, yz strains, Column 8: TLS index
!         2) result_slope: Slope values using linear fitting for each TLS
!                        Column 1: TLS index, Columns 2-7: Slope for xx, yy, zz, xy, xz, yz strains
  logical f1,f2,f3,f
  real, parameter :: pi =3.14159265
  character(len=30) :: input,input1,output,output2,input2,input3,command
  character(len=4) :: atom_name,dir
  real :: counter,sumx,sumx2,sumxy,sumy,sumy2
  real, dimension(3):: atom,center,min_1,min_2,saddle
  real(8),allocatable :: e1(:,:),e2(:,:),delta(:,:),eps1(:,:),eps2(:,:),gama(:,:),slope(:)
  integer :: nmol, emp, natoms,mm=0
  integer :: int1=0, int2=1,i,j,nfiles
  integer :: n,flag

! Input number of atoms and number of TLS configurations'
  write(6,*) 'number of atoms'
  read(5,*) natoms
  write(6,*) 'number of configurations'
  read(5,*) nfiles

  input='./pre_proc.dat'
  output='result_cpl'
  output2='result_slope'
  open(10,file=input,status='old')
  open(16,file=output,status='replace')
  open(17,file=output2,status='replace')

  write(17,"(A60)") "N   XX   YY   ZZ   XY   XZ   YZ"

  allocate(e1(6,21),e2(6,21),eps1(6,21),eps2(6,21))
  allocate(delta(6,21),gama(6,21))
  allocate(slope(6))

! Calculate slope for each strain direction and for each TLS (nfiles)
  do j=1,nfiles  
  read(10,*)n,flag,de,dd,de1,de2
  write(6,*) n
  if(flag.eq.1)then
    mm=mm+1
    if(j<10)then
      write(dir,"(i1)")j
    elseif(j<100)then
      write(dir,"(i2)")j
    elseif(j<1000)then
      write(dir,"(i3)")j
    elseif(j<10000)then
      write(dir,"(i4)")j
    endif

    ! Read values of energy vs strain from '1eng.dat' and '2eng.dat' from each TLS subdirectory   
    input1 ='./' // trim(dir)//trim('/1eng.dat')
    input2 ='./' // trim(dir)//trim('/2eng.dat')
    INQUIRE(FILE =input1,EXIST=f1)
    INQUIRE(FILE =input2,EXIST=f2)
    f=f1.and.f2
    if(f)then 
      open(11,file=input1,status='old')
      open(12,file=input2,status='old')
        read(11,*)  ! skip first line of 1eng.dat and 2eng.dat because it is energy at zero strain
        read(12,*)
      ! Read in value of strain (eps1) and corresponding energy for min1 (e1) and min2 (e2)
      do ii=1,6
        do i = 1,21
          read(11,"(1F20.11,E20.12)")eps1(ii,i),e1(ii,i)
          read(12,"(1F20.11,E20.12)")eps2(ii,i),e2(ii,i)
          ! Calculate asymmetry (delta) at each value of strain
          delta(ii,i)=e1(ii,i)-e2(ii,i)
       enddo
      enddo
      ! Write out strain and energy values to 'result_cpl'
      do i=1,21
         write(16,"(f20.11,2x)",advance='no')eps1(1,i)
         do ii=1,6
            write(16,"(E20.12,4x)",advance='no')delta(ii,i)
         enddo
         ! Perform linear regression to find slope of asymmetry versus strain
         write(16,"(i6)")n
      enddo
      write(17,"(I6,2x)",advance='no')n
      do ii=1,6
         slope(ii)=0.0d0
         counter = 0.0d0
         sumx=0.0d0
         sumx2=0.0d0
         sumxy=0.0d0
         sumy=0.0d0
         sumy2=0.0d0
         do i=1,21
            counter = counter + 1.0d0
            sumx = sumx + eps1(ii,i)
            sumx2 = sumx2 + eps1(ii,i)*eps1(ii,i)
            sumxy = sumxy + eps1(ii,i)*delta(ii,i)
            sumy = sumy + delta(ii,i)
            sumy2 = sumy2 + delta(ii,i)*delta(ii,i)
         enddo
         slope(ii) = (counter*sumxy - sumx*sumy) / (counter*sumx2 - sumx**2)
         ! Write out slope values for each strain direction to result_slope
         write(17,"(E20.12,2x)",advance='no')abs(0.5*slope(ii))
      enddo
      write(17,*)
    else
      write(16,"(i6,4x,'one or more file dont exist')")n
    endif
  else
  endif 
  close(11)
  close(12)
  enddo
end program transform
