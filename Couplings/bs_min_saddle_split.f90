program transform 
! Description: This code splits the min_saddle file from the barrier search run into subfolders for each TLS.  Within each subfolder, 'min1', 'min2', and 'saddle' files are created, corresponding to the configuration of minimum 1, minimum 2, and the saddle point, respectively.  These results can be used to calculate parameters such as the relaxation time and coupling constant.
! Note 1: Must attach 3x3 lattice vector matrix to beginning of min_saddle for variable e1 to be read in
! Note 2: Input files must be called 'min_saddle' and 'pre_proc.dat'
  real, parameter :: pi =3.14159265
  character(len=40) :: input,output1,output,output2,output3,command,input1
  character(len=5) :: atom_name,dir
  real, dimension(3,3) :: e1
  real, dimension(3):: atom,center,min_1,min_2,saddle
  real ::p,d,sd1,sd2,a
  integer :: nmol, emp, natoms,mm=0
  integer :: int1=0, int2=3,i,j,nfiles,nindex
  integer :: n,flag
  real :: de,dd,de1,de2

  write(6,*) 'number of atoms:'
  read(5,*) natoms
  write(6,*) 'number of configurations:'
  read(5,*) nfiles

  input ='./min_saddle'
  open(10,file=input,status='old')
  input1='./pre_proc.dat'
  open(9,file=input1,status='old')

! Read in 3x3 lattice matrix from top of min_saddle file
  read(10,*)e1  

  do j=1,nfiles  
  read(9,*)n,flag,de,dd,de1,de2
  if(flag.eq.1)then  ! only use those configurations that pass all criteria from pre_processing code
    if(j<10)then
      write(dir,"(i1)")j
    elseif(j<100)then
      write(dir,"(i2)")j
    elseif(j<1000)then
      write(dir,"(i3)")j
    elseif(j<10000)then
      write(dir,"(i4)")j
    endif

! Create subdirectories for each 'good' TLS based on pre-preocessing screening (flag = 1)  
    command ='mkdir ' // trim(dir)
    CALL system(command) 
  
    output1='./'// trim(dir)//'/min1'
    output2='./'// trim(dir)//'/min2'
    output3='./'// trim(dir)//'/saddle'
    write(6,*)output1
    open(11,file=output1,status='replace')
    open(12,file=output2,status='replace')
    open(13,file=output3,status='replace')

    read(10,*) nindex

    write(11,*) 'CELL_PARAMETERS'
    write(11,'(3I10)')int1,int2,n
    write(11,'(3f20.14)') e1

    write(12,*) 'CELL_PARAMETERS'
    write(12,'(3I10)')int1,int2,n
    write(12,'(3f20.14)') e1

    write(13,*) 'CELL_PARAMETERS'
    write(13,'(3I10)')int1,int2,n
    write(13,'(3f20.14)') e1

    do i = 1,natoms
      read(10,*)atom_name, nindex
      read(10,*)min_1
      read(10,*)min_2
      read(10,*)saddle
      write(11,'(1A4,1I10)')atom_name, i
      write(11,'(3F20.15)')min_1
      write(12,'(1A4,1I10)')atom_name, i
      write(12,'(3F20.15)')min_2
      write(13,'(1A4,1I10)')atom_name, i
      write(13,'(3F20.15)')saddle
    enddo
      !print *, j
      !print *, atom_name, nindex
      !print *, min_1
    close(11)
    close(12)
    close(13)
  else
     read(10,*) nindex
    do i = 1,natoms
      read(10,*)atom_name
      read(10,*)min_1
      read(10,*)min_2
      read(10,*)saddle
    enddo

      !print *, j
      !print *, atom_name, nindex
      !print *, min_1

  endif 
  enddo
end program transform
