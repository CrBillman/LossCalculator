program transform 
! Code goes through all barriers listed in bisc.out file and categorizes them to determine which to keep for subsequent analyses
!     1) first, it codes all '1' and '5' barriers (those that found a good saddle point) and separates them from '2', '3', '4', and '0' codes.  All these latter codes are given a '-1' or '-5' and will no longer be used.
!     2) then, the code goes through all the good '1' and '5' codes and determines if they are duplicates by comparing their barrier heights and distances between minima using a threshold input below (a and b variables).  These are then given a '1' for a good double well pair (flag variable)
!     3) finally, these good double well pairs are filtered based on if their asymmetry is more than asym_thr - if they are, the flag is changed to '2'.  These barriers with their final codes are then output to pre_proc.dat, which is used as input for min_saddle_split.f90 to split each double well into the min1, min2, and saddle configurations for subsequent calculations.
  real, parameter :: pi =3.14159265
  character(len=30) :: input,output
  character(len=3) :: atom_name,dir
  real :: de_thr, dd_thr, asym_thr
  real, dimension(3,3) :: e1
  real, allocatable :: de(:), dd(:), de1(:),de2(:)
  integer, allocatable :: flag(:)
  integer :: nmol, emp, natoms,mm=0,n
  integer :: int1=0, int2=1,i,j,nfiles,best=0,best2=0

!-----FILTER CRITERIA------!
  de_thr = 0.0005  ! lower bound for asymmetry to compare duplicates
  dd_thr = 0.0005  ! lower bound for distance between minima to compare duplicates
  asym_thr = 1.0   ! upper bound for barrier asymmetry
!-----END FILTER CRITERIA---!

  write(6,*) 'number of configurations'
  read(5,*) nfiles
  input ='bisc.out'
  open(10,file=input,status='old')
  output ='pre_proc.dat'
  open(11,file=output,status='replace')

  natoms=1008
  allocate(de(nfiles),de1(nfiles),de2(nfiles),dd(nfiles),flag(nfiles))
! read bisection file
  do j=1,nfiles  
    read(10,*)mm,n,de(j),dd(j),flag(j),de1(j),de2(j) 
  enddo
  goodcount=0 
  do j=1,nfiles  
    if(flag(j).ne.5.and.flag(j).ne.1)then
! wrong flag bad structure skip
      flag(j)=-1
    elseif(de(j).le.0.0001.or.dd(j).le.0.0001)then
      flag(j)=-5
    else
      goodcount=goodcount+1
! compare de and dd to previous structures
      if(goodcount.gt.1)then 
        do i=1,j-1
          a=abs(de(j)-de(i)) 
          b=abs(dd(j)-dd(i))
          if(a.gt.de_thr.or.b.gt.dd_thr)then
!   accepted state 
            flag(j)=1
          else
!   rejected duplicate state based on de_thr and dd_thr thresholds
            flag(j)=0
            exit
          endif
        enddo
      else
! first state accept
        flag(j)=1
      endif
! rejected state based on asymmetry greater than asym_thr
      if(flag(j).eq.1.and.de(j).gt.asym_thr)then
        flag(j)=2
      endif
    endif
      write(11,"(2i6,4f20.15)")j,flag(j),de(j),dd(j),de1(j),de2(j)
      if(flag(j).eq.1)then
        best=best+1
      endif
      if(flag(j).eq.2)then
        best2=best2+1
      endif
  enddo
  write(6,*)goodcount,best,best2
end program transform
