!***************************************
!*                                     *
!*  ANT1D - geom.f90                   *
!*                                     *
!*  Atomic structure                   *
!*                                     *
!***************************************
!*                                     *
!*  This source file is part of the    *
!*  ANT1D project.                     *
!*                                     *
!*  Copyright (c) 2006 - 2015 by       *
!*                                     *
!*  David Jacob                        *
!*                                     *
!*  MPI fuer Mikrostrukturphysik       *
!*  Weinberg 2                         *
!*  06120 Halle                        *
!*  Germany                            *
!*                                     *
!***************************************
module geom
  use constants
  use filemaster
  use parameters


  implicit none

  !
  ! Data structure holding
  ! information of an atom:
  ! * index to atomic data
  ! * position of the atom
  ! * first and last atomic orbital
  !
  type TAtom
     integer :: index
     real(double) :: x, y, z
     integer :: firstao, lastao
  end type TAtom

  integer :: NDA
  type(TAtom), dimension(MaxDAtoms) :: DevAtoms
  integer :: NDevAOrbs
  
  integer, dimension(2) :: NLA
  type(TAtom), dimension(2,MaxDAtoms) :: LeadAtoms

contains

  subroutine InitGeom
    use parameters
    use messages
    implicit none
    
    integer :: iunit, ios, ia

!sf only read in geometry data of the lead if absorbing boundary conditions are disabled
  if (.not. abc) then
    ! Read geometry for electrode 1
    print *, "Reading xyz file for electrode No. 1:", Lead1XYZ
    iunit =fopen( Lead1XYZ, 'old', ios )
    if( ios /= 0 )call FileErr( "Geom/InitGeom", Lead1XYZ )
    NLA(1) = ReadXYZFile( iunit, LeadAtoms(1,:) )
    call fclose( iunit )

    ! Read geometry for electrode 2
    print *, "Reading xyz file for electrode No. 1:", Lead2XYZ
    iunit =fopen( Lead2XYZ, 'old', ios )
    if( ios /= 0 )call FileErr( "Geom/InitGeom", Lead2XYZ )
    NLA(2) = ReadXYZFile( iunit, LeadAtoms(2,:) )
    call fclose( iunit )
  end if

    ! Read device geometry
    print *, "Reading xyz file for device:", DevXYZ
    iunit =fopen( DevXYZ, 'old', ios )
    if( ios /= 0 )call FileErr( "Geom/InitGeom", DevXYZ )
    NDA = ReadXYZFile( iunit, DevAtoms )
    call fclose( iunit )
    NDevAOrbs=0
    do ia=1,NDA
       NDevAOrbs=NDevAOrbs+(DevAtoms(ia)%lastao-DevAtoms(ia)%firstao+1)
    end do

  end subroutine InitGeom

  !
  ! Reads xyz-file containing Number of atoms, atomic numbers and their positions
  !  returns number of atoms
  !
  ! Format:
  !
  ! <NAtoms>
  ! <AN1>  <x1> <y1> <z1>
  ! <AN2>  <x2> <y2> <z2>
  ! <AN3>  <x3> <y3> <z3>
  ! ...
  !
  integer function ReadXYZfile( iunit, AtomStruc )
    use AtomData
    use messages
    use util
    implicit none

    integer, intent(in) :: iunit
    type(TAtom), dimension(:),intent(out) :: AtomStruc

    integer :: NAtoms, AN, ia, firstao, ios
    character(len=100) :: ANStr
    
    read(unit=iunit,fmt=*,iostat=ios), NAtoms
    if( ios /= 0 ) call ReadErr( "Geom/ReadXYZFile", "NDAtoms" )
    print *, "Atoms read from xyz-file:"
    print *, "      Atom#          AN       index      1st AO     Last AO"

    firstao = 1

    do ia=1,NAtoms
       read(unit=iunit,fmt=*,iostat=ios), AN, AtomStruc(ia)%x, AtomStruc(ia)%y, AtomStruc(ia)%z
       if( ios /= 0 ) call ReadErr( "Geom/ReadXYZFile", "AN, AtomStruc(ia)%x, AtomStruc(ia)%y, AtomStruc(ia)%z" )
       !
       ! Find index to Atomic Data array 
       ! 
       AtomStruc(ia)%index = FindAtomData( AN )
       if( AtomStruc(ia)%index <= 0 )then
          call int2str( AN, ANStr )
          call ErrMessage( "Geom/ReadParameters", "No atomic data found for AN "//trim(ANStr)//".", .true. )
       end if
       !
       ! First and last atomic orbital on Atom
       !
       AtomStruc(ia)%firstao=firstao
       firstao = firstao + AtomDataArr(AtomStruc(ia)%index)%NAO
       AtomStruc(ia)%lastao= firstao - 1
       print *, ia, AN, AtomStruc(ia)%index, AtomStruc(ia)%firstao, AtomStruc(ia)%lastao
    end do

    ReadXYZFile = NAtoms
    
    !!print *, "Number of device atoms = ", NDAtoms
    !!print *, "Atom indices = ", AtomStruc(1:NDAtoms)%index

  end function ReadXYZfile

end module geom
