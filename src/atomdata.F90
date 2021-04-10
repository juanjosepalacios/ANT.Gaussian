!***************************************
!*                                     *
!*  ANT1D - atomdata.f90               *
!*                                     *
!*  Atomic data of each atom species   *
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
module AtomData
  use constants
  implicit none
  save

  ! Number of atoms
  !!integer :: NAtomData

  type TAtomData
     ! Conventional atomic number
     integer :: AN
     ! Number of Atom shells per atom
     integer :: NShells
     ! Number of Atomic orbitals per atom
     integer :: NAO
     ! Atomic Shell types for each atom (s,p,d,f,g,h...)
     character, dimension(MaxShAtom) :: AtShells
     ! Spin-orbit coupling for each Shell
     real(double), dimension(MaxShAtom) :: SOC
     ! Atomic orbitals for each atom (0,1,2,3,4,...)
     integer, dimension(MaxAOAtom) :: AOT
  end type TAtomData

  ! Array of atomic data
  type(TAtomData),dimension(MaxAtomData) :: AtomDataArr

contains

  !************************
  !  Reads in atomic data  
  !************************
  subroutine InitAtomData( iunit )
    use parameters, only: NAtomData
    use constants
    use messages
    implicit none

    integer, intent(in) :: iunit
    integer :: ios, iatom, ish, nshao, ishao, iao
   
    character(len=100) :: Keyword

    integer :: AN
    character, dimension(MaxShAtom) :: AtShells
    namelist/AtomData/AN,AtShells

    print *
    print *, "***************************************"
    print *, "* Reading atomic data from input file *" 
    print *, "***************************************"
    print *
    write (*,'(A,I2)'), " # of expected atomic data blocks: ", NAtomData

    do iatom=1,NAtomData
       print *
       write(*,'(A,I2)'), " Reading atomic data block #", iatom
       print *
       AN  = 0
       AtShells = '0'  
       read( unit=iunit, nml=AtomData )
       !!write(unit=*,nml=AtomData)
       AtomDataArr(iatom)%AN = AN
       AtomDataArr(iatom)%AtShells = AtShells       
       AtomDataArr(iatom)%NShells = 0
       AtomDataArr(iatom)%NAO = 0
       do ish=1,MaxShAtom
          if( AtShells(ish) == '0' )exit
          AtomDataArr(iatom)%Nshells = AtomDataArr(iatom)%Nshells + 1
          nshao = 0
          select case( AtomDataArr(iatom)%AtShells(ish) )
          case("s"); nshao = 1;  iao=0
          case("p"); nshao = 3;  iao=1
          case("d"); nshao = 5;  iao=4
          case("f"); nshao = 7;  iao=9
          case("g"); nshao = 9;  iao=15
          case("h"); nshao = 11; iao=25
          end select
          do ishao=1,nshao
             AtomDataArr(iatom)%NAO = AtomDataArr(iatom)%NAO + 1
             AtomDataArr(iatom)%AOT(AtomDataArr(iatom)%NAO) = iao; iao=iao+1
          end do
       end do
       write (*,'(A,I3)'), " AN ", AtomDataArr(iatom)%AN
       write (*,'(A,I3)'), " #Atomic Shells = ", AtomDataArr(iatom)%NShells       
       write (*,'(A,I3)'), " #Atomic Orbitals = ", AtomDataArr(iatom)%NAO
       write (*,'(A,100(I2))'), " Atomic Shells: ", AtomDataArr(iatom)%AOT(1:AtomDataArr(iatom)%NAO)
    end do
  end subroutine InitAtomData

  !***************************************************************
  ! Finds index of atomic data in data array for an atomic number
  !***************************************************************
  integer function FindAtomData( AN )
    use parameters
    implicit none

    integer, intent(in) :: AN
    integer :: idata

    do idata=1,NAtomData
       FindAtomData = idata
       if( AtomDataArr(idata)%AN ==  AN ) exit
    end do
    if( idata > NAtomData )FindAtomData = -1
    
  end function FindAtomData

end module AtomData

