!*********************************************************!
!*********************  ANT.G-2.5.2  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   David Jacob                                           !
!                                                         !
!      Max-Planck-Institute for Microstructure Physics    !
!      Halle, 06120 (GERMANY)                             !
!                                                         !
!*********************************************************!
  MODULE AntCommon
!*********************************************************!
! module for sharing a *few* common variables among       !
! several modules                                         !
! i.e a kind of common block for ANT                      !
!*********************************************************!
  implicit none

  ! Job name from title section of Gaussian input file
  ! Used for output file names
  CHARACTER(len=50) :: jobname
  ! basic file name for ant1d input files
  CHARACTER(len=50) :: ant1dname
  
  ! *** Unique file unit identifiers ***
  ! To avoid chaos every input/output file has a unique identifier

  integer, parameter :: ifu_log = 6

  integer, parameter :: ifu_nam = 10
  integer, parameter :: ifu_ini = 11

  integer, parameter :: ifu_xyz = 101
  integer, parameter :: ifu_tra = 102
  integer, parameter :: ifu_red = 103
  integer, parameter :: ifu_dos = 104
  integer, parameter :: ifu_ham = 105
  integer, parameter :: ifu_diagham = 110
  integer, parameter :: ifu_mul = 106
  integer, parameter :: ifu_msh = 107
  integer, parameter :: ifu_hyb = 108
  integer, parameter :: ifu_ac  = 109

  integer, parameter :: ifu_bl  = 20

  integer, parameter :: ifu_dm  = 30
  integer, parameter :: ifu_dmx = 31
  integer, parameter :: ifu_fm  = 32

  integer, parameter :: ifu_ant = 40
  integer, parameter :: ifu_lead = 41



  
end module AntCommon
