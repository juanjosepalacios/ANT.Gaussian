!*********************************************************!
!*********************  ANT.G-2.3.4  *********************!
!*********************************************************!
!                                                         !
!   Copyright (c) by                                      !
!                                                         !
!   Juan Jose Palacios                                    !
!   Maria Soriano                                         !
!                                                         !
!   Departamento de Fisica de la Materia Condensada       !
!   Universidad Autonoma de Madrid.                       !
!   Spain.                                                !
!                                                         !
!   Octubre 2012                                          !
!                                                         !
!*********************************************************!
  MODULE MolMod
!*********************************************************!
  use util
  use antcommon
#ifdef G03ROOT
    use g03Common, only: GetNAtoms, GetAtmChg
#endif
#ifdef G09ROOT
    use g09Common, only: GetNAtoms, GetAtmChg
#endif

  implicit none

contains

subroutine Mol_Sub(HD,SD,PD,shift)

use parameters, only: EW1, EW2
use parameters, only: NCorrBlM, CorrBegM, CorrEndM,IP,EA,ForIntCh
use numeric, only: RMatPow,RSDiag
use constants
use Cluster, only : hiaorbno, loaorbno

implicit none

real*8, dimension(:,:,:),intent(inout) :: HD 
real*8, intent(in) :: shift
real*8, dimension(:,:,:),intent(in) :: PD
real*8, dimension(:,:),intent(in) :: SD

integer :: NHyb,info,n,i,j,ispin,CBeg,CEnd,l,k,NSpin
real*8, allocatable :: HD_M(:,:),eigenval(:),XD_M(:,:),SPN_M(:,:),e_PD(:,:)
real*8, allocatable :: XD_MU(:,:),eu_PD(:,:),HD_Mu(:,:),SPP_M(:,:),ICh(:,:)
real*8, dimension(size(HD,2),size(HD,3)) :: HD_A
real*8 :: G_h,G_e,E_H,E_L,E_Lu,E_Hu
real*8 :: suma

NSpin = size(HD,1)

if (IP == EA) print*, "War: IP = EA", IP,EA

do l = 1,NCorrBlM !Molecular Correlated Blocks
  CBeg = loaorbno(CorrBegM(l))
  CEnd = hiaorbno(CorrEndM(l))

print *
print *, "----------------------------------------------"
print *, "--- Molecular Gap Correction    ---"
print *, "----------------------------------------------"
print *

nhyb = (CEnd-CBeg) + 1

print *, " Correlated Sub-Shell dimension = ", nhyb,"Beg_basis function = ",CBeg &
& ,"End_basis function = ", CEnd
print *, "----------------------------------------------"

allocate(HD_M(NHyb,NHyb),SPP_M(NHyb,NHyb),ICh(Nspin,NHyb))
allocate(XD_M(nhyb,nhyb),SPN_M(NHyb,NHyb))
allocate(e_PD(NHyb,NHyb),XD_Mu(nhyb,nhyb),eu_PD(NHyb,NHyb),HD_Mu(NHyb,NHyb))

  print*, "Orthogonalization of correlated sub-shell"
  print *, "----------------------------------------------"
  call RMatPow( SD(CBeg:CEnd,CBeg:CEnd), -0.5d0, SPN_M )
  call RMatPow( SD(CBeg:CEnd,CBeg:CEnd), 0.5d0, SPP_M )

do ispin=1,NSpin !spin loop

print*, "Spin =", ispin

!  HD_A = HD(ispin,:,:)

  e_PD = matmul(PD(ispin,CBeg:CEnd,CBeg:CEnd),SD(CBeg:CEnd,CBeg:CEnd)) !Charge Matrix

  HD_M = HD_A(cBeg:CEnd,cBeg:CEnd) 

  HD_M = matmul( SPN_M, HD_M )
  HD_M = matmul( HD_M, SPN_M ) 

!---------------------------------------------------------------------------
! Obtain the molecular Coefficient 
! We have print the coefficient to study the molecular orbitals.
!-------------------------------------------------------------------------

  allocate( eigenval(nhyb))

  info = 0
  call RSDiag( HD_M, eigenval, info ) 

  if (info /= 0 ) print*, "Warning, INFO /= ", info
  
  ! Charge Molecular Matrix
  e_PD = matmul(transpose(HD_M),matmul(e_PD,HD_M)) 

   XD_M = HD_M

  do i=1,nhyb
    do j=1,nhyb
      XD_M(i,j) = XD_M(i,j)*XD_M(i,j)
    end do
  end do

  print *, "--------------------------------------------------------------------"
  print *, "--------- Print Norm. Atomic Composition near the Fermi Level   ----"
  do i = 1, NHyb
    if (eigenval(i)+shift <= 0.0) ICh(ispin,i)=d_one
    if (eigenval(i)+shift > 0.0) ICh(ispin,i)=d_zero
    if (i < NHyb) then
      if (eigenval(i)+shift <= 0.0 .and. eigenval(i+1)+shift > 0.0) then
        E_H = eigenval(i)
        E_L = eigenval(i+1)
        if (ispin > 1) then
          if (E_H < E_Hu) E_H = E_Hu
          if (E_L > E_Lu) E_L = E_Lu
        end if
        print *, "--------------------- HOMO and LUMO Located ---------------------"
        write(6,'(A20,F20.10,A20,F20.10)'), "E(HOMO) - \mu = ", E_H+shift, "E(LUMO) - \mu= ", E_L+shift
      end if
    end if
    if (eigenval(i)+shift > EW1 .and. eigenval(i)+shift < EW2) then
      print *, "----------------------------------------------"
      write(6, '(A20,I6,A20,F20.10,A20,F20.10,A20,F20.10)') "Mol. Orb. = ",i,"Eigenvalue - \mu = ",eigenval(i)+shift, "ne = ", &
       &(2.0d0/NSpin)*e_PD(i,i), "nh = ", (2.0d0/NSpin)*(1.0d0-e_PD(i,i))
      print *, "----------------------------------------------"
      do j=1,nhyb
        write(6, '(I6,10F20.10)') j, XD_M(j,i)*100
      end do
    end if
  end do

  XD_M = HD_M 
  
  HD_M = HD_A(cBeg:CEnd,cBeg:CEnd) 

  HD_M = matmul( SPN_M, HD_M )
  HD_M = matmul( HD_M, SPN_M )  

  HD_M = matmul(transpose(XD_M),matmul(HD_M,XD_M)) 

  suma = 0.0d0
  do i=1,Nhyb
   suma = suma + ((2.0d0/NSpin)*e_Pd(i,i))
  end do
  print*, "Total electron Charge in Molecule = ", suma

  suma = 0.0d0
  do i=1,Nhyb
   suma = suma + ((2.0d0/NSpin)*(1.0d0-e_Pd(i,i)))
  end do
  print*, "Total Hole Charge in Molecule =  = ", suma

  ! Print Molecular Orbital
  if (ispin == 1)  open (7,file=trim(jobname)//'.MOA.dat')
  if (ispin == 2)  open (7,file=trim(jobname)//'.MOB.dat')
  Do J=0,nhyb
    if (J == 0) then
      write (7,*) (eigenval(i), i=1,nhyb)
    else
      call flush(7)
      write (7,*) (XD_M(j,i), i=1,nhyb)
    end if
  END DO
  deallocate(eigenval)
  close(7)

  if (ispin == 1) then
    HD_Mu = HD_M
    XD_Mu = XD_M
    eu_PD = e_PD
    E_Lu = E_L
    E_Hu = E_H
  end if
  
  If (ispin == 2 .or. NSpin == 1) then !Evaluation IF

    G_h = E_H - IP
    G_e = EA - E_L

    print *, "-------------Energy Parameters -------------------------------------------------------"
    write(6,'(A20,F20.10,A20,F20.10,A20,F20.10,A20,F20.10)') "E_H = ", E_H, "E_L = ", E_L, "IP = ", IP, "EA = ", EA
    write(6,'(A20,F20.10,A20,F20.10,A20,F20.10)') "G_h = ", G_h, "G_e = ", G_e, "Gap =", G_h+G_e 
    write(6,'(A20,F20.10)') "Chemical Potential = ", -shift
    print *, "--------------------------------------------------------------------------------------"

    do k = 1,Nspin

      HD_A(cBeg:CEnd,cBeg:CEnd)=d_zero
      if (ForIntCh) then
        do i=1,nhyb
          if (k == 1) HD_A(i+CBeg-1,i+CBeg-1) = HD_Mu(i,i) - G_h*ICh(k,i) + G_e*(d_one-ICh(k,i))
          if (k == 2) HD_A(i+CBeg-1,i+CBeg-1) = HD_M(i,i) - G_h*ICh(k,i) + G_e*(d_one-ICh(k,i))
        end do
      else
        do i=1,nhyb
!         print*, "HD B = ", i,HD_M(i,j)+shift
          if (k == 1) HD_A(i+CBeg-1,i+CBeg-1) = HD_Mu(i,i) - G_h*eu_PD(i,i) + G_e*(d_one-eu_PD(i,i))
          if (k == 2) HD_A(i+CBeg-1,i+CBeg-1) = HD_M(i,i) - G_h*e_PD(i,i) + G_e*(d_one-e_PD(i,i))
!         print*, "HD A = ", i,HD_A(i+CBeg-1,j+CBeg-1)+shift
        end do 
      end if

      if (k == 1) then
        HD_Mu = matmul(XD_Mu,matmul(HD_A(cBeg:CEnd,cBeg:CEnd),transpose(XD_Mu)))
        HD_Mu = matmul(SPP_M,matmul(HD_Mu,SPP_M))
        HD_A(cBeg:CEnd,cBeg:CEnd) = HD_Mu
      else
        HD_M = matmul(XD_M,matmul(HD_A(cBeg:CEnd,cBeg:CEnd),transpose(XD_M)))
        HD_M = matmul(SPP_M,matmul(HD_M,SPP_M))  
        HD_A(cBeg:CEnd,cBeg:CEnd) = HD_M
      end if

      HD(k,cBeg:CEnd,cBeg:CEnd) = HD_A(cBeg:CEnd,cBeg:CEnd)

    end do ! K loop
  end if ! Evaluation IF
end do !spin loop

print *, "----------------------------------------------"
print*, " End - Molecular transformation"
print *, "----------------------------------------------"

deallocate(HD_M,XD_M,SPN_M,e_PD,ICh)
deallocate(HD_Mu,XD_Mu,eu_PD,SPP_M)

end do !Molecular correlated Bolcks

end subroutine Mol_Sub

!!   Subroutine MO_input(A,S,shift)
!!   use constants
!!   use Cluster, only : hiaorbno, loaorbno
!!   use parameters, only: EW1, EW2, MONortho,MOOrtho,ELtype
!!   use parameters, only: SubSPBEG,SUBSPEND,AORB
!!   use numeric, only: RSDiag,RM_Ortho,R_Ortho,RNormCond
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:,:), intent(in) :: A
!!   real*8, dimension(:,:), intent(in) :: S
!!   real*8, intent(in) :: shift
!!   
!!   integer :: dimi, info, i,j, NSpin,k,l,is,CBeg,CEnd
!!   real*8, allocatable :: AT(:,:),ST(:,:),XT(:,:),eigenval(:)
!!   
!!   NSpin = size(A,1)
!!   
!!   !if (NCorrBlM < 1) NCorrBlM = 1
!!   
!!   do is=1,NSpin
!!   !  do l = 1,NCorrBlM !Molecular Correlated Blocks
!!   !    IF (ElType(1) == "GHOST" .and. ElType(2) == "GHOST") then
!!   if (AORB) then
!!         CBeg = subspbeg
!!         CEnd = subspend
!!       ELSE
!!   !      CBeg = loaorbno(CorrBegM(l))
!!   !      CEnd = hiaorbno(CorrEndM(l))
!!         CBeg = loaorbno(subspbeg)
!!         CEND = hiaorbno(subspend)
!!       END IF
!!   
!!       dimi=CEnd-CBeg+1
!!   
!!       allocate(AT(dimi,dimi),ST(dimi,dimi),XT(dimi,dimi),eigenval(dimi))
!!   
!!       AT = A(is,CBeg:CEnd,CBeg:CEnd)
!!       ST = S(CBeg:CEnd,CBeg:CEnd)
!!   
!!       call R_Ortho(0,AT,ST)
!!   
!!       if (MOOrtho) then
!!         print *, "----------------------------------------------"
!!         write(6, '(A100)') "Orthogonal Molecular Orbital"
!!         print *, "----------------------------------------------"
!!   
!!         call dH_MO(AT,XT,is,shift)
!!       end if
!!       if (MONortho) then
!!         print *, "----------------------------------------------"
!!         write(6, '(A100)') "No-Orthogonal Molecular Orbital"
!!         print *, "----------------------------------------------"
!!   
!!         AT = A(is,CBeg:CEnd,CBeg:CEnd)
!!         ST = S(CBeg:CEnd,CBeg:CEnd)
!!   
!!         call R_Ortho(0,AT,ST)
!!   
!!         XT = AT
!!   
!!         call RSDiag( XT, eigenval, info )
!!         if (info /= 0 ) print*, "Warning, INFO /= ", info,  "in MO_input, RSDiag subroutine"
!!   
!!         call RM_Ortho(0,XT,ST)
!!         call RNormCond(XT)
!!   
!!         AT = XT
!!   
!!         do i=1,dimi
!!           do j=1,dimi
!!             XT(i,j) = AT(i,j)*AT(i,j)
!!           end do
!!         end do
!!   
!!         do i=1,dimi
!!           If (eigenval(i)+shift > EW1*100 .and. eigenval(i)+shift < EW2*100) then
!!             print *, "----------------------------------------------"
!!             write(6, '(A20,I6,A20,2F20.10)') "Mol. Orb. = ",i,"Eigenvalue - \mu = ",eigenval(i)-shift,eigenval(i)
!!             print *, "----------------------------------------------"
!!             do j=1,dimi
!!               write(6, '(I6,10F20.10)') j, XT(j,i)*100
!!             end do
!!           end if
!!         end do
!!   
!!         XT = AT
!!         ! Print Molecular Orbital
!!         if (is == 1)  open (7,file=trim(jobname)//'.MONU.dat')
!!         if (is == 2)  open (7,file=trim(jobname)//'.MOND.dat')
!!         Do J=0,dimi
!!           if (J == 0) then
!!             write (7,*) (eigenval(i), i=1,dimi)
!!           else
!!             call flush(7)
!!             write (7,*) (XT(j,i), i=1,dimi)
!!           end if
!!         END DO
!!         close(7)
!!       
!!       end if !MO, NMO
!!       deallocate(AT,ST,XT,eigenval)
!!     !end do !Ncorr
!!   end do !Spin loop
!!   
!!   end subroutine MO_input
!!   Subroutine MO_input3(A,S,chempota,chempotd,shift,shiftemp)
!!   use constants
!!   use Cluster, only : hiaorbno, loaorbno
!!   use parameters, only: EW1, EW2, MONortho,MOOrtho,ELtype
!!   use parameters, only: SubSPBEG,SUBSPEND,AORB
!!   use numeric, only: RSDiag,RM_Ortho,R_Ortho,RNormCond
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:,:), intent(in) :: A
!!   real*8, dimension(:,:), intent(in) :: S
!!   real*8, intent(in) :: shift,chempota,chempotd,shiftemp
!!   
!!   integer :: dimi, info, i,j, NSpin,k,l,is,CBeg,CEnd
!!   real*8, allocatable :: AT(:,:),ST(:,:),XT(:,:),eigenval(:)
!!   real*8 :: realshift
!!   
!!   NSpin = size(A,1)
!!   
!!   !if (NCorrBlM < 1) NCorrBlM = 1
!!   
!!   do is=1,NSpin
!!   !  do l = 1,NCorrBlM !Molecular Correlated Blocks
!!   !    IF (ElType(1) == "GHOST" .and. ElType(2) == "GHOST") then
!!   if (AORB) then
!!         CBeg = subspbeg
!!         CEnd = subspend
!!       ELSE
!!   !      CBeg = loaorbno(CorrBegM(l))
!!   !      CEnd = hiaorbno(CorrEndM(l))
!!         CBeg = loaorbno(subspbeg)
!!         CEND = hiaorbno(subspend)
!!       END IF
!!   
!!       dimi=CEnd-CBeg+1
!!   
!!   if (is == 1) realshift = shift -(chempota - shiftemp)
!!   if (is == 2) realshift = shift -(chempotd - shiftemp)
!!   
!!       allocate(AT(dimi,dimi),ST(dimi,dimi),XT(dimi,dimi),eigenval(dimi))
!!   
!!       AT = A(is,CBeg:CEnd,CBeg:CEnd)
!!       ST = S(CBeg:CEnd,CBeg:CEnd)
!!   
!!       call R_Ortho(0,AT,ST)
!!   
!!       if (MOOrtho) then
!!         print *, "----------------------------------------------"
!!         write(6, '(A100)') "Orthogonal Molecular Orbital"
!!         print *, "----------------------------------------------"
!!   
!!         call dH_MO(AT,XT,is,realshift)
!!       end if
!!       if (MONortho) then
!!         print *, "----------------------------------------------"
!!         write(6, '(A100)') "No-Orthogonal Molecular Orbital"
!!         print *, "----------------------------------------------"
!!   
!!         AT = A(is,CBeg:CEnd,CBeg:CEnd)
!!         ST = S(CBeg:CEnd,CBeg:CEnd)
!!   
!!         call R_Ortho(0,AT,ST)
!!   
!!         XT = AT
!!   
!!         call RSDiag( XT, eigenval, info )
!!         if (info /= 0 ) print*, "Warning, INFO /= ", info,  "in MO_input, RSDiag subroutine"
!!   
!!         call RM_Ortho(0,XT,ST)
!!         call RNormCond(XT)
!!   
!!         AT = XT
!!   
!!         do i=1,dimi
!!           do j=1,dimi
!!             XT(i,j) = AT(i,j)*AT(i,j)
!!           end do
!!         end do
!!   
!!         do i=1,dimi
!!           If (eigenval(i)+shift > EW1*10 .and. eigenval(i)+shift < EW2*10) then
!!             print *, "----------------------------------------------"
!!             write(6, '(A20,I6,A20,2F20.10)') "Mol. Orb. = ",i,"Eigenvalue - \mu =",eigenval(i)+realshift,eigenval(i)
!!             print *, "----------------------------------------------"
!!             do j=1,dimi
!!               write(6, '(I6,10F20.10)') j, XT(j,i)*100
!!             end do
!!           end if
!!         end do
!!   
!!         XT = AT
!!         ! Print Molecular Orbital
!!         if (is == 1)  open (7,file=trim(jobname)//'.MONU.dat')
!!         if (is == 2)  open (7,file=trim(jobname)//'.MOND.dat')
!!         Do J=0,dimi
!!           if (J == 0) then
!!             write (7,*) (eigenval(i), i=1,dimi)
!!           else
!!             call flush(7)
!!             write (7,*) (XT(j,i), i=1,dimi)
!!           end if
!!         END DO
!!         close(7)
!!   
!!       end if !MO, NMO
!!       deallocate(AT,ST,XT,eigenval)
!!     !end do !Ncorr
!!   end do !Spin loop
!!   
!!   end subroutine MO_input3
!!   
!!   Subroutine MO_input2(A,S,shift,chempot)
!!   use constants
!!   use Cluster, only : hiaorbno, loaorbno
!!   use parameters, only: EW1, EW2, MONortho,MOOrtho,ELtype
!!   use parameters, only: NCorrBlM, CorrBegM, CorrEndM,NCorrBl, CorrBeg, CorrEnd
!!   use numeric, only: RSDiag,RM_Ortho,R_Ortho,RNormCond,DVDiag,RMAtPow,CND2
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:,:), intent(in) :: A
!!   real*8, dimension(:,:), intent(in) :: S
!!   real*8, intent(in) :: shift,chempot
!!   
!!   real*8, dimension(size(A,2),size(A,3)) :: HS,InvSD,VSh,VSt
!!   integer :: dimi, info, i,j, NSpin,k,l,is,CBeg,CEnd
!!   real*8, allocatable :: AT(:,:),ST(:,:),XT(:,:),eigenval(:),w(:),wi(:),VR(:,:),VL(:,:),X(:,:)
!!   
!!   NSpin = size(A,1)
!!   
!!   if (NCorrBlM < 1) NCorrBlM = 1
!!   
!!   do is=1,NSpin
!!     do l = 1,NCorrBlM !Molecular Correlated Blocks
!!       IF (ElType(1) == "GHOST" .and. ElType(2) == "GHOST") then
!!         CBeg = CorrBeg(l)
!!         CEnd = CorrEnd(l)
!!       ELSE
!!         CBeg = loaorbno(CorrBegM(l))
!!         CEnd = hiaorbno(CorrEndM(l))
!!       END IF
!!   
!!       dimi=CEnd-CBeg+1
!!   
!!   InvSD = 0.0d0
!!   call RMatPow(S,-1.0d0,InvSD)
!!   
!!   VSh = 0.0d0
!!   VSh= shift*S
!!   Vst = 0.0d0
!!   Vst(CBeg:CEnd,CBeg:CEnd) = (chempot-shift)*S(CBeg:CEnd,CBeg:CEnd)
!!   
!!   HS = A(is,:,:)+VSh!+Vst
!!   HS = matmul(HS(:,:),InvSD)
!!   
!!       allocate(AT(dimi,dimi),ST(dimi,dimi),XT(dimi,dimi),eigenval(dimi),x(dimi,dimi))
!!       allocate(w(dimi),wi(dimi),vr(dimi,dimi),vl(dimi,dimi))
!!   
!!       AT = HS(CBeg:CEnd,CBeg:CEnd)
!!       ST = S(CBeg:CEnd,CBeg:CEnd)
!!   
!!   !    call R_Ortho(0,AT,ST)
!!   
!!       if (MOOrtho) then
!!         print *, "----------------------------------------------"
!!         write(6, '(A100)') "Orthogonal Molecular Orbital"
!!         print *, "----------------------------------------------"
!!   
!!   w = 0.0d0
!!   wi=0.0d0
!!   vr=0.0d0
!!   vl=0.0d0
!!       call DVDiag( AT, w,wi,VR,VL, info )
!!   if (info /= 0 ) print*, "Warning, INFO /= ", info,  "in dH_MO, RSDiag subroutine"
!!   
!!   !    call RNormCond(XT)
!!   
!!   !X = matmul(transpose(VL),VL)
!!   
!!   !X=transpose(VL)
!!   X=VR
!!   
!!   !do i=1,dimi
!!   !do j=1,dimi
!!   !if (I/=j) then
!!   !if (dabs(X(i,j)-X(j,i)) > 1.0d-10) print*, "ERROR 01",X(i,j),X(j,i)
!!   !end if
!!   !end do
!!   !if (X(i,i) /= 1.0d0) print*, "ERROR 02", i, X(i,i)
!!   !end do
!!   
!!       call RNormCond(X)
!!   !call CND2(dimi,VR,VL)
!!   
!!   !X = matmul(transpose(VL),VR)
!!   !do i=1,dimi
!!   !do j=1,dimi
!!   !if (I/=j) then
!!   !if (dabs(X(i,j)-X(j,i)) > 1.0d-10) print*, "ERROR 04"
!!   !end if
!!   !end do
!!   !if (X(i,i) /= 1.0d0) print*, "ERROR 03", i, X(i,i)
!!   !end do
!!   
!!   
!!   !AT = matmul(transpose(VL),matmul(AT,VR))
!!   !AT = 0.0d0
!!   !do i=1,dimi
!!   !  AT(i,i) = w(i)
!!   !end do
!!   !AT=0.0d0
!!   !AT = matmul(transpose(VL),matmul(AT,VR))
!!   
!!   !do i=1,dimi
!!   !  w(i) = AT(i,i) 
!!   !end do
!!   
!!   do i=1,dimi
!!     do j=1,dimi
!!       XT(i,j) = X(i,j)*X(i,j)
!!     end do
!!   end do
!!   
!!   do i=1,dimi
!!     If (w(i)+shift > EW1*1.0 .and. w(i)+shift < EW2*1.0) then
!!       print *, "----------------------------------------------"
!!       write(6, '(A20,I6,A20,10F20.10)') "Mol. Orb. = ",i,"Eigenvalue - \muX = ",w(i),w(i)+(chempot-shift),wi(i)
!!       print *, "----------------------------------------------"
!!       do j=1,dimi
!!         write(6, '(I6,10F20.10)') j, XT(j,i)*100
!!       end do
!!     end if
!!   end do
!!   
!!   ! Print Molecular Orbital
!!   if (is == 1)  open (7,file=trim(jobname)//'.MOOU.dat')
!!   if (is == 2)  open (7,file=trim(jobname)//'.MOOD.dat')
!!   Do J=0,dimi
!!     if (J == 0) then
!!       write (7,*) (w(i), i=1,dimi)
!!     else
!!       call flush(7)
!!       write (7,*) (X(j,i), i=1,dimi)
!!     end if
!!   END DO
!!   close(7)
!!   
!!       end if
!!   
!!       if (MONortho) then
!!         print *, "----------------------------------------------"
!!         write(6, '(A100)') "No-Orthogonal Molecular Orbital"
!!         print *, "----------------------------------------------"
!!   
!!   HS = A(is,:,:)+VSh!+Vst
!!         AT = HS(CBeg:CEnd,CBeg:CEnd)
!!         ST = S(CBeg:CEnd,CBeg:CEnd)
!!   
!!         call R_Ortho(0,AT,ST)
!!   
!!         XT = AT
!!   
!!         call RSDiag( XT, eigenval, info )
!!         if (info /= 0 ) print*, "Warning, INFO /= ", info,  "in MO_input, RSDiag subroutine"
!!   
!!         call RM_Ortho(0,XT,ST)
!!         call RNormCond(XT)
!!   
!!         AT = XT
!!   
!!         do i=1,dimi
!!           do j=1,dimi
!!             XT(i,j) = AT(i,j)*AT(i,j)
!!           end do
!!         end do
!!   
!!         do i=1,dimi
!!           If (eigenval(i)+shift > EW1*1 .and. eigenval(i)+shift < EW2*1) then
!!             print *, "----------------------------------------------"
!!             write(6, '(A20,I6,A20,2F20.10)') "Mol. Orb. = ",i,"Eigenvalue - \mu = ",eigenval(i),eigenval(i)+(chempot-shift)
!!             print *, "----------------------------------------------"
!!             do j=1,dimi
!!               write(6, '(I6,10F20.10)') j, XT(j,i)*100
!!             end do
!!           end if
!!         end do
!!   
!!         XT = AT
!!         ! Print Molecular Orbital
!!         if (is == 1)  open (7,file=trim(jobname)//'.MONU.dat')
!!         if (is == 2)  open (7,file=trim(jobname)//'.MOND.dat')
!!         Do J=0,dimi
!!           if (J == 0) then
!!             write (7,*) (eigenval(i), i=1,dimi)
!!           else
!!             call flush(7)
!!             write (7,*) (XT(j,i), i=1,dimi)
!!           end if
!!         END DO
!!         close(7)
!!   
!!       end if !MO, NMO
!!       deallocate(AT,ST,XT,eigenval)
!!     end do !Ncorr
!!   end do !Spin loop
!!   
!!   end subroutine MO_input2
!!   
!!   subroutine dH_MO(A,X,spin,shift)
!!   use parameters, only: EW1, EW2
!!   use numeric, only: RSDiag
!!   use constants
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:), intent(inout) :: A
!!   real*8, dimension(:,:), intent(out) :: X
!!   integer, intent(in) :: spin
!!   real*8, intent(in) :: shift
!!   
!!   real*8, dimension(size(A,1)) :: eigenval
!!   integer :: dimi, info, i,j
!!   real*8, dimension(size(A,1),size(A,2)) :: XT
!!   
!!   dimi = size(A,1)
!!   call RSDiag( A, eigenval, info )
!!   if (info /= 0 ) print*, "Warning, INFO /= ", info,  "in dH_MO, RSDiag subroutine"
!!   
!!   X = A
!!   
!!   A = d_zero
!!   do i=1,dimi
!!     A(i,i) = eigenval(i)
!!   end do
!!   
!!   do i=1,dimi
!!     do j=1,dimi
!!       XT(i,j) = X(i,j)*X(i,j)
!!     end do
!!   end do
!!   
!!   do i=1,dimi
!!     If (eigenval(i)+shift > EW1*100.0 .and. eigenval(i)+shift < EW2*100.0) then
!!       print *, "----------------------------------------------"
!!       write(6, '(A20,I6,A20,2F20.10)') "Mol. Orb. = ",i,"Eigenvalue - \muX = ",eigenval(i)-shift,eigenval(i)
!!       print *, "----------------------------------------------"
!!       do j=1,dimi
!!         write(6, '(I6,10F20.10)') j, XT(j,i)*100
!!       end do
!!     end if
!!   end do
!!   
!!   ! Print Molecular Orbital
!!   if (spin == 1)  open (7,file=trim(jobname)//'.MOOU.dat')
!!   if (spin == 2)  open (7,file=trim(jobname)//'.MOOD.dat')
!!   Do J=0,dimi
!!     if (J == 0) then
!!       write (7,*) (eigenval(i), i=1,dimi)
!!     else
!!       call flush(7)
!!       write (7,*) (X(j,i), i=1,dimi)
!!     end if
!!   END DO
!!   close(7)
!!   
!!   end subroutine dH_MO
!!   
!!   subroutine Frozen_Dens_Core(HD,PD,SD)
!!   use parameters, only: EW1, EW2,Eltype
!!   use parameters, only: NCorrBlM, CorrBegM, CorrEndM,NCorrBl,CorrBeg,CorrEnd
!!   use parameters, only: NQPCore, QPC,coreend
!!   use constants
!!   use Cluster, only : hiaorbno, loaorbno
!!   use numeric, only : RM_Ortho,RSDiag,R_ortho
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:,:),intent(in) :: HD
!!   real*8, dimension(:,:,:),intent(in) :: PD
!!   real*8, dimension(:,:),intent(in) :: SD
!!   
!!   integer :: Nspin,i,j,k,l,CBeg,CEnd,is,dimi,info
!!   real*8, allocatable :: AT(:,:),ST(:,:),PT(:,:),XT(:,:),W(:)
!!   real*8 :: tChar0,tChar,tchar1,dimi2
!!   real*8,dimension(size(PD,2),size(PD,3)) :: PTT
!!   
!!   NSpin = size(HD,1)
!!   dimi2 = size(HD,2)
!!   
!!   open (7,file=trim(jobname)//'.PCore.dat')
!!   
!!   if (NCorrBlM < 1) NCorrBlM = 1
!!   
!!   do is=1,NSpin
!!     do l = 1,NCorrBlM !Molecular Correlated Blocks
!!   
!!       IF (ElType(1) == "GHOST" .and. ElType(2) == "GHOST") then
!!         CBeg = CorrBeg(l)
!!         CEnd = CorrEnd(l)
!!       ELSE
!!         CBeg = loaorbno(CorrBegM(l))
!!         CEnd = hiaorbno(CorrEndM(l))
!!       END IF
!!   
!!   print*, "CBeg, CEnd = ", CBeg, CEnd
!!   
!!       dimi=CEnd-CBeg+1
!!   
!!       allocate(PT(dimi,dimi),AT(dimi,dimi),ST(dimi,dimi),XT(dimi,dimi),W(dimi))
!!       AT = HD(is,cBeg:CEnd,cBeg:CEnd)
!!       PT = PD(is,cBeg:CEnd,cBeg:CEnd)
!!       ST = SD(cBeg:CEnd,cBeg:CEnd)
!!   
!!   !    PT = matmul(PT,ST) !Charge Matrix
!!       PTT = matmul(PD(is,:,:),SD) !Charge Matrix
!!   
!!       tchar0 = d_zero
!!       do i =1,dimi2
!!         tchar0 = tchar0 + PTT(i,i)
!!       end do
!!   
!!       write(6,'(A20,F20.10)') "----- Tr Total = ", tchar0
!!   
!!       XT = AT
!!   
!!       call R_Ortho(0,XT,ST)
!!       call RSDiag(XT,W,info)
!!       call RM_Ortho(0,XT,ST)
!!   
!!       ! Computing the density Matrix
!!   
!!       do j=1,dimi
!!         do i=1,dimi
!!           W=0.0d0
!!           do k=1,dimi
!!             W(k)=XT(i,k)*XT(j,k)
!!           end do
!!           PT(i,j)=0.0d0
!!           do k=1,QPC(1)-1
!!             PT(i,j)=PT(i,j)+W(k)
!!           end do
!!         end do
!!       end do
!!   
!!       XT = PT
!!       PT = matmul(PT,ST) !Charge Matrix
!!   
!!       tchar0 = d_zero
!!       do i =1,dimi
!!         tchar0 = tchar0 + PT(i,i)
!!       end do
!!   
!!       write(6,'(A20,F20.10)') "----- Tr 0 = ", tchar0
!!   
!!   !    PT = matmul(transpose(XT),matmul(PT,XT))   
!!   !    PT = matmul(PT,ST) !Charge Matrix
!!   
!!       tchar = d_zero
!!       write(6,'(A100)') "--------------------- Density Matrix in Active Subshell ---------------"
!!       do j =1,NQPcore
!!         tchar = tchar + PT(QPC(j),QPC(j))
!!         write(6,'(A20,I6,A20,F20.10)') "Mol. Orb. ", QPC(j), "Den. Mat. =  ", PT(QPC(j))
!!       end do
!!   
!!       PT = XT
!!   !    tchar1 = tchar0-tchar
!!   
!!   !    do i=1,dimi
!!   !      do j = 1, NQPCore
!!   !        if (i == QPC(j) )  then
!!   !          PT(i,i) = d_zero
!!   !        end if
!!   !      end do
!!   !    end do
!!   
!!   !    PT = matmul(XT,matmul(PT,transpose(XT)))
!!   
!!   !    tchar0 = d_zero
!!   !    do i =1,dimi
!!   !      tchar0 = tchar0 + PT(i,i)
!!   !    end do
!!   
!!   !    write(6,'(A20,F20.10)') "----- Tr END = ", tchar0
!!   
!!   !    If (tchar0 /= tchar1) print*, "War: tchar0 /= tchar1 = ", tchar0, tchar1
!!   
!!       do i=1,dimi
!!         do j=1,i
!!   !        if (i  < CBeg .or. j < CBeg .or. i > CEnd .or. j > CEnd) then
!!   !          call flush(7)
!!   !          write(7, '(3I6, E20.10)') is, i,j,PD(is,i,j)
!!   !        else
!!             call flush(7)
!!             write(7, '(3I6, E20.10)') is, i,j,PT(i,j)
!!   !        end if
!!         end do
!!       end do
!!   
!!       deallocate(PT,AT,XT,ST,W)
!!     end do
!!   end do
!!   
!!   close(7)
!!   
!!   !if (Coreend) STOP
!!   end subroutine Frozen_Dens_core
!!   
!!   Subroutine XFrozen(HD,SD,shift,jjcy)
!!   use parameters, only: EW1, EW2,Eltype
!!   use parameters, only: NCorrBlM, CorrBegM, CorrEndM,NCorrBl,CorrBeg,CorrEnd
!!   use parameters, only: NQPCore, QPC,coreend,COREA
!!   use constants
!!   use Cluster, only : hiaorbno, loaorbno
!!   use numeric, only : RM_Ortho,RSDiag,R_ortho
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:,:),intent(inout) :: HD
!!   real*8, dimension(:,:),intent(in) :: SD
!!   real*8, intent(in) :: shift
!!   integer, intent(in) :: jjcy
!!   
!!   integer :: Nspin,i,j,k,l,CBeg,CEnd,is,dimi,info
!!   real*8, allocatable :: AT(:,:),ST(:,:),PT(:,:),XT(:,:),W(:),XT0(:,:),W0(:)
!!   real*8 :: shift2
!!   
!!   shift2=4.50474
!!   
!!   NSpin = size(HD,1)
!!   
!!   if (NCorrBlM < 1) NCorrBlM = 1
!!   
!!   do is=1,NSpin
!!     do l = 1,NCorrBlM !Molecular Correlated Blocks
!!   
!!       IF (ElType(1) == "GHOST" .and. ElType(2) == "GHOST") then
!!         CBeg = CorrBeg(l)
!!         CEnd = CorrEnd(l)
!!       ELSE
!!         CBeg = loaorbno(CorrBegM(l))
!!         CEnd = hiaorbno(CorrEndM(l))
!!       END IF
!!   
!!       dimi=CEnd-CBeg+1
!!   
!!       allocate(PT(dimi,dimi),AT(dimi,dimi),ST(dimi,dimi),XT(dimi,dimi),W(dimi))
!!       allocate(XT0(dimi,dimi),W0(dimi))
!!   
!!       AT = HD(is,cBeg:CEnd,cBeg:CEnd)
!!       ST = SD(cBeg:CEnd,cBeg:CEnd)
!!   
!!       ! Read Molecular Orbital
!!       if (is == 1)  open (8,file=trim(jobname)//'.XCoreU.dat')
!!       if (is == 2)  open (8,file=trim(jobname)//'.XCoreD.dat')
!!       Do J=0,dimi
!!         if (J == 0) then
!!           read(8,*) (W0(i), i=1,dimi)
!!         else
!!           read(8,*) (XT0(j,i), i=1,dimi)
!!         end if
!!       END DO
!!       close(8)
!!   
!!       XT = AT
!!   
!!       call R_Ortho(0,XT,ST)
!!       call RSDiag(XT,W,info)
!!   
!!       PT = d_zero
!!       do i =1,dimi
!!         if (COREA) then
!!           if (i < QPC(1)) then
!!             PT(i,i) = W0(i)+shift2-shift
!!           else if (i <= QPC(NQPCore) .and. jjcy <= 6) then
!!             PT(i,i) = W0(i)+shift2-shift
!!             print*, i,W0(i)+shift2, W(i)+shift
!!   !          PT(i,i) = W(i)
!!           else
!!             PT(i,i) = W(i)
!!           end if
!!         else
!!           PT(i,i) = W(i)
!!         end if
!!       end do
!!   
!!       AT = matmul(XT0,matmul(PT,transpose(XT0)))
!!   
!!       call R_Ortho(1,AT,ST)
!!       
!!       HD(is,cBeg:CEnd,cBeg:CEnd) = AT
!!   
!!       deallocate(PT,AT,XT,ST,W,XT0,W0)
!!     end do
!!   end do
!!   
!!   end subroutine XFrozen
!!   
!!   subroutine dH_channel(A,X,spin,B,Y,N0,energy)
!!   use constants
!!   use numeric, only : RSNCond
!!   
!!   implicit none
!!   
!!   complex*16, dimension(:,:), intent(in) :: A,B
!!   real*8, dimension(:), intent(in) :: X,Y
!!   integer, intent(in) :: spin,N0
!!   real*8, intent(in) :: energy
!!   
!!   integer :: dimi, i,j,k,l
!!   complex*16, dimension(size(A,1),size(A,2)) :: XT,YT!,YB
!!   real*8, dimension(size(A,1),size(A,2)) :: YA
!!   
!!   dimi=size(A,1)
!!   
!!   do i=1,dimi
!!     do j=1,dimi
!!       XT(i,j) = dabs(A(i,j))*dabs(A(i,j))
!!       YT(i,j) = dabs(B(i,j))*dabs(B(i,j))
!!       YA(i,j) = dreal(XT(i,j))+dreal(YT(i,j))
!!     end do
!!   end do
!!   
!!   call RSNCond(YA)
!!   
!!   if (N0 == 0) then
!!     do i=1,dimi
!!       If (x(i) >= 1.0d-10 .or. x(i) <= -1.0d-10) then
!!         print *, "-----------------------------------------------------------------------------"
!!         write(6, '(I6,A20,I6,A20,2E20.10)') spin,"EigenChannel = ",i,"Transmission= ",x(i),y(i)
!!         print *, "-----------------------------------------------------------------------------"
!!         write(6, '(A8,A10,2A20,2A10)') " -N- ","-%AC_TCh- ","  ------ TS ------  ","  ------ TSI ------  "," -%AC-  ","-%AC(TSi)-"
!!         do j=1,dimi
!!           write(6, '(I6,100F10.5)') j, YA(j,i)*100,A(j,i),B(j,i),dreal(XT(j,i))*100,dreal(YT(j,i))*100
!!         end do
!!       end if
!!     end do
!!   end if
!!   
!!   write(169,'(3000F10.5)') energy, (YA(j,dimi)*100, j=1,dimi),(YA(j,dimi-1)*100, j=1,dimi),(YA(j,dimi-2)*100, j=1,dimi),(YA(j,dimi-3)*100, j=1,dimi)
!!   
!!   end subroutine dH_channel
!!   
!!   
!!   subroutine Red_Den(HD,PD,SD,nel)
!!   !use parameters, only: EW1, EW2,Eltype
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     !use parameters, only: NQPCore, QPC,coreend
!!     use constants
!!     use numeric, only : RM_Ortho,RSDiag,R_ortho
!!   
!!   implicit none
!!   
!!   real*8, dimension(:,:,:),intent(in) :: HD
!!   real*8, dimension(:,:,:),intent(in) :: PD
!!   real*8, dimension(:,:),intent(in) :: SD
!!   real*8, intent(in) :: nel
!!   
!!   integer :: Nspin,i,j,k,l,MBeg,MEnd,is,dimi,info,neli,neli0
!!   real*8, allocatable :: AT(:,:),ST(:,:),PT(:,:),XT(:,:),W(:)
!!   real*8 :: tChar0,tChar,tchar1,dimi2
!!   real*8,dimension(size(PD,2),size(PD,3)) :: PTT
!!   
!!   NSpin = size(HD,1)
!!   dimi2 = size(HD,2)
!!   neli0 = nel*NSpin/2.0
!!   
!!   print*, "neli = ", neli0
!!   
!!   do is=1,NSpin
!!   
!!       if (AORB ) then
!!         MBeg = SubSpBeg
!!         MEnd = SubSpEnd
!!       else
!!         MBeg = loaorbno(SubSpBeg)
!!         MEnd = hiaorbno(SubSpEnd)
!!       end if
!!   
!!       dimi=MEnd-MBeg+1
!!   
!!       allocate(PT(dimi,dimi),AT(dimi,dimi),ST(dimi,dimi),XT(dimi,dimi),W(dimi))
!!       AT = HD(is,mBeg:mEnd,mBeg:mEnd)
!!       PT = PD(is,mBeg:mEnd,mBeg:mEnd)
!!       ST = SD(mBeg:mEnd,mBeg:mEnd)
!!   !    PT = matmul(PT,ST) !Charge Matrix
!!       PTT = matmul(PD(is,:,:),SD) !Charge Matrix
!!   
!!       tchar0 = d_zero
!!       do i =1,dimi2
!!         tchar0 = tchar0 + PTT(i,i)
!!       end do
!!   
!!       write(6,'(A20,F20.10)') "----- Tr Total = ", tchar0*(2.0/Nspin)
!!   
!!       XT = AT
!!   
!!       call R_Ortho(0,XT,ST)
!!       call RSDiag(XT,W,info)
!!       call RM_Ortho(0,XT,ST)
!!   
!!       ! Computing the density Matrix
!!   
!!   do l=1,7
!!   neli=neli0-4+l
!!       do j=1,dimi
!!         do i=1,dimi
!!           W=0.0d0
!!           do k=1,dimi
!!             W(k)=XT(i,k)*XT(j,k)
!!           end do
!!           PT(i,j)=0.0d0
!!           do k=1,neli
!!             PT(i,j)=PT(i,j)+W(k)
!!           end do
!!         end do
!!       end do
!!   
!!       AT = d_zero
!!       AT = matmul(PT,ST) !Charge Matrix
!!   
!!       tchar0 = d_zero
!!       do i =1,dimi
!!         tchar0 = tchar0 + AT(i,i)
!!       end do
!!   
!!       write(6,'(A100)') "--------------------- Density Matrix in Active Subshell---------------"
!!       do j =1,dimi
!!   !      tchar = tchar + AT(j,j)
!!         write(6,'(A20,I6,A20,F20.10)') "Mol. Orb. ", j, "Den. Mat. = ",AT(j,j)*2.0/Nspin
!!       end do
!!       write(6,'(A20,F20.10,I6)') "----- Tr 0 = ", (2.0d0/Nspin)*tchar0,neli0-neli
!!   
!!   
!!   end do
!!   
!!   
!!       deallocate(PT,AT,XT,ST,W)
!!   end do
!!   
!!   end subroutine Red_Den
!!   
!!   !Compute the partial density matrix.
!!   !Que devuleva la P en la nueva métrica lista para sumar
!!   !Devuelve PS en la nueva metrica
!!   subroutine OrthProjMatMix(SD,A,B,NAOrbs) !Para A tilde, Mol Tilde
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB
!!     use constants, only: d_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: RMatPow,R_DDbase,CMatPow
!!   
!!     implicit none
!!   
!!     real*8, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: B
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, intent(in) :: NAOrbs
!!   
!!     real*8, dimension(size(A,1),size(A,2)) :: D
!!     real*8, dimension(size(B,1),size(B,2)) :: E
!!     real*8, dimension(NAOrbs,NAOrbs) :: C,SDnh,SDph,InvSD
!!   !  complex*16, dimension(NAOrbs,NAOrbs) :: S1
!!     integer :: MBeg,MEnd,i,j,k,l,IntN,N
!!   
!!     real*8 :: QMolX0
!!   
!!       if (AORB ) then
!!         MBeg = SubSpBeg
!!         MEnd = SubSpEnd
!!       else
!!         MBeg = loaorbno(SubSpBeg)
!!         MEnd = hiaorbno(SubSpEnd)
!!       end if
!!   
!!     N = IntTyp
!!     D = A
!!     E = B
!!   
!!     If (N == 0) IntN = 1
!!     If (N == 1) IntN = 0
!!     If (N == -1) IntN = -10
!!     If (N == -10) IntN = -1
!!     If (N == 10) IntN = 10
!!   
!!     call RMatPow( SD, -1.0d0, InvSD )
!!     call RMatPow( SD, -0.5d0, SDnh )
!!     call RMatPow( SD, 0.5d0, SDph )
!!   
!!     if (N == 10) then
!!   
!!       D = matmul(D,SD)
!!       E = matmul(E,SD)
!!       C = D
!!       D = 0.0d0
!!       if (MBeg > 1 .and. NAORBS > MEnd) then
!!         D(1:Mbeg-1,1:Mbeg-1) = E(1:Mbeg-1,1:Mbeg-1)
!!         D(Mend+1:NAOrbs,Mend+1:NAOrbs) = E(Mend+1:NAOrbs,Mend+1:NAOrbs)
!!         D(Mbeg:Mend,Mbeg:Mend) = C(Mbeg:Mend,Mbeg:Mend)
!!       else if (MBeg > 1 .and. NAORBS == MEnd) then
!!         D(Mbeg:Mend,Mbeg:Mend) = C(Mbeg:Mend,Mbeg:Mend)
!!         D(1:Mbeg-1,1:Mbeg-1) = E(1:Mbeg-1,1:Mbeg-1)
!!       else if (MBeg == 1 .and. NAORBS > MEnd) then
!!         D(Mbeg:Mend,Mbeg:Mend) = C(Mbeg:Mend,Mbeg:Mend)
!!         D(Mend+1:NAOrbs,Mend+1:NAOrbs) = E(Mend+1:NAOrbs,Mend+1:NAOrbs)
!!       end if
!!   
!!     else
!!   
!!       if (N == 0) D = matmul(SD,matmul(D,SD)) !out
!!       if (N == 0) E = matmul(SD,matmul(E,SD)) !Forz
!!   
!!       call R_DDbase(N,D,SD,MBeg,MEnd)
!!   
!!       C = D !out
!!       call R_DDbase(N,E,SD,MBeg,MEnd)
!!   
!!       D = E !Froz
!!       do i = Mbeg,Mend
!!         do j = Mbeg,Mend
!!           D(i,j) = C(i,j)
!!         end do
!!       end do
!!   
!!       call R_DDbase(IntN,D,SD,MBeg,MEnd)
!!   
!!       if (N == 0) D = matmul(InvSD,matmul(D,InvSD))
!!   
!!       A = D
!!   
!!     end if
!!   
!!   end subroutine OrthProjMatMix
!!   subroutine PALSORMix(SD,A,B,PSV) !Para A tilde, Mol Tilde
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB
!!     use constants, only: d_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: RMatPow,R_DDbase,R_PALSORT,CMatPow
!!   
!!     implicit none
!!   
!!     real*8, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: B
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, dimension(:), intent(in) :: PSV
!!   
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnew,InvSDnew,SDphnew,SDnhnew
!!     real*8, dimension(size(B,1),size(B,2)) :: SDnh,SDph,InvSD,C,D,E
!!     integer :: i,j,k,l,Tdim
!!     complex*16, dimension(size(A,1),size(A,2)) :: S1
!!   
!!     Tdim = size(A,2)
!!   
!!     call RMatPow( SD, -1.0d0, InvSD )
!!     call RMatPow( SD, -0.5d0, SDnh )
!!     call RMatPow( SD, 0.5d0, SDph )
!!   
!!     ! SDnew in Thy Base
!!     SDnew = SD
!!     call R_PALSORT(SDnew,SD,PSV)
!!   
!!     S1 = Sdnew
!!     call CMatPow( S1, -1.0d0, S1 ) !Inv
!!     InvSDnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, 0.5d0, S1 ) !+0.5
!!     SDphnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, -0.5d0, S1 ) !+0.5
!!     SDnhnew = dreal(S1)
!!   
!!     C = A !new dual
!!     D = B !frozen dual
!!   
!!     C = matmul(SDph,matmul(C,SDph)) 
!!     C = matmul(SDphnew,matmul(C,SDphnew))
!!     
!!     D = matmul(SDph,matmul(D,SDph))
!!     D = matmul(SDphnew,matmul(D,SDphnew))
!!   
!!     E = D
!!     do i = 1,Tdim
!!       do j = 1,Tdim
!!         if (PSV(i) == 1 .and. PSV(j) == 1) E(i,j) = C(i,j)
!!       end do
!!     end do
!!   
!!     E = matmul(SDnhnew,matmul(E,SDnhnew)) !ortho
!!     E = matmul(SDnh,matmul(E,SDnh)) !dual
!!   
!!     A = E
!!   
!!   end subroutine PALSORMix
!!   subroutine OrthProjMat(SD,A,NAOrbs) !Para A tilde, Mol Tilde
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB,ProjQ
!!     use constants, only: d_zero,c_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: CMatPow,R_DDbase,RMatPow
!!   
!!     implicit none
!!   
!!     real*8, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, intent(in) :: NAOrbs
!!   
!!     real*8, dimension(NAOrbs,NAOrbs) :: InvSDnew,SDphnew,SDnew,InvSD,B
!!     integer :: MBeg,MEnd,N
!!     complex*16, dimension(NAOrbs,NAOrbs) :: S1
!!   
!!       if (AORB ) then
!!         MBeg = SubSpBeg
!!         MEnd = SubSpEnd
!!       else
!!         MBeg = loaorbno(SubSpBeg)
!!         MEnd = hiaorbno(SubSpEnd)
!!       end if
!!   
!!     N = IntTyp
!!     call RMatPow( SD, -1.0d0, InvSD )
!!   
!!     if (IntTyp >= 0 .and. IntTyp /= 10) then
!!   
!!       SDnew = SD
!!       call R_DDbase(N,SDnew,SD,MBeg,MEnd)
!!   
!!       S1 = c_zero
!!       S1 = Sdnew
!!       call CMatPow( S1, -1.0d0, S1 ) !Inv
!!   !    call RMatPow( SDnew, -1.0d0, InvSDnew ) !Inv
!!       InvSDnew = dreal(S1)
!!       S1 = c_zero
!!       S1 = Sdnew
!!       call CMatPow( S1, 0.5d0, S1 ) !+0.5
!!   !    call RMatPow( SDnew, 0.5d0, SDphnew ) !+0.5
!!       SDphnew = dreal(S1)
!!     
!!     else if (IntTyp /= 10 .and. IntTyp < 0) then
!!   
!!       SDnew = InvSD
!!       call R_DDbase(N,SDnew,SD,MBeg,MEnd)
!!   
!!       S1 = c_zero
!!       S1 = Sdnew
!!       call CMatPow( S1, -1.0d0, S1 ) ! SD
!!       InvSDnew = dreal(S1)
!!   !    S1 = Sdnew
!!       call CMatPow( S1, 0.5d0, S1 ) ! 0.5
!!       SDphnew = dreal(S1)
!!   
!!     end if
!!   
!!     if (N == 10) then
!!   
!!       A = matmul(A,SD)
!!   
!!     else
!!   
!!       if (N == 0) A = matmul(SD,matmul(A,SD))
!!   
!!       call R_DDbase(N,A,SD,MBeg,MEnd)
!!   
!!       if (N == 0) A = matmul(InvSDnew,matmul(A,InvSDnew))
!!   
!!       A = matmul(SDphnew,matmul(A,SDphnew))
!!   
!!     end if
!!   
!!   end subroutine OrthProjMat
!!   subroutine OrthProjQ(SD,A,NAOrbs) !Para A tilde, Mol Tilde
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB,ProjQ
!!     use constants, only: d_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: CMatPow,R_DDbase,RMatPow
!!   
!!     implicit none
!!   
!!     real*8, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, intent(in) :: NAOrbs
!!   
!!     real*8, dimension(NAOrbs,NAOrbs) :: InvSDnew,SDphnew,SDnew,InvSD,B
!!     integer :: MBeg,MEnd,N
!!     complex*16, dimension(NAOrbs,NAOrbs) :: S1
!!   
!!       if (AORB ) then
!!         MBeg = SubSpBeg
!!         MEnd = SubSpEnd
!!       else
!!         MBeg = loaorbno(SubSpBeg)
!!         MEnd = hiaorbno(SubSpEnd)
!!       end if
!!   
!!   B = 0.0d0
!!       B(Mbeg:Mend,Mbeg:Mend)=A(Mbeg:Mend,Mbeg:Mend)
!!   A = 0.0d0
!!       A(Mbeg:Mend,Mbeg:Mend)=SD(Mbeg:Mend,Mbeg:Mend)
!!   
!!   A = matmul(B,A)
!!   !return
!!   
!!   !  N = IntTyp
!!   !  call RMatPow( SD, -1.0d0, InvSD )
!!   
!!   !    SDnew = SD
!!   !    call R_DDbase(N,SDnew,SD,MBeg,MEnd)
!!   
!!   !    S1 = Sdnew
!!   !    call CMatPow( S1, -1.0d0, S1 ) !Inv
!!   !    InvSDnew = dreal(S1)
!!   !    S1 = Sdnew
!!   !    call CMatPow( S1, 0.5d0, S1 ) !+0.5
!!   !    SDphnew = dreal(S1)
!!   
!!   !    B=A
!!   !    if (N == 0) N=-1!A = matmul(SD,matmul(A,SD))
!!   !    call R_DDbase(-1,B,SD,MBeg,MEnd)
!!   !    N=IntTyp
!!   !    if (N == 0) N=0!A = matmul(InvSDnew,matmul(A,InvSDnew))
!!   !    B = matmul(SDphnew,matmul(B,SDphnew))
!!   
!!   !    A = matmul(SD,matmul(A,SD))
!!   !    call R_DDbase(N,A,SD,MBeg,MEnd)
!!   
!!   !    A = matmul(InvSDnew,matmul(A,InvSDnew))
!!   
!!   !    A = matmul(SDphnew,matmul(A,SDphnew))
!!   
!!   !    A(Mbeg:Mend,Mbeg:Mend)=B(Mbeg:Mend,Mbeg:Mend)
!!   
!!   end subroutine OrthProjQ
!!   subroutine PALSORB(SD,A,PSV) !Para A tilde, Mol Tilde
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB
!!     use constants, only: d_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: CMatPow,R_DDbase,RMatPow,R_PALSORT
!!   
!!     implicit none
!!   
!!     real*8, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, dimension(:), intent(in) :: PSV
!!   
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnew,InvSDnew,SDphnew,SDnhnew
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnh,SDph,InvSD,C
!!     integer :: i,j,k,l,Tdim
!!     complex*16, dimension(size(A,1),size(A,2)) :: S1
!!   
!!     call RMatPow( SD, -1.0d0, InvSD )
!!     call RMatPow( SD, -0.5d0, SDnh )
!!     call RMatPow( SD, 0.5d0, SDph )
!!   
!!     ! SDnew in Thy Base
!!     SDnew = SD
!!     call R_PALSORT(SDnew,SD,PSV)
!!   
!!     S1 = Sdnew
!!     call CMatPow( S1, -1.0d0, S1 ) !Inv
!!     InvSDnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, 0.5d0, S1 ) !+0.5
!!     SDphnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, -0.5d0, S1 ) !+0.5
!!     SDnhnew = dreal(S1)
!!   
!!     C = A !dual
!!   
!!     C = matmul(SDph,matmul(C,SDph))
!!     C = matmul(SDphnew,matmul(C,SDphnew)) !A NOortho
!!     A = matmul(SDnhnew,matmul(C,SDnhnew)) !A ortho
!!   
!!   end subroutine PALSORB
!!   subroutine C_PALSOR(N1,N2,SD,A,PSV) !Para A tilde, Mol Tilde
!!   !N1 = 0 A in NO(in) =N2 (out)
!!   !N1 =1 A in dual
!!   !N1 = 2 ortho
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB
!!     use constants, only: d_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: CMatPow,R_DDbase,RMatPow,R_PALSORT
!!   
!!     implicit none
!!   
!!     complex*16, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, dimension(:), intent(in) :: PSV
!!     integer, intent(in) :: N1,N2
!!   
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnew,InvSDnew,SDphnew,SDnhnew
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnh,SDph,InvSD
!!     complex*16, dimension(size(A,1),size(A,2)) :: C
!!     integer :: i,j,k,l,Tdim
!!     complex*16, dimension(size(A,1),size(A,2)) :: S1
!!   
!!     call RMatPow( SD, -1.0d0, InvSD )
!!     call RMatPow( SD, -0.5d0, SDnh )
!!     call RMatPow( SD, 0.5d0, SDph )
!!   
!!     ! SDnew in Thy Base
!!     SDnew = SD
!!     call R_PALSORT(SDnew,SD,PSV)
!!   
!!     S1 = Sdnew
!!     call CMatPow( S1, -1.0d0, S1 ) !Inv
!!     InvSDnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, 0.5d0, S1 ) !+0.5
!!     SDphnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, -0.5d0, S1 ) !+0.5
!!     SDnhnew = dreal(S1)
!!   
!!     C = A
!!   
!!     if (N1 == 0 ) then
!!       C = matmul(SDnh,matmul(C,SDnh)) !C ortho
!!       if (N2 == 0) A = matmul(SDphnew,matmul(C,SDphnew)) !A NO
!!       if (N2 == 1) A = matmul(SDnhnew,matmul(C,SDnhnew)) !A Dual
!!       if (N2 == 2) A = C !A ortho
!!     else if (N1 == 1) then
!!       C = matmul(SDph,matmul(C,SDph)) !C ortho
!!       if (N2 == 0) A = matmul(SDphnew,matmul(C,SDphnew)) !A NO
!!       if (N2 == 1) A = matmul(SDnhnew,matmul(C,SDnhnew)) !A Dual
!!       if (N2 == 2) A = C !A ortho
!!     else if (N2 == 2) then
!!       if (N2 == 0) A = matmul(SDphnew,matmul(C,SDphnew)) !A NO
!!       if (N2 == 1) A = matmul(SDnhnew,matmul(C,SDnhnew)) !A Dual
!!       if (N2 == 2) A = C !A ortho
!!     end if
!!   
!!   end subroutine C_PALSOR
!!   subroutine R_PALSOR(N1,N2,SD,A,PSV) !Para A tilde, Mol Tilde
!!   !N1 = 0 A in NO(in) =N2 (out)
!!   !N1 =1 A in dual
!!   !N1 = 2 ortho
!!     use parameters, only: SubSpBeg,SubSpEnd,IntTyp,AORB
!!     use constants, only: d_zero
!!     use cluster, only: LoAOrbNo, HiAOrbNo
!!     use numeric, only: CMatPow,R_DDbase,RMatPow,R_PALSORT
!!   
!!     implicit none
!!   
!!     real*8, dimension(:,:), intent(inout) :: A
!!     real*8, dimension(:,:), intent(in) :: SD
!!     integer, dimension(:), intent(in) :: PSV
!!     integer, intent(in) :: N1,N2
!!   
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnew,InvSDnew,SDphnew,SDnhnew
!!     real*8, dimension(size(A,1),size(A,2)) :: SDnh,SDph,InvSD
!!     real*8, dimension(size(A,1),size(A,2)) :: C
!!     integer :: i,j,k,l,Tdim
!!     complex*16, dimension(size(A,1),size(A,2)) :: S1
!!   
!!     call RMatPow( SD, -1.0d0, InvSD )
!!     call RMatPow( SD, -0.5d0, SDnh )
!!     call RMatPow( SD, 0.5d0, SDph )
!!   
!!     ! SDnew in Thy Base
!!     SDnew = SD
!!     call R_PALSORT(SDnew,SD,PSV)
!!   
!!     S1 = Sdnew
!!     call CMatPow( S1, -1.0d0, S1 ) !Inv
!!     InvSDnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, 0.5d0, S1 ) !+0.5
!!     SDphnew = dreal(S1)
!!     S1 = Sdnew
!!     call CMatPow( S1, -0.5d0, S1 ) !+0.5
!!     SDnhnew = dreal(S1)
!!   
!!     C = A
!!   
!!     if (N1 == 0 ) then
!!       C = matmul(SDnh,matmul(C,SDnh)) !C ortho
!!       if (N2 == 0) A = matmul(SDphnew,matmul(C,SDphnew)) !A NO
!!       if (N2 == 1) A = matmul(SDnhnew,matmul(C,SDnhnew)) !A Dual
!!       if (N2 == 2) A = C !A ortho
!!     else if (N1 == 1) then
!!       C = matmul(SDph,matmul(C,SDph)) !C ortho
!!       if (N2 == 0) A = matmul(SDphnew,matmul(C,SDphnew)) !A NO
!!       if (N2 == 1) A = matmul(SDnhnew,matmul(C,SDnhnew)) !A Dual
!!       if (N2 == 2) A = C !A ortho
!!     else if (N2 == 2) then
!!       if (N2 == 0) A = matmul(SDphnew,matmul(C,SDphnew)) !A NO
!!       if (N2 == 1) A = matmul(SDnhnew,matmul(C,SDnhnew)) !A Dual
!!       if (N2 == 2) A = C !A ortho
!!     end if
!!   
!!   end subroutine R_PALSOR
END MODULE MOLMOD
