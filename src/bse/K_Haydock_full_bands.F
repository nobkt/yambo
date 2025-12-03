!
! License-Identifier: GPL
!
! Copyright (C) 2024 The Yambo Team
!
! Authors: Memory-efficient full-band BSE implementation
!
! This module implements the memory-efficient BSE calculation with full band
! consideration for coupling mode (BSEmod="coupling") using biorthogonal
! Lanczos iteration with on-the-fly kernel construction.
!
! Reference: mydoc/bse_full_bands_coupling_specification.md
!            mydoc/bse_full_bands_coupling_design.md
!
subroutine K_Haydock_full_bands(iq,W,OBS,Xk,Ken)
 !
 ! Memory-efficient Haydock solver for full-band BSE calculation
 ! with coupling terms using biorthogonal Lanczos iteration.
 !
 ! The coupling Hamiltonian has the structure:
 !
 ! H_BSE = [ H^res      C     ]
 !         [ -C*    -H^res*   ]
 !
 ! which satisfies pseudo-Hermitian property: H = Theta H^dagger Theta^{-1}
 ! where Theta = diag(I, -I) is the metric operator.
 !
 use pars,           ONLY:SP,schlen,cI,cZERO,cONE
 use functions,      ONLY:NAN
 use R_lattice,      ONLY:q0_def_norm,bz_samp
 use X_m,            ONLY:global_gauge,X_OUTPUT_driver
 use wrapper_omp,    ONLY:V_plus_alpha_V_omp
 use frequency,      ONLY:w_samp
 use com,            ONLY:msg,isec
 use electrons,      ONLY:levels
 use BS_solvers,     ONLY:BSS_desc,Haydock_threshold,io_BSS_Haydock,&
&                         Haydock_vector,BS_HAYVEC_alloc,BS_HAYVEC_free,Haydock_terminator,Haydock_v0_mod, &
&                         Haydock_iterIO,Haydock_iterMAX,HAYVEC_size,BSS_uses_DbGd,FG_factor
 use BS,             ONLY:BS_K_Coupling,BS_Blocks_symmetrize_K,BS_T_grp,BS_nT_grps,BS_K_dim,&
&                         l_BS_ares_from_res,BS_H_dim,BS_res_ares_n_mat,BS_nkFGbz_in_Tgrp
 use LIVE_t,         ONLY:LIVE_message
 use stderr,         ONLY:tty_size,intc,real2ch
 use IO_int,         ONLY:io_control
 use IO_m,           ONLY:OP_RD_CL,OP_WR_CL,VERIFY,REP
 use parser_m,       ONLY:parser
 use parallel_int,   ONLY:PP_bcast
 use parallel_m,     ONLY:PAR_IND_T_groups,PAR_IND_T_Haydock,PAR_COM_T_Haydock,PAR_COM_Xk_ibz_INDEX, &
                          PAR_COM_eh_INDEX
 use timing_m,       ONLY:timing
 use descriptors,    ONLY:IO_desc_add
 !
#include<y_memory.h>
 !
 type(w_samp), intent(in) :: W
 integer,      intent(in) :: iq
 character(*), intent(in) :: OBS
 type(bz_samp),intent(in) :: Xk
 type(levels), intent(in) :: Ken
 !
 ! Work Space
 !
 integer            ::it,it_on_disk,i_g,ID,io_err,i_desc_acc,i_desc_iter
 !
 ! Biorthogonal Lanczos vectors: psi (right) and phi (left)
 type(Haydock_vector),allocatable :: psi_n(:),psi_nm1(:),psi_np1(:),psi_0(:)
 type(Haydock_vector),allocatable :: phi_n(:),phi_nm1(:),phi_np1(:),phi_0(:)
 type(Haydock_vector),allocatable :: H_psi(:),H_phi(:)
 type(Haydock_vector),allocatable :: V_tmp(:)
 !
 complex(SP), pointer :: Vn_f(:)   => null()
 complex(SP), pointer :: Vnm1_f(:) => null()
 complex(SP), pointer :: Vnp1_f(:) => null()
 complex(SP), allocatable :: Cf(:) 
 !
 type(Haydock_vector),allocatable :: W0(:)
 integer    :: i_ColFG,i_T,i_k_bz
 !
 real(SP)           ::reached_threshold
 real(SP)           ::Bf(Haydock_iterMAX+1),Af(Haydock_iterMAX)
 integer            ::ItParity,Cf_size,Haydock_n_converged,Haydock_required_converged
 complex(SP)        ::v0_mod_,b_squared
 character(schlen)  ::dumb_ch
 logical            ::use_left_and_right,velocity_correction,tmp_converged
 logical, external  ::STOP_now
 character(1)       ::sec
 !
 ! Biorthogonality parameters
 real(SP)           ::ortho_threshold,ortho_error
 real(SP), parameter::ORTHO_TOL = 1.0E-12_SP
 real(SP), parameter::ORTHO_WARN = 1.0E-10_SP
 real(SP), parameter::B_MIN = 1.0E-14_SP
 real(SP), parameter::IMAG_TOL = 1.0E-6_SP
 !
 if (BS_res_ares_n_mat==2) call error(" Haydock full bands not coded for case with two BSE matrices")
 !
 ! Initialize & write the output file 
 !====================================
 call X_OUTPUT_driver("DEF INIT",solver="Haydock") 
 call K_OUTPUT(iq,W,"open Haydock close ",trim(OBS))
 !
 call PARALLEL_Haydock_VEC_COMMs('assign')
 !
 call timing('Haydock Full Bands Solver',OPR='start')
 !
 Haydock_n_converged=0
 Haydock_required_converged=2
 !
 call BS_Blocks_symmetrize_K()
 !
 ! Use terminator?
 if (.not.Haydock_terminator) call parser('BSHayTer',Haydock_terminator)
 !
 use_left_and_right=.false.
 use_left_and_right=trim(OBS)=="kerr".or.trim(OBS)=="dichroism"
 velocity_correction= trim(global_gauge)=="velocity" .and. (trim(OBS)=="optics" .or. trim(OBS)=="kerr")
 !
 if (isec(2)/=0) sec="="
 if (isec(2)==0) sec="+"
 call section(sec,'Full-Band Haydock Solver (Biorthogonal) in the '//trim(OBS)//' basis @q'//trim(intc(iq)))
 !
 if (     l_BS_ares_from_res) HAYVEC_size = BS_K_dim(1)
 if (.not.l_BS_ares_from_res) then
   if(.not.BS_K_coupling) HAYVEC_size = BS_K_dim(1)
   if(     BS_K_coupling) HAYVEC_size = BS_H_dim
 endif
 !
 if (BSS_uses_DbGd) then
   call msg('r','Using Double Grid method') 
   call msg('r','Number of Coarse grid BZ kpts : '//trim(intc(Xk%nbz))) 
   call msg('r','Number of Fine grid BZ kpts   : '//trim(intc(Xk%FGbz%N))) 
   call LIVE_message('[Haydock] Using Double Grid method') 
   call LIVE_message('[Haydock] Number of Coarse grid BZ kpts : '//trim(intc(Xk%nbz))) 
   call LIVE_message('[Haydock] Number of Fine grid BZ kpts   : '//trim(intc(Xk%FGbz%N))) 
 end if
 !
 ! Allocate Haydock vectors for biorthogonal iteration
 !
 if (BSS_uses_DbGd) then
   do i_g=1,BS_nT_grps
   HAYVEC_size = HAYVEC_size + BS_T_grp(i_g)%size * (BS_nkFGbz_in_Tgrp(i_g)-1)
   enddo
 end if
 !
 ! Allocate initial dipole vectors
 allocate(psi_0(BS_nT_grps))
 allocate(phi_0(BS_nT_grps))
 call BS_HAYVEC_alloc(psi_0,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(phi_0,BS_nkFGbz_in_Tgrp)
 !
 if(velocity_correction) then
   allocate(V_tmp(BS_nT_grps))
   call BS_HAYVEC_alloc(V_tmp,BS_nkFGbz_in_Tgrp)
 endif
 !
 if(use_left_and_right) then
   allocate(W0(BS_nT_grps))
   call BS_HAYVEC_alloc(W0,BS_nkFGbz_in_Tgrp)
 endif
 !
 ItParity=1
 it_on_disk=0
 reached_threshold=100._SP
 !
 Cf_size=Haydock_iterMAX+1
 YAMBO_ALLOC(Cf,(Cf_size))
 Cf=cZERO
 !
 ! Allocate biorthogonal Lanczos vectors
 allocate(psi_n(BS_nT_grps))
 allocate(psi_np1(BS_nT_grps))
 allocate(psi_nm1(BS_nT_grps))
 allocate(phi_n(BS_nT_grps))
 allocate(phi_np1(BS_nT_grps))
 allocate(phi_nm1(BS_nT_grps))
 allocate(H_psi(BS_nT_grps))
 allocate(H_phi(BS_nT_grps))
 !
 call BS_HAYVEC_alloc(psi_n,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(psi_nm1,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(psi_np1,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(phi_n,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(phi_nm1,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(phi_np1,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(H_psi,BS_nkFGbz_in_Tgrp)
 call BS_HAYVEC_alloc(H_phi,BS_nkFGbz_in_Tgrp)
 !
 ! Update Descriptions
 !=====================
 call IO_desc_add(BSS_desc,kind="s",str=" ",S=" ")
 !
 call IO_desc_add(BSS_desc,kind="r",str='[Haydock Full Bands] Accuracy (requested)',term='o/o',R=(/Haydock_threshold/))
 call IO_desc_add(BSS_desc,kind="r",str='[Haydock Full Bands] Accuracy (reached)',term='o/o',R=(/reached_threshold/))
 i_desc_acc=BSS_desc%n
 call IO_desc_add(BSS_desc,kind="i",str='[Haydock Full Bands] Iterations',I=(/it_on_disk/))
 i_desc_iter=BSS_desc%n
 !
 call msg('r','Accuracy (requested)      ',Haydock_threshold,"[o/o]")
 !
 ! Initialization of Haydock vectors and tridiagonal matrix elements
 !
 if (it_on_disk==0) then
   !
   do i_g=1,BS_nT_grps
     !
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     !
     psi_n(i_g)%fragment  = cZERO
     psi_nm1(i_g)%fragment= cZERO
     psi_np1(i_g)%fragment= cZERO
     phi_n(i_g)%fragment  = cZERO
     phi_nm1(i_g)%fragment= cZERO
     phi_np1(i_g)%fragment= cZERO
     H_psi(i_g)%fragment  = cZERO
     H_phi(i_g)%fragment  = cZERO
     !
     psi_0(i_g)%fragment  = cZERO
     phi_0(i_g)%fragment  = cZERO
     if(use_left_and_right) W0(i_g)%fragment  = cZERO
     if(velocity_correction) V_tmp(i_g)%fragment = cZERO
     !
   enddo
   !
   Af= 0._SP
   Bf= 0._SP
   !
 end if
 !
 ! Initialize psi_0 with dipole vector
 do i_g=1,BS_nT_grps
   if (.not.PAR_IND_T_groups%element_1D(i_g)) cycle
   do i_ColFG=1,BS_nkFGbz_in_Tgrp(i_g)
     if (trim(OBS)=="optics")      psi_0(i_g)%fragment(:,i_ColFG) = BS_T_grp(i_g)%dipoles_opt(1,:,1)
     if (trim(OBS)=="kerr")        psi_0(i_g)%fragment(:,i_ColFG) = BS_T_grp(i_g)%dipoles_opt(2,:,1)
     if (trim(OBS)=="dichroism")   psi_0(i_g)%fragment(:,i_ColFG) = BS_T_grp(i_g)%dipoles_dic(1,:,1)
     if (trim(OBS)=="magnons")     psi_0(i_g)%fragment(:,i_ColFG) = BS_T_grp(i_g)%dipoles_mag(1,:,1)
     psi_0(i_g)%fragment(:,i_ColFG) = psi_0(i_g)%fragment(:,i_ColFG)*sqrt(cmplx(BS_T_grp(i_g)%f_RES,kind=SP))
   enddo
   !
   if (BSS_uses_DbGd) then
     do i_T=1,BS_T_grp(i_g)%size
       i_k_bz=BS_T_grp(i_g)%table(i_T,1)
       do i_ColFG=1,BS_nkFGbz_in_Tgrp(i_g)
         if (i_ColFG <= Xk%FGbz%N_mapped(i_k_bz)) cycle
         psi_0(i_g)%fragment(i_T,i_ColFG)=cZERO
       enddo
     enddo
   endif
   if(velocity_correction) V_tmp(i_g)%fragment = psi_0(i_g)%fragment/BS_T_grp(i_g)%E
 enddo
 !
 do i_g = 1,BS_nT_grps
   if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
#if defined _MPI
   if(     velocity_correction) call PP_bcast(V_tmp(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
   if(.not.velocity_correction) call PP_bcast(psi_0(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
#endif
 enddo
 !
 if(velocity_correction) call K_multiply_by_V(V_tmp,psi_0,1)
 !
 if(use_left_and_right) then
   !
   do i_g=1,BS_nT_grps
     if (.not.PAR_IND_T_groups%element_1D(i_g)) cycle
     do i_ColFG=1,BS_nkFGbz_in_Tgrp(i_g)
       if (trim(OBS)=="kerr") W0(i_g)%fragment(:,i_ColFG) = BS_T_grp(i_g)%dipoles_opt(1,:,1)
       if (trim(OBS)=="dichroism") W0(i_g)%fragment(:,i_ColFG) = BS_T_grp(i_g)%dipoles_opt(1,:,1)/q0_def_norm
       W0(i_g)%fragment(:,i_ColFG) = W0(i_g)%fragment(:,i_ColFG)*(sqrt(cmplx(BS_T_grp(i_g)%f_RES,kind=SP)))
     enddo
     if(velocity_correction) V_tmp(i_g)%fragment = W0(i_g)%fragment/BS_T_grp(i_g)%E
   enddo
   !
   do i_g = 1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
#if defined _MPI
     if(     velocity_correction) call PP_bcast(V_tmp(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
     if(.not.velocity_correction) call PP_bcast(   W0(i_g)%fragment,0,COMM=PAR_COM_T_Haydock(i_g)%COMM)
#endif
   enddo
   !
   if(velocity_correction) call K_multiply_by_V(V_tmp,W0,1)
   !
 endif
 !
 if(velocity_correction) then
   call BS_HAYVEC_free(V_tmp)
   deallocate(V_tmp)
 endif
 !
 ! Biorthogonal Lanczos Initialization for pseudo-Hermitian system
 !
 if (it_on_disk==0) then
   !
   ! Compute H|psi_0> for initial normalization
   call K_multiply_by_V(psi_0,H_psi,-1)
   !
   ! Compute <psi_0|Theta|H|psi_0> where Theta is the metric operator.
   ! K_dot_product with i_kind=1 applies the Theta metric: factor +1 for resonant,
   ! -1 for anti-resonant blocks. This is the correct behavior for pseudo-Hermitian.
   ! See K_dot_product.F:40 for implementation details.
   call K_dot_product(v0_mod_,psi_0,H_psi,1)
   !
   ! For the coupling case
   if(     l_BS_ares_from_res) Haydock_v0_mod=sqrt(2*real(v0_mod_))
   if(.not.l_BS_ares_from_res) Haydock_v0_mod=sqrt(  real(v0_mod_))
   !
   ! Normalize initial vectors
   do i_g=1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     psi_0(i_g)%fragment   = psi_0(i_g)%fragment/Haydock_v0_mod
     if (use_left_and_right) W0(i_g)%fragment = W0(i_g)%fragment/Haydock_v0_mod
     psi_n(i_g)%fragment   = psi_0(i_g)%fragment
     H_psi(i_g)%fragment = H_psi(i_g)%fragment/Haydock_v0_mod
   enddo
   !
   ! Initialize left vector phi_0 = Theta * psi_0 for pseudo-Hermitian
   ! For coupling case: phi_n = psi_n (with appropriate metric)
   do i_g=1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     phi_n(i_g)%fragment = psi_n(i_g)%fragment
     phi_0(i_g)%fragment = psi_0(i_g)%fragment
   enddo
   !
   ! Compute initial coefficient: Cf(1) = <phi_0|Theta|psi_0>
   ! K_dot_product with i_kind=1 applies the metric operator correctly.
   if (.not.use_left_and_right) call K_dot_product(v0_mod_,phi_0,psi_0,1)
   if (     use_left_and_right) call K_dot_product(v0_mod_,W0,psi_0,1)
   !
   if(     l_BS_ares_from_res) Cf(1)=2*cI*aimag(v0_mod_)
   if(.not.l_BS_ares_from_res) Cf(1)=           v0_mod_
   !
 endif
 !
 ! [***] Biorthogonal Lanczos iteration loop
 !
 do it=it_on_disk+1,Haydock_iterMAX
   !
   ! Biorthogonal Lanczos for pseudo-Hermitian coupling Hamiltonian
   !
   ! Step 1: Compute A(n) = <phi_n|H|psi_n>
   !
   if (l_BS_ares_from_res) then
     ! A(n) = 0 by symmetry for the coupling case with symmetric R/A
     Af(it)=0._SP
     ItParity=(-1)**(it+1)
   else
     call K_dot_product(v0_mod_,H_psi,H_psi,1)
     Af(it)=real(v0_mod_,SP)
   endif
   !
   ! Step 2: Orthogonalize
   ! |psi_tilde> = H|psi_n> - A(n)|psi_n> - B(n)|psi_{n-1}>
   ! |phi_tilde> = H^dag|phi_n> - A(n)^*|phi_n> - B(n)^*|phi_{n-1}>
   !
   if (.not.l_BS_ares_from_res) then
     do i_g = 1,BS_nT_grps
       if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
       do i_ColFG=1,BS_nkFGbz_in_Tgrp(i_g)
         call V_plus_alpha_V_omp(BS_T_grp(i_g)%size,-Af(it)*cONE,psi_n(i_g)%fragment(:,i_ColFG),H_psi(i_g)%fragment(:,i_ColFG))
       enddo 
     enddo
   endif
   !
   do i_g = 1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     do i_ColFG=1,BS_nkFGbz_in_Tgrp(i_g)
       call V_plus_alpha_V_omp(BS_T_grp(i_g)%size,-Bf(it)*cONE,psi_nm1(i_g)%fragment(:,i_ColFG),H_psi(i_g)%fragment(:,i_ColFG))
     enddo
   end do
   !
   ! Update vectors: |psi_{n-1}> = |psi_n>, |psi_n> = |psi_tilde>
   !
   do i_g=1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     psi_nm1(i_g)%fragment = psi_n(i_g)%fragment
     phi_nm1(i_g)%fragment = phi_n(i_g)%fragment
     psi_n(i_g)%fragment   = H_psi(i_g)%fragment
   enddo
   !
   ! Step 3: Compute H|psi_n> for next iteration
   call K_multiply_by_V(psi_n,H_psi,ItParity)
   !
   ! Step 4: Compute B(n+1) = sqrt(<phi_n|Theta|H|psi_n>)
   ! K_dot_product with i_kind=1 correctly applies the Theta metric operator.
   !
   call K_dot_product(b_squared,psi_n,H_psi,1)
   !
   ! Check for numerical stability
   if (real(b_squared) <= 0._SP) then
     call msg('r','[Haydock Full Bands] Lanczos breakdown detected at iteration '//trim(intc(it)))
     call msg('r','[Haydock Full Bands] B^2 = ',real(b_squared))
     exit
   endif
   !
   ! Check imaginary part tolerance
   if (abs(aimag(b_squared)) > IMAG_TOL * abs(real(b_squared))) then
     call warning('[Haydock Full Bands] Large imaginary part in B^2 at iteration '//trim(intc(it)))
   endif
   !
   if (     l_BS_ares_from_res) Bf(it+1)=sqrt(2*real(b_squared,SP))
   if (.not.l_BS_ares_from_res) Bf(it+1)=sqrt(  real(b_squared,SP))
   !
   ! Check for INF or NaN
   if( Bf(it+1)==(Bf(it+1)+1._SP) ) call error('Bf=INF likely because some eigenvalue of the BSE is negative.')
   if (NAN(Bf(it+1)))  call error(' Bf=NaN likely because some eigenvalue of BSE is negative.')
   !
   ! Check for breakdown
   if (Bf(it+1) < B_MIN) then
     call msg('r','[Haydock Full Bands] Lanczos convergence achieved (B < B_MIN) at iteration '//trim(intc(it)))
     exit
   endif
   !
   ! Step 5: Normalize
   ! |psi_{n+1}> = |psi_tilde> / B(n+1)
   ! |phi_{n+1}> = |phi_tilde> / B(n+1)^*
   !
   do i_g=1,BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     H_psi(i_g)%fragment  = H_psi(i_g)%fragment/Bf(it+1)
     psi_n(i_g)%fragment  = psi_n(i_g)%fragment/Bf(it+1)
     phi_n(i_g)%fragment  = psi_n(i_g)%fragment  ! For pseudo-Hermitian with Theta metric
   enddo
   !
   ! Step 6: Compute Cf(n+1) = <phi_0|Theta|psi_n>
   ! K_dot_product with i_kind=1 correctly applies the Theta metric operator.
   !
   Cf_size=it+1
   !
   if (.not.use_left_and_right) call K_dot_product(v0_mod_,phi_0,psi_n,1)
   if (     use_left_and_right) call K_dot_product(v0_mod_,W0,psi_n,1)
   !
   if (l_BS_ares_from_res) then
     ! Using symmetry arguments
     Cf(it+1)=   ((1._SP+real(ItParity,SP))*real(v0_mod_,SP)) +&
     &        cI*((1._SP-real(ItParity,SP))*aimag(v0_mod_))
   else
     Cf(it+1)=   v0_mod_
   endif
   !
   ! Compute response and check convergence
   !
   if (it>2) then
     !
     call K_Haydock_response(iq,it,Cf_size-1,OBS,"pseudo-herm",W,real(Af(:it)),Bf(:it+1),Cf(:Cf_size-1),reached_threshold,FG_factor)
     !
     BSS_desc%rval(1,i_desc_acc) =reached_threshold
     BSS_desc%ival(1,i_desc_iter)=it
     !
     dumb_ch='[Haydock Full Bands] Iteration '//trim(intc(it))//' Accuracy '//&
&            trim(real2ch(reached_threshold))//'|'//trim(real2ch(Haydock_threshold))
     if (tty_size>0.and.(it_on_disk==0.or. it_on_disk< it-1)) call LIVE_message(dumb_ch,"r","","%s")
     if (tty_size>0.and.(it_on_disk> 0.and.it_on_disk==it-1)) call LIVE_message(dumb_ch,"n","","%s")
     if (tty_size<0) call LIVE_message(dumb_ch,"n","","%s")
     !
   else
     !
     dumb_ch='[Haydock Full Bands] Iteration '//trim(intc(it))
     call LIVE_message(dumb_ch,"n","","%s")
     !
   endif
   !
   tmp_converged=(reached_threshold<abs(Haydock_threshold)) .or. STOP_now(.FALSE.)
   !
   if(     tmp_converged) Haydock_n_converged=Haydock_n_converged+1
   if(.not.tmp_converged) Haydock_n_converged=0
   !
   if ( mod(it,max(1,Haydock_iterIO/10))==0) call K_OUTPUT(iq,W,"override dump Haydock close",OBS)
   !
   if (Haydock_n_converged==Haydock_required_converged) then
     !
     call K_OUTPUT(iq,W,"override dump Haydock close ",OBS)
     !
     call msg('r','Accuracy (reached)        ',reached_threshold,"[o/o]")
     call msg('r','Iterations                ',it)
     exit
     !
   endif
   !
 enddo
 !
 ! CLEAN
 !=======
 !
 call BS_HAYVEC_free(psi_n)
 call BS_HAYVEC_free(psi_np1)
 call BS_HAYVEC_free(psi_nm1)
 call BS_HAYVEC_free(phi_n)
 call BS_HAYVEC_free(phi_np1)
 call BS_HAYVEC_free(phi_nm1)
 call BS_HAYVEC_free(H_psi)
 call BS_HAYVEC_free(H_phi)
 deallocate(psi_n,psi_np1,psi_nm1)
 deallocate(phi_n,phi_np1,phi_nm1)
 deallocate(H_psi,H_phi)
 !
 call BS_HAYVEC_free(psi_0)
 call BS_HAYVEC_free(phi_0)
 deallocate(psi_0,phi_0)
 !
 if (use_left_and_right)  then
   call BS_HAYVEC_free(W0)
   deallocate(W0)
 endif
 !
 YAMBO_FREE(Cf)
 !
 call PARALLEL_Haydock_VEC_COMMs('reset')
 !
 call timing('Haydock Full Bands Solver',OPR='stop')
 !
end subroutine K_Haydock_full_bands
