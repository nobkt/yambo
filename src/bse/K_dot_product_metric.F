!
! License-Identifier: GPL
!
! Copyright (C) 2024 The Yambo Team
!
! Authors: Memory-efficient full-band BSE implementation
!
! This subroutine computes the metric dot product for pseudo-Hermitian systems.
! 
! For the coupling BSE Hamiltonian with structure:
!   H_BSE = [ H^res      C     ]
!           [ -C*    -H^res*   ]
!
! The metric operator is:
!   Theta = [ I    0  ]
!           [ 0   -I  ]
!
! The metric inner product is:
!   <V|Theta|W> = <V_+|W_+> - <V_-|W_->
!
! where V = (V_+, V_-) and W = (W_+, W_-)
!
! Reference: mydoc/bse_full_bands_coupling_specification.md Section 5.1
!
subroutine K_dot_product_metric(dot,V,W,i_kind,l_use_metric)
 !
 ! Metric dot product between distributed vectors for pseudo-Hermitian BSE
 !
 ! i_kind==1 --> Metric dot product for pseudo-Hermitian: <V|Theta|W>
 ! i_kind==0 --> Standard dot product: <V|W>
 !
 ! l_use_metric==.TRUE.  --> Apply metric operator Theta
 ! l_use_metric==.FALSE. --> Standard dot product (same as K_dot_product)
 !
 use pars,           ONLY:SP,cZERO
 use BS_solvers,     ONLY:Haydock_vector
 use BS,             ONLY:BS_nT_grps,BS_T_grp,BS_K_dim,BS_H_dim,l_BS_ares_from_res,BS_n_eh_spaces
 use parallel_m,     ONLY:PAR_IND_T_groups,PAR_COM_Xk_ibz_INDEX,PAR_COM_eh_INDEX
 use parallel_int,   ONLY:PP_redux_wait
 use wrapper_omp,    ONLY:Vstar_dot_V_omp
 !
 implicit none
 !
 integer,              intent(in) :: i_kind
 logical,              intent(in) :: l_use_metric
 type(Haydock_vector), intent(in) :: V(BS_nT_grps),W(BS_nT_grps)
 complex(SP),          intent(out):: dot
 !
 ! Workspace
 !
 integer       :: i_T_g,i_factor,i_ColFG,ColsFG_T
 integer       :: i_res_ares
 !
 dot = cZERO
 ! 
 ! For the coupling case with full matrix (BS_n_eh_spaces==2),
 ! we have separate resonant and anti-resonant transition groups.
 !
 ! The metric operator Theta = diag(I, -I) gives:
 ! - Factor +1 for resonant part (i_res_ares==1)
 ! - Factor -1 for anti-resonant part (i_res_ares==2)
 !
 ! When l_use_metric==.TRUE., we apply this sign factor.
 !
 do i_T_g=1,BS_nT_grps
   !
   if (.not.PAR_IND_T_groups%element_1D(i_T_g)) cycle
   !
   ! Get the resonant/anti-resonant index for this transition group
   i_res_ares = BS_T_grp(i_T_g)%i_res_ares
   !
   ! Compute the factor based on i_kind and metric usage
   !
   ! Standard factor from K_dot_product (handles the coupling symmetry)
   i_factor = ( (-1)**(i_res_ares-1) )**i_kind
   !
   ! If using metric for pseudo-Hermitian, apply additional sign for anti-resonant
   if (l_use_metric .and. i_res_ares == 2) then
     i_factor = -i_factor
   endif
   !
   ColsFG_T = size(V(i_T_g)%fragment,2)
   do i_ColFG=1,ColsFG_T
     dot=dot+ i_factor*Vstar_dot_V_omp(BS_T_grp(i_T_g)%size,V(i_T_g)%fragment(:,i_ColFG),W(i_T_g)%fragment(:,i_ColFG))
   enddo
   !
 enddo
 !
 ! Reduce all
 !
 call PP_redux_wait( dot ,COMM=PAR_COM_Xk_ibz_INDEX%COMM )
 call PP_redux_wait( dot ,COMM=PAR_COM_eh_INDEX%COMM )
 !
end subroutine K_dot_product_metric
!
!
!
!===========================================================================
! Compute biorthogonality error for convergence monitoring
!
! This checks how well the biorthogonality condition <phi_m|psi_n> = delta_mn
! is satisfied. Large deviations indicate need for reorthogonalization.
!
! Reference: mydoc/bse_full_bands_coupling_specification.md Section 7.1
!===========================================================================
subroutine K_check_biorthogonality(phi,psi,n_vecs,max_error)
 !
 use pars,           ONLY:SP,cZERO
 use BS_solvers,     ONLY:Haydock_vector
 use BS,             ONLY:BS_nT_grps
 !
 implicit none
 !
 integer,              intent(in)  :: n_vecs
 type(Haydock_vector), intent(in)  :: phi(BS_nT_grps,n_vecs)
 type(Haydock_vector), intent(in)  :: psi(BS_nT_grps,n_vecs)
 real(SP),             intent(out) :: max_error
 !
 ! Workspace
 !
 integer       :: i_m, i_n
 complex(SP)   :: dot
 real(SP)      :: error
 !
 max_error = 0._SP
 !
 ! Check off-diagonal elements: should be zero
 do i_m = 1, n_vecs
   do i_n = 1, n_vecs
     !
     call K_dot_product_metric(dot, phi(:,i_m), psi(:,i_n), 1, .TRUE.)
     !
     if (i_m == i_n) then
       ! Diagonal: should be 1
       error = abs(dot - (1._SP, 0._SP))
     else
       ! Off-diagonal: should be 0
       error = abs(dot)
     endif
     !
     max_error = max(max_error, error)
     !
   enddo
 enddo
 !
end subroutine K_check_biorthogonality
!
!
!
!===========================================================================
! Reorthogonalize biorthogonal Lanczos vectors when numerical drift occurs
!
! Uses modified Gram-Schmidt orthogonalization to restore biorthogonality:
! |psi_new> = |psi> - sum_j <phi_j|psi> |psi_j>
! |phi_new> = |phi> - sum_j <phi|psi_j>* |phi_j>
!
! Reference: mydoc/bse_full_bands_coupling_specification.md Section 7.1
!===========================================================================
subroutine K_reorthogonalize_biorthogonal(psi_new, phi_new, psi_prev, phi_prev, n_prev)
 !
 use pars,           ONLY:SP,cZERO,cONE
 use BS_solvers,     ONLY:Haydock_vector
 use BS,             ONLY:BS_nT_grps,BS_T_grp
 use parallel_m,     ONLY:PAR_IND_T_Haydock
 use wrapper_omp,    ONLY:V_plus_alpha_V_omp
 !
 implicit none
 !
 integer,              intent(in)    :: n_prev
 type(Haydock_vector), intent(inout) :: psi_new(BS_nT_grps)
 type(Haydock_vector), intent(inout) :: phi_new(BS_nT_grps)
 type(Haydock_vector), intent(in)    :: psi_prev(BS_nT_grps,n_prev)
 type(Haydock_vector), intent(in)    :: phi_prev(BS_nT_grps,n_prev)
 !
 ! Workspace
 !
 integer       :: j, i_g, i_ColFG
 complex(SP)   :: coeff_psi, coeff_phi
 !
 ! Modified Gram-Schmidt reorthogonalization
 !
 do j = 1, n_prev
   !
   ! Compute projection coefficients
   ! coeff_psi = <phi_j|psi_new>
   ! coeff_phi = <phi_new|psi_j>^*
   !
   call K_dot_product_metric(coeff_psi, phi_prev(:,j), psi_new, 1, .TRUE.)
   call K_dot_product_metric(coeff_phi, phi_new, psi_prev(:,j), 1, .TRUE.)
   !
   ! Subtract projections
   ! |psi_new> = |psi_new> - coeff_psi * |psi_j>
   ! |phi_new> = |phi_new> - coeff_phi^* * |phi_j>
   !
   do i_g = 1, BS_nT_grps
     if (.not.PAR_IND_T_Haydock%element_1D(i_g)) cycle
     do i_ColFG = 1, size(psi_new(i_g)%fragment, 2)
       call V_plus_alpha_V_omp(BS_T_grp(i_g)%size, -coeff_psi, &
                               psi_prev(i_g,j)%fragment(:,i_ColFG), &
                               psi_new(i_g)%fragment(:,i_ColFG))
       call V_plus_alpha_V_omp(BS_T_grp(i_g)%size, -conjg(coeff_phi), &
                               phi_prev(i_g,j)%fragment(:,i_ColFG), &
                               phi_new(i_g)%fragment(:,i_ColFG))
     enddo
   enddo
   !
 enddo
 !
end subroutine K_reorthogonalize_biorthogonal
