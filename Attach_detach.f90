module Attach_detach
  
  use CommonParams
  use Rnd_seed
  ! use Minima_switch
  use Actomyos_Energy_Params, only : d_min, n_m
  implicit none
  
  real(8), parameter :: Q10 = 1.5**((temp_cels-37.d0)/10)
  
  !for  attachment_detachment

  real(8), parameter :: gamma_coop = 40.0 !gamma value for the cooperativity between myosin heads
  real(8), parameter :: MS_fact=1.9
  real(8), parameter :: WtoS_fact = 2.95d0
  real(8), parameter :: dopeWD = 4.d-4*Q10
  real(8), parameter :: geom_fac_d = 1.d0
  real(8), parameter :: Xtr = 10.d0 !nm distance threshold for the attachment rate 
  real(8), parameter :: dope = 0.125d0*Q10 !term to fast the process (used to met the time scaling simulation/experiment if needed)
  real(8), parameter :: kwd = MS_fact*dopeWD*(60e0*(gamma_coop**2))*(1.d-9) !ns-1 weak to detached constant rate
  real(8), parameter :: kdw = dopeWD * (85.d0) * (10.d-9) !ns-1 detached to weak constant rate
  real(8), parameter :: probs = WtoS_fact*4.5d0*(10.d-9)*dope !ns-1 nm-1 attachment rate slope
  real(8) :: probd = geom_fac_d*(1/WtoS_fact)*probs*0.27d0 !ns-plot ++.eq.ind_coord = 1 nm-1 detachment rate slope
  real(8), parameter :: probd_N = 490.0*(10.d-9) ! ns-1 high detachment rate in negative region
  real(8) ::  csi_R = 0.0, csi_L = 0.0
  real(8), parameter :: ON_max_old = 42*MS_fact*120.d0
  real(8), parameter :: T_on_max = 500.d0, T_on_min = 33.5d0 !60% and 10% of 300 kPa
 
 
  
  private

  public probd, probs
  
  public :: attachment_detachment
 

contains

  subroutine attachment_detachment()
    implicit none
    
    integer :: i, j, n_coop, ind_TT, k, min_temp(1), XB_fil
    real(8) :: y, rnd, rnd_s, stod, Ca_mu,x_SL, k_for, k_bac, MI_AD_fact, MI_AT_fact, MI_WD_fact, T_fil, ON_fact
    real(8) :: F_titin, pos_fact
    real(8), parameter :: alpha_AD_fact = 4.d0, alpha_AT_fact = 6.0d0, alpha_WD_fact = 1.d0
    real(8), parameter :: a_titin=9.d0, b_titin = 2*a_titin/(150**2)

    ON_min = piperine*MS_fact*20.d0
    ON_max = (ON_max_old)*T_on_max/1200.d0 !to adapt the new T_on_man to a more elegant limit
  

    if (z .le. -50.d0) then
       F_titin = 0.d0
    else
       F_titin = 2*0.5*b_titin*(z+50)**2
    end if
    F_titin = 0*titin_tension
     
    !$omp parallel do default(shared) &
    !$omp private(j, i, y, k, rnd, rnd_s, stod, n_coop, x_SL, csi_R, csi_L, Ca_mu, ind_TT, k_for, k_bac, min_temp, &
    MI_AD_fact, MI_AT_fact, MI_WD_fact, T_fil, XB_fil, ON_fact, pos_fact)&
         !$omp reduction(+: ATP_cons)
    

    do j = 1, Nfil

      

       !Xb based cooperativity
    !   if (XB_fil .le. 2) then
    !      ON_fact = ON_min
    !   else if (XB_fil .le. 8) then
    !      ON_fact = ON_max-(ON_max-ON_min)/(8-2)*(8-XB_fil)
    !   else if (XB_fil .gt. 8) then
    !      ON_fact=ON_max
    !   end if

       !Tension based cooperativity
       
    !   if (T_fil .gt.0) then
    !      ON_fact = ON_min + (ON_max-ON_min)*exp(-ON_tau/(max(0.d0,(T_fil+F_titin))**2))
    !   else
    !      ON_fact=ON_min
    !   end if
  
    T_fil=0.d0

    
       
       do i = Nxb,1,-1
          if (1 .eq. 1) then
             MI_AD_fact = alpha_AD_fact
             MI_AT_fact = alpha_AT_fact
             MI_WD_fact = alpha_WD_fact
          else
             MI_AD_fact = 1.d0
             MI_AT_FACT = 1.d0
             MI_WD_fact = 1.d0
          end if

          !analysis of the tension
          
          if (AttachFlag(i,j) .ge. a_state_1) then
             if ((x(i, j) + y0(i, j)) .ge. 0) then
                T_fil = T_fil + Kplus*(x(i,j) + y0(i, j))
             else
                T_fil = T_fil + Kminus*(x(i,j) + y0(i, j))
             end if
          end if

          if (i .le. Nxb/2) then 
             pos_fact=0.!ON_min/3.
          else
             pos_fact=2.d0
          end if

          pos_fact = 1.0 !one is for an uniform distribution of the ON
          
          !analysis of the ON_fact
          ON_fact = 0.d0
          if (T_fil .ge. T_on_max) then 
             ON_fact = ON_max
          else if (T_fil+F_titin .gt. T_on_min) then
             ON_fact = (ON_max-pos_fact*ON_min)/(T_on_max-T_on_min)*(T_fil-T_on_min)+pos_fact*ON_min
          else if (T_fil .le. T_on_min) then
             ON_fact = pos_fact*ON_min
          end if



!!!str Kon(T)
!          if (mod(TotalTime,10000.) .eq. 0.) then
!             if ((j==900) .and. (myId==0) .and. (i==1)) then
!                write(3,'(6e14.5)') T_fil, ON_fact,  alpha_WD_fact*kdw*(ON_fact)*DT*1e6
!             end if
!          end if
!          
!          if (mod(TotalTime,10000.) .eq. 0.) then
!             if ((j==900) .and. (myId==2) .and. (i==1)) then
!                write(4,'(6e14.5)') T_fil, ON_fact,  alpha_WD_fact*kdw*(ON_fact)*DT*1e6
!             end if
!          end if          
!!!!end kON(T)
          
          
          n_coop=0 !co-operativity n factor
          x_SL = x(i,j)+LB/2.d0+i*XB_s !relative position of MH on AF
          if (i .gt. 1) then
             if (AttachFlag(i-1, j) .ne. d_state) n_coop = n_coop + 1
          end if
          if (i .lt. Nxb) then
             if ( AttachFlag(i+1, j) .ne. d_state) n_coop = n_coop + 1
          end if
          rnd = unifrd__(rndSeed(i, j))
         
          if  (AttachFlag(i, j) .eq. d_state) then 
             csi_R = SL_R(x_SL) !probability factor for SL dependence (right)
             csi_L = SL_L(x_SL)!probability factor for SL dependence (overlapping)
             if (rnd<MI_WD_fact*kdw*(ON_fact)*csi_R*csi_L*DT) then
                AttachFlag(i, j) = w_state 
             end if
          else if (AttachFlag(i, j) .eq. w_state) then
             Ca_mu = mu_Ca
             ind_TT = int((x_SL-z)/d_TT)+1 
             if (ind_TT <= n_tt .and. ind_TT >= 1) then
                if (Ca_Flag(ind_TT, j) .eq. Ca_on) Ca_mu = 1.d0 
             end if
             csi_R = SL_R(x_SL) !probability factor for SL dependence (right)
             csi_L = SL_L(x_SL) !probability factor for SL dependence (overlapping)
             if (rnd<MI_WD_fact*kwd*DT) then
                AttachFlag(i, j)=d_state
             else if  (x(i, j)>0.d0 .and. x(i,j)<Xtr) then
                if (rnd<(MI_WD_fact*kwd+x(i, j)*probs*MI_AT_fact*csi_R*csi_L/Ca_mu)*DT) then
                   AttachFlag(i, j) = a_state_1
                   xAttached(i, j) = x(i, j)
                   zAttached(i, j) = z 
                   y0(i, j) = y0_min+(y0_max-y0_min)*unifrd__(rndSeed_y0(i,j))
                   ind_TT = int((x_SL-z)/d_TT)+1 
                   if (ind_TT < 0 .or. ind_TT > n_tt) then
                      write(6, *) "Error ind_TT=", ind_TT, i, csi_R, csi_L, z, x(i,j), X_SL
                 !     stop
                   end if
                   TT_b(ind_TT, j) = TT_b(ind_TT, j) + 1
                   !  write(10,'(e14.5)') xAttached(i, j)
                end if
             end if
          else if (AttachFlag(i, j) .ge. a_state_1) then
           
             ind_TT = int((x_SL-z)/d_TT)+1 
             !  cycle
             if ((x(i,j)+y0(i,j)) .ge. 0) then
                stod = (x(i,j)+y0(i,j))*probd*MI_AD_fact
            !    stod = 4*probd
             else 
                stod =  probd_N*(1+abs(x(i, j)+y0(i,j)))
             end if
             if (AttachFlag(i,j) .eq. a_state_3) then
                stod=stod*1.d0
             end if
             if (abs(x(i,j)+y0(i,j)) .ge. 15.d0) then
                stod = 1.d0
             end if
             y=-(xAttached(i,j)+ (z - zAttached(i, j)) + y0(i, j))
             if (-y .le. delta_min) then
                stod = 1.d0
             end if

             if (rnd<stod*DT) then
                AttachFlag(i, j) = w_state
                ATP_cons = ATP_cons + 1
                x_SL = xAttached(i, j) + LB/2.d0 + i*XB_s
                ind_TT = int((x_SL-zAttached(i,j))/d_TT) + 1
                TT_b(ind_TT, j) = TT_b(ind_TT, j) - 1
                if (TT_b(ind_TT, j) .lt. 0) then
                   write(6,*) "TT counter fails!!!",myId,  i, z, ind_TT
                   stop
                end if
             else
                
                if ((-y .gt. delta_max) .or. (-y .lt. delta_min)) then
                   write(*,*) "Error exceeding minima" , xAttached(i,j), z, zAttached(i, j), y0(i, j), i, j
                   
                   stop
                end if
                
                k = min(max(1,NINT(1+(y-delta_min)/delta_int)),int((delta_max-delta_min)/delta_int))
             
                
                if (AttachFlag(i,j) .eq. a_state_1) then
                   k_for = const_rate(k,2)
                else if (AttachFlag(i,j) .eq. a_state_2) then
                   k_for = const_rate(k,4)
                   k_bac = const_rate(k,3)
                else if (AttachFlag(i,j) .eq. a_state_3) then
                   k_bac = const_rate(k,5)
                end if
                 
                if (AttachFlag(i,j) .eq. a_state_1) then
                   if (rnd<(stod*DT+k_for*DT)) then
                      AttachFlag(i,j)=a_state_2
                   end if
                else if (AttachFlag(i,j) .eq. a_state_2) then
                   if (rnd<(stod*DT+k_for*DT)) then
                      AttachFlag(i,j)=a_state_3
                   else if (rnd<(stod+k_for+k_bac)*DT) then
                      AttachFlag(i,j)=a_state_1
                   end if
                else if  (AttachFlag(i,j) .eq. a_state_3) then
                   if (rnd<(stod*DT+k_bac*DT)) then
                      AttachFlag(i,j)=a_state_2
                   end if
                end if
             end if
          end if
       end do
    end do
       
  end subroutine attachment_detachment
  
  
  real(8) function SL_R(x_SL)
    implicit none
    real(8), intent(in) :: x_SL
    if (x_SL .le. z) then
       SL_R = 0
    else if (X_SL .le. z+2) then
       SL_R = dexp(-(((z+2)-x_SL)**2)/(a_R**2))
    else if ((x_SL .gt. z+2) .and. (x_SL .lt. z+LA-2)) then
       SL_R = 1
    else if ((x_SL .ge. z+LA-2) .and. (x_SL .le. z+LA)) then
       SL_R =  dexp(-((x_SL-(z+LA-2))**2)/(a_R**2))
    else if (X_SL .ge. z+LA) then
       SL_R = 0
    end if
  end function SL_R
  
  
  real(8) function SL_L(x_SL)
    implicit none
    real(8), intent(in) :: x_SL
    if (x_SL .le. -z) then
       SL_L = dexp(-((-z-x_SL)**2)/(a_L**2))
    else
       SL_L = 1
    end if
  end function SL_L

  
end module Attach_detach
