
module CommonParams
  implicit none

  real(8), parameter :: temp_cels = 25.d0
  real(8), parameter :: KB_T  = 0.0138*(temp_cels+273.15d0)!4.27d0 ! (pN*nm) (4 celsius degrees-> kT=3.82; 24 cd -> kT=4.06; 27 c.d.-> kT=4.14 pN nm)
  real(8), parameter :: pi = 4.d0*datan(1.d0)


  real(8), parameter :: LB = 100.d0, LA = 1200.d0, LM = 1650.d0, a_R = 10.d0, a_L = 10.d0 !(nm) parameter for the SL dependence of A/D rates
  real(8), parameter :: d_TT = 36 ! distance between two consecutive TT units
  integer, parameter ::  NU_MF = 147 !number of myosin heads for half myosin filament
  real(8), parameter  :: XB_s = ((LM-LB)/2.d0)/NU_MF
  integer, parameter :: Nxb = NU_MF
  integer, parameter :: Nfil = 1800
  integer, parameter ::  w_state = 1, d_state = 0  !state of the strongly or weakly attached and detached state
  integer, parameter ::  a_state_1 =2, a_state_2 = 3, a_state_3 = 4
  integer, parameter :: n_TT = LA/d_tt
  real(8), parameter :: y0_min = -2.25d0, y0_max = 2.25d0
  real(8), parameter :: mu_Ca=555.d0 !inibition parameter for no Ca in TT unit

  real(8), parameter :: Kplus = 1.7d0 ! myosin arm spring force constant for positive (pN/nm)
  real(8), parameter :: Kminus = 1.7d0 ! myosin arm spring force constant for minus (pN/nm)

  
  real(8), parameter :: DT = 1.d3  !nano-second

  
  real(8), parameter :: steady_time = 2.d8!, start_time = steady_time +0.507d8  !fo the analysis of Fusi, + 110 ms +st (200ms)

   real(8) :: piperine !Effect of the piperine
   real(8) :: start_time !for the analysis of Fusi, + 110 ms +st (200ms)
    real(8) :: ON_min = 0.d0
  real(8) :: ON_max = 0.d0!to adapt the new T_on_man to a more elegant limit

  
  integer, parameter :: Ca_off = 0, Ca_on_x = 1, Ca_on = 2  !state of the TT unit

  real(8) :: TotalTimeMax
  integer :: AttachFlag(Nxb, Nfil), rndSeed(Nxb, Nfil), rndSeedT(n_TT, Nfil), rndSeed_s(Nxb,nFil), rndSeed_y0(Nxb, NFil)
  real(8) :: x(Nxb, Nfil), dxdt(Nxb, Nfil), y0(Nxb, Nfil)
  real(8) :: xAttached(Nxb, Nfil), zAttached(Nxb, Nfil) 
  integer :: Ca_Flag(n_TT, Nfil), TT_b(n_TT, Nfil)
  integer ::  rndSeedz
  real(8) :: z, dzdt, z_clamp, dz_clampdt, z_clamp_ini=0
  real(8) :: Ktrap ! (pN/real)
 

 real(8), parameter :: etaMyosinHead = 70.d0 ! Viscosity for Myosin Head (pN*ns/nm)
  real(8), parameter :: etaActinFilament = etaMyosinHead * Nxb * Nfil! Viscosity for Actin Filament (pN*ns/nm)
  
  real(8) :: F_appl = 0.d0, Total_Force = 0.d0, F_ext_save = 0.d0, titin_tension
  real(8) :: TotalTime !nano-second
  integer :: ATP_cons = 0
 
  integer :: set_up = 3
  real(8) :: perc = 0.1d0, Tjump=0.d0,  delta_step = 0.d0, Tjump_Ca1=1e6,  Tjump_Ca2=1e6

  logical :: isometric_ini = .FALSE., rate_call = .FALSE.
 

  !MPI
  integer :: myId, nProcs, ierr
  character(len=2) fid

  real(8) :: Ca_conc = 1.d-3, Ca_shift_in




!!!!!!!!Rate reading parameters

 integer :: len_det, len_rate
 integer,dimension(2) :: len_at_1, len_at_2, len_at_3
 integer, parameter :: L=10000
 real(8), dimension(L) :: detach_dist


 real(8) :: k_for, k_bac, y_rate_min, y_rate_max, dy_rate, y_att_min, y_att_max, dy_att


 real(8), parameter :: delta_min = -80, delta_max  = 20, delta_int = 0.1
 integer, parameter :: Stretch_ind = int((delta_max-delta_min)/delta_int)
 real(8), dimension(Stretch_ind,L) ::  attach_fst_dist, attach_snd_dist, attach_trd_dist
 real(8) ::  const_rate(Stretch_ind,5), detach_dist_t(L)
 
 real(8) :: V_mean, ATP_mean
 
 
end module CommonParams
