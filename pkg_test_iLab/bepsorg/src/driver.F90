

!******************************************************
!! main program for BEPS
!! Editted : J.Wang
!! Date    : June2017
!******************************************************

program main
use shr_kind_mod, only: r8 =>shr_kind_r8
use controlInput_mod
use beps_time_manager
use bepstype
use bepstypeInit
use ecoRespMod
use mid_results
!use mpi_mod
use beps_soilMod
use beps_par
use outputMod
use restart
use esmf
implicit none

type(climatedata)    :: meteo
type(results)        :: mid_res
real(r8)             :: CosZs,hr_loc,hr_arc             !! solar zenith angle, local time, local time arc
integer              :: i,j,k,llll,jj,ii,kk
type(soil)           :: soilp                !! at single point
real(r8)             :: Ccd(0:4),Cssd(0:4),Csmd(0:4),Cfsd(0:4),Cfmd(0:4),Csm(0:4),Cm(0:4),Cs(0:4),Cp(0:4)
real(r8)             :: param(0:49),var_o(0:40),var_n(0:40),coef(0:49)
real(r8)             :: lai
type(para),pointer   :: ppar                  !! initial parameters
type(surf),pointer   :: bfields               !! boundary fields
type(forc),pointer   :: climate               !! climate fieldsii =
type(soils),pointer  :: psoil                 !! GLobally
type(res),pointer    :: pp                    !! for output

real(r8)             :: ratio_cloud,shortRad_df,shortRad_dir
!real(r8)             :: NPP_yr_acc(npoints,PFT)            !! for storing yearly accumulated NPP, Mh/ha/s
! real(r8)             :: agb2vod
! real(r8)             :: D0(1:9)
! real(r8)             :: taweff(1:9)
integer              :: kount,rst_nstep
integer              :: yr,mn,dy,tod,caldy,n_meteo
!-- iLab::converted from real(r8) to integer
integer              :: doys
integer              :: ierr
real(r8)             :: daylen
!-- iLab::for revised interface to 'av_output'
character(len=len('YYYY-MM-DDTHH:MM:SS')), save :: ref_date = ''
integer :: yr_ref, mn_ref, dy_ref, tod_ref
real(r8)             :: secs_elapsed,secs_meteo
! .. Parameters for model-internal use
real(r8)             :: pio180 = PI/180.
! variables for daily input, used in climin and getmonth
real(r8)             :: spds, cpds
! .. Local Arrays
real(r8)             :: atmean, atrange
! .. Local Scalars ..^M
real(r8)             :: r
real(r8)             :: rdaymid, delta, arg, h0, h1, sd, sd1, dhour, tmin, tmp1
real(r8)             :: a, b, sunset_arc
integer              :: nd                   !! for counting the time-step number, i.e. ith step
! .. Intrinsic Functions ..
intrinsic ACOS, COS, SIN, MOD,atan, REAL,int
! parameters used for calculating VOD,@MOUSONG.WU,2019-11
!NPP_yr_acc(:,:) = 0.
!agb2vod = 0.9517
!D0 = (/0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05/)
!taweff = (/0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006/)

!! setting up MPI enviroments with the namelists
!   call Initmpi()

call rdnamelist()
if(nscale == 1) then
  nlp = n_site
  npoints = nlp
  write(*,*) 'site points check', npoints
end if

!! Initialize the beps types
   call Initbepstype()

!! Initialize output
   call Init_output

!--------------------------------------------------------------------------
! Initialize ESMF.  This is done outside of the ESMF_INTERFACE ifdef
! because it is needed for the time manager, even if the ESMF_INTERFACE
! is not used.
!--------------------------------------------------------------------------
call ESMF_Initialize()

!! setting time manager
   if(nsrest == nsrStartup) then
   !!! calling time_manager set init
     call set_timemgr_init(calendar_in  = calendar,&
                         start_ymd_in = icdate  ,&
                         start_tod_in = icsec   ,&
                         nelapse_in   = sim_duration,&
                         dtime_in     = step)

   else if(nsrest == nsrContinue) then
      call restart_io("read")

      call set_timemgr_init(calendar_in  = calendar,   &
                            start_ymd_in = rst_icdate, &
                            start_tod_in = rst_icsec , &
                            nelapse_in   = sim_duration,&
                            dtime_in     = step)
   end if

   call timemgr_init()

if (nscale == 0) then     ! nscale = 0 for global simulation, 1 for site simulation
   print *, 'BEPS run at global scale!'
    !! Reading boundary fields, yearly data,and soil Cpools for BEPS @J.Wang (note: if yearly and C pools data fields will change year by year, these datasets should be read in the time looping
   call read_boundary()
   call read_yrdata()
   call read_cpools()
else
    print *, 'BEPS run at site scale!'
    call read_boundary_site()   ! read site data, including yrdata, boundary data, and carbon pools
end if
   bfields  => bound
   climate  => clim
   psoil    => soilstat
   pp       => output
   ppar     => assim

call read_prior_para()       ! put the parameters to be optimized in a NETCDF
                                !file and read them as well their
                                !uncertainties,@MOUSONG.WU,2019-11
nd = 0

do     !! start time looping
   nd = nd + 1
   !-- iLab::inserted in order to pass 'calday' downstream
   caldy = get_curr_calday()
   call get_curr_date(yr,mn,dy,tod)
   doys = get_doys(yr)
   !-- iLab::inserted in order to pass ref_date downstream
   if( nd.eq.1 ) then
      yr_ref = yr
      mn_ref = mn
      dy_ref = dy
      tod_ref = tod
      write(ref_date(1:4),   '(i4.4)') yr_ref
      write(ref_date(5:5),   '(a)')    '-'
      write(ref_date(6:7),   '(i2.2)') mn_ref
      write(ref_date(8:8),   '(a)')    '-'
      write(ref_date(9:10),  '(i2.2)') dy_ref
      write(ref_date(11:19), '(a)')    'T00:00:00'
   endif
   !-- iLab::determine seconds elapsed since reference time
   call timemgr_diff_secs(yr_ref*10000+mn_ref*100+dy_ref, tod_ref, yr*10000+mn*100+dy, tod,&
            secs_elapsed)
    call get_CO2_concentration(yr,CO2_air)
    call get_COS_concentration(yr,COS_air)
!! change hourly input into daily input for further using model for long-term simulations, @MOUSONG.WU, 201905
if (meteo_input >= 0) then  ! call hourly meteo. input
    if (nscale == 0) then
        call read_meteo_hourly(yr, mn, dy, tod)
    else
        call timemgr_diff_secs(2010*10000+1*100+1,0,yr*10000+mn*100+dy,tod,&
            secs_meteo)
        n_meteo = int(secs_meteo/3600 + 1)
        call read_meteo_site(n_meteo)
    end if
else
    if(is_first_step() .or. is_end_curr_day()) then
      call read_meteo_daily(yr, mn, dy, tod)
    end if
end if

if (lai_input >=0) then
   if (is_first_step()) print *, 'lai is input!'
   if (is_first_step() .or. is_end_curr_day()) then
    if (nscale == 0) then
        call read_lai(yr, mn, dy, tod, caldy)
    else
        call read_lai_site(caldy)
    end if
   end if
else
  
 if (is_first_step()) then
    print *, 'lai is simulated with phenology scheme!'
 end if
end if

! Reading the Vcmax here for assimilation usage @J.Wang
!if(is_first_step() .or. is_end_curr_month()) then
!   if (nscale == 0) then
!    call read_Vcmax(yr, mn, dy, tod)
!   else
!    call read_Vcmax_site(yr, mn, dy, tod)
!   end if
!end if

!    call mpi_barrier(mpi_comm_world,ierr)
    do i = 1,npoints    !! spatial iteration

    !! calculate the solar zenith
    call s_coszs(yr, mn, dy, tod, caldy, doys, &
         bfields%latitude(i),bfields%longitude(i),CosZs,hr_loc,hr_arc)
    !!if(myid == 0) write(*,*) "hr_loc=",hr_loc

    !! retrieve meteo for this point
               meteo%Srad                   = climate%Srad(i)
               meteo%wind                   = climate%Wind(i)
               meteo%rainfall               = climate%Rain(i)
               meteo%snow                   = climate%Snow(i)
               meteo%rh                     = climate%Rh(i)
               if (meteo_input < 0) then
                meteo%tempmx                = climate%Tempmx(i)    ! read daily max and min temperatures instead
                meteo%tempmn                = climate%Tempmn(i)
               else
                meteo%temp                  = climate%Temp(i)
               end if

! .. compute daily course of temperature and daylength
               rdaymid = REAL (sim_duration+1) / 2.
               delta = -23.4*COS(2.*PI*(rdaymid+10.)/365.)
               spds = SIN(bfields%latitude(i)*pio180)*SIN(delta*pio180)
               cpds = COS(bfields%latitude(i)*pio180)*COS(delta*pio180)
               arg = -spds/cpds
               IF (arg>1.) THEN
               !        polar night:
                  daylen = 0.
               ELSE IF (arg<-1) THEN
               !        polar day:
                  daylen = 24.
               ELSE
               !        normal day / night:
                  daylen = ACOS(arg)/PI*24.
               END IF

!###########Compute subdaily temperature based on daily input,@MOUSONG.WU,201905#####################

if(meteo_input < 0) then
               ! .. compute average conditions
               atmean = (meteo%tempmx + meteo%tempmn)/ 2.
               atrange = meteo%tempmx - meteo%tempmn

               !  hour angle at sunset, added by MOUSONG.WU, 201905
               sunset_arc   = (daylen/2.)*2.0*PI/24.0

               IF (daylen>=4. .AND. daylen<=20.) THEN
               !        sunrise
                  h0 = 12. - daylen/2.
               !        sundown
                  h1 = 12. + daylen/2.
               !        at sundown:
                  sd1 = SIN(PI*(2.*h1+(daylen-52.)/2.)/(daylen+4.))
               !! unroll zum vektorisieren
                  IF (hr_loc>h0 .AND. hr_loc<h1) THEN
                     sd = SIN(PI*(2.*hr_loc+(daylen-52.)/2.)/(daylen+4.))
                     meteo%temp = atmean + atrange/2.*sd
                  ELSE
                  ! temperature at sundown
                    tmp1 = atmean + atrange/2.*sd1
                  ! hours since sundown
                    dhour = MOD(hr_loc-h1+24.,24.)
                    tmin = atmean - atrange/2.
                    meteo%temp = tmp1 - (tmp1-tmin)*(dhour/(24.-daylen))
                  END IF
               ELSEIF (daylen>20.) THEN
                    sd = COS(PI*(hr_loc-14.)/(daylen/2.+2.))
                    meteo%temp = atmean + atrange/2.*sd
               ELSE
                    meteo%temp = atmean
               END IF
               climate%Temp(i) = meteo%temp
! convert daily solar radiation into hourly using the method by M. Collares-Pereira and A. Rabl,
! �The average distribution of solar radiation-correlations between diffuse and hemispherical
!and between daily and hourly insolation values,� Solar Energy,vol. 22, no. 2, pp. 155�164, 1979.
               a = 0.409 + 0.5016 * SIN(sunset_arc - 60.)
               b = 0.6609 - 0.4767 * SIN(sunset_arc - 60.)
               meteo%Srad = meteo%Srad*(a+b*COS(hr_arc))*(PI/24.)*(COS(hr_arc)-COS(sunset_arc))/ &
                            &(SIN(sunset_arc)-(2*PI*sunset_arc/360.)*COS(sunset_arc))
end if

! Calculate cloud fraction, separate shortwave radiation
               if(CosZs < 0.001) then
                     ratio_cloud=0.
               else
                     ratio_cloud=meteo%Srad/(1367.*CosZs)
               end if
               if (ratio_cloud > 0.8) then
                     shortRad_df  = 0.13*meteo%Srad
               else
                     shortRad_df  = (0.943+0.734*ratio_cloud-4.9*ratio_cloud**2+1.796*ratio_cloud**3+2.058*ratio_cloud**4)*&
                                   & meteo%Srad
               end if
               shortRad_df   = min(shortRad_df,meteo%Srad)
               shortRad_df   = max(shortRad_df,0.)

               shortRad_dir  = meteo%Srad - shortRad_df

               meteo%S_dff   = shortRad_df
               meteo%S_dir   = shortRad_dir
               climate%Swdr(i)   = shortRad_dir
               climate%Swdf(i)   = shortRad_df

      do j = 1,PFT    !! PFT iteration
           if(bfields%lcno(i,j) > 0 .and. bfields%sw(i) >0. .and. bfields%stext(i) >0 .and. bfields%clumping(i) > 0.5) then
               call readparam(bfields%lcno(i,j),param)
               if(lai_input >= 0) then
                 lai = bfields%lai(i,j)
               else 
                 if (is_first_step()) then
                    lai = bfields%laiyr(i,j)
                    mid_res%lai_old = lai
                 else
                    lai = mid_res%lai_new
                 end if
               end if
               lai = lai*param(2)/bfields%clumping(i)
               call readcoef(bfields%lcno(i,j),bfields%stext(i),coef)

               if(nsrest == nsrStartup .and. is_first_step()) then
                  call Init_soil_parameters(bfields%lcno(i,j),bfields%stext(i), param(27),soilp)
                  soilp%r_drainage = param(26)
                  !soilp%r_drainage = ppar%p_drainage(j)            ! read this
                                                                    !para from NC file,@MOUSONG.WU,2019-11
                  ii = bfields%stext(i)
                 ! write(*,*) 'p_Ksat = ',ppar%p_Ksat(ii) 
                 ! write(*,*) 'Ksat_old = ',soilp%Ksat(0)
                  do kk = 0,4
                      soilp%Ksat(kk) = ppar%p_Ksat_scalar(ii)*soilp%Ksat(kk)
                      soilp%b(kk)    = ppar%p_b_scalar(ii)*soilp%b(kk)
                  end do
                  ! replace these three para. above with values from nc
                  ! file,@MOUSONG.WU,2019-11
                  call Init_soil_status(soilp,bfields%st(i),climate%Temp(i),bfields%sw(i),bfields%sdp(i))
                  do k = 0,40
                       var_o(k)   = 0.
                  end do
                  do k = 3,8
                       var_o(k)   = climate%Temp(i)
                  end do
                  do k = 9,14
                       var_o(k)   = soilp%temp_soil_p(k-9)
                  end do
                  do k = 21,26
                       var_o(k)   = soilp%thetam_prev(k-21)
                  end do
                  do k = 27,32
                       var_o(k)    = soilp%ice_ratio(k-27)
                  end do
               else if(nsrest == nsrContinue .and. is_first_step()) then
                  call Init_soil_parameters(bfields%lcno(i,j),bfields%stext(i),param(27),soilp)
                  soilp%r_drainage = param(26)
                  !soilp%r_drainage = ppar%p_drainage(j)  ! read from nc
                                                          ! file,MOUSONG.WU@2019-11
                  ii = bfields%stext(i)
                  do kk = 0,4
                      soilp%Ksat(kk) = ppar%p_Ksat_scalar(ii)*soilp%Ksat(kk)
                      soilp%b(kk)    = ppar%p_b_scalar(ii)*soilp%b(kk)
                  end do
                  !Replace with NC values,for optimization
                  !purpose,@MOUSONG.WU,2019-11
                  var_o(:)         = v2last(i,:,j)
                  do k= 9,14
                     soilp%temp_soil_c(k-9)  = var_o(k)
                  end do
                  do k= 21,26
                     soilp%thetam(k-21)      = var_o(k)
                  end do
!                 do k = 27,32
!                    soilp%ice_ratio(k-27)   = var_o(k)
!                 end do

               else
                  var_o(:)         = v2last(i,:,j)
                  call retrive_soilp(soilp,i,j,0)
               end if

               !!! simulating photosynthesis
 !              do llll = 0,40
 !                 write(*,*)  "DG004: Var_o = ",llll, var_o(llll)
 !              end do

               call inter_prg(yr, mn, dy, tod, &
                    lai,lai_input,bfields%lcno(i,j),bfields%clumping(i),ppar%p_Vcmax(j),ppar%p_VJ_slope(j),ppar%p_f_leaf,&
                    ppar%p_kc25,ppar%p_ko25,ppar%p_tau25,ppar%p_sif_alpha(j),ppar%p_sif_beta(j),param,meteo,CosZs,&
                    var_o,var_n,soilp,mid_res,daylen)
! CHANGE Vcmax read from Vcmax file with the Vcmax read from initial para. NC
! file, for optimization purpose,@MOUSONG.WU,2019-11
 !              do llll =0,40
  !                write(*,*)  "DG004: Var_n = ",llll,var_n(llll)
  !             end do
               v2last(i,:,j)  = var_n(:)
               call retrive_soilp(soilp,i,j,1)
               !!! simluating Ra
               call plant_resp(ppar%p_q10(j),bfields%lcno(i,j),mid_res,bfields%laiyr(i,j),lai,meteo%temp,&
                               soilp%temp_soil_c(1),CosZs)
!USE p_q10 here to adjust q10,p_q10 is read from initial para. NC file, for
!optimization purpose,@MOUSONG.WU,2019-11
               !!! simulating Rh
               Ccd(0)       = bfields%Ccd(i,j)
               Cssd(0)      = bfields%Cssd(i,j)
               Csmd(0)      = bfields%Csmd(i,j)
               Cfsd(0)      = bfields%Cfsd(i,j)
               Cfmd(0)      = bfields%Cfmd(i,j)
               Csm(0)       = bfields%Csm(i,j)
               Cm(0)        = bfields%Cm(i,j)
               Cs(0)        = bfields%Cm(i,j)
               Cp(0)        = bfields%Cp(i,j)
               ! to get soil texture for this point,@MOUSONG.WU,2019-11
               jj = bfields%stext(i)

               call soil_resp(Ccd,Cssd,Csmd,Cfsd,Cfmd,Csm,Cm,Cs,&
                             Cp,bfields%nppyr(i,j),coef,bfields%stext(i),soilp,mid_res)
               
               !! for output variables
               pp%GPPpft(i,j)   = mid_res%GPP*bfields%PCT_PFT(i,j)/100.
               pp%SIFpft(i,j)   = mid_res%SIF*bfields%PCT_PFT(i,j)/100.
               pp%NPPpft(i,j)   = mid_res%NPP*bfields%PCT_PFT(i,j)/100.
               pp%NEPpft(i,j)   = mid_res%NEP*bfields%PCT_PFT(i,j)/100.
               pp%SHpft(i,j)    = mid_res%SH*bfields%PCT_PFT(i,j)/100.
               pp%LHpft(i,j)    = mid_res%LH*bfields%PCT_PFT(i,j)/100.
               pp%Transpft(i,j) = mid_res%Trans*bfields%PCT_PFT(i,j)/100.
               pp%Evappft(i,j)  = mid_res%Evap*bfields%PCT_PFT(i,j)/100.
               pp%Net_Radpft(i,j)  = mid_res%Net_Rad*bfields%PCT_PFT(i,j)/100.
               pp%LAIpft(i,j)   = lai*bfields%PCT_PFT(i,j) /100.
               pp%Thetampft(i,j)   = mid_res%thetam_surf*bfields%PCT_PFT(i,j)/100.
               pp%fAPARpft(i,j) = mid_res%fAPAR*bfields%PCT_PFT(i,j)/100.
               pp%COS_fluxpft(i,j) = (mid_res%COS_plant+mid_res%COS_grnd)*bfields%PCT_PFT(i,j)/100.               
               
               pp%NPP_yr_acc(i,j) = pp%NPP_yr_acc(i,j) + mid_res%NPP*bfields%PCT_PFT(i,j)/100.*1.e-2*step       ! convert NPP to Mg/ha for calculation of VOD
               
               if (is_end_curr_year()) then
               ! calculate VOD (vegetation optical depth) with results derived from SMOS-IC product, @Mousong.Wu, 201905,taweff is a PFT specific parameter
                  pp%VODpft(i,j)   = ppar%p_agb2vod*atan(ppar%p_taweff(j)*pp%NPP_yr_acc(i,j)) + ppar%p_D0(j)*lai
                  pp%NPP_yr_acc(i,j) = 0.
!                  write(*,*) 'VOD= ',pp%VODpft(i,j)
               else
                  pp%VODpft(i,j) = 0.
               end if
              
               ! write(*,*) 'hr_loc = ', hr_loc
               !! calculate the OCO-2 SIF   across at 1:30pm
               !write(*,*) 'SIFpft= ',mid_res%SIF*bfields%PCT_PFT(i,j) 
               if(hr_loc >= 13. .and. hr_loc <14.) then
                   pp%SIFpft_sat(i,j)  = mid_res%SIF*bfields%PCT_PFT(i,j)/100.
               else
                   pp%SIFpft_sat(i,j)  = 0.
               end if
               !pp%SIFpft_sat(i,j) = max(pp%SIFpft_sat(i,j),0.)
               !write(*,*) 'SIFpft_sat = ', pp%SIFpft_sat(i,j)
           end if
      end do   !! end PFT loop

            pp%GPP(i)   = sum(pp%GPPpft(i,:))
            pp%SIF(i)   = sum(pp%SIFpft(i,:))
            pp%SIF_sat(i)  = sum(pp%SIFpft_sat(i,:))
            pp%NPP(i)   = sum(pp%NPPpft(i,:))
            pp%NEP(i)   = sum(pp%NEPpft(i,:))
            pp%SH(i)    = sum(pp%SHpft(i,:))
            pp%LH(i)    = sum(pp%LHpft(i,:))
            pp%Trans(i) = sum(pp%Transpft(i,:))
            pp%Evap(i)  = sum(pp%Evappft(i,:))
            pp%Net_Rad(i) = sum(pp%Net_Radpft(i,:))
            pp%LAI(i)     = sum(pp%LAIpft(i,:))
            pp%Thetam(i)  = sum(pp%Thetampft(i,:))
            pp%fAPAR(i)   = sum(pp%fAPARpft(i,:))
            pp%VOD(i)     = sum(pp%VODpft(i,:))
            pp%COS_flux(i)     = sum(pp%COS_fluxpft(i,:))

    end do      !! end spatial loop
!          call mpi_barrier(mpi_comm_world,ierr)
          !! advance time
          call advance_timestep()

          !!! write oout data and restart fields
          call av_output(yr, mn, dy, tod, get_nstep(), is_end_curr_month(), ref_date, secs_elapsed)

          !-- iLab::restart requires next time-step
          call get_curr_date(yr,mn,dy,tod)
          if(restart_frq <0) then   !! <0 ndays
             kount = get_nstep()
             rst_nstep  = -restart_frq*86400/step
             if(mod(kount,rst_nstep) ==0) call restart_io("write",yr,mn,dy,tod)
          else if(restart_frq ==0) then
             if(is_end_curr_month()) call restart_io("write",yr,mn,dy,tod)
          end if

          !!! determine whether it's the last step
           if(is_last_step()) exit

end do   !! end time loop

!! For clean up
call ESMF_Finalize()

!if(myid==0) 
write(6,*)   "BEPS run finished successfully!"
!call MPI_Finalize(ierr)
!stop
end program
