program brute_force_integration
implicit none

integer, parameter :: ik = 4 
integer, parameter :: rk = 8
integer(ik), parameter :: NVMAX = 10001
integer(ik), parameter :: NMATERIALMAX = 10
integer(ik), parameter :: NTEMPMAX = 10
real(rk), parameter :: pi = 3.1415926535897932384626433832795028841972_rk
real(rk), parameter :: c = 1.602176565E-19_rk ! [J/eV]
real(rk), parameter :: kb = 1.3806503E-23_rk ! Boltzmann Constant [J/K]
real(rk), parameter :: kc = (1.0E+03_rk)*(kb/c) ! [meV/K]
!-------------------------------------------------------------------------------------
integer(ik) :: i, j, k, l, ne, nv, nmaterials, ntemp
real(rk) :: vmin, vmax, e, de, z, v(NVMAX), dv, emin, emax
real(rk) :: temp(NMATERIALMAX,NTEMPMAX), delta_real(NMATERIALMAX,NTEMPMAX), delta_imag(NMATERIALMAX,NTEMPMAX)
real(rk) :: integral(NMATERIALMAX,NTEMPMAX,NVMAX)
character(LEN = 32) :: material_name(NMATERIALMAX), filename

! Get simulation parameters
call param_input(ne, emin, emax, z, nv, vmin, vmax)

! Get material parameters
call material_param_input(nmaterials, ntemp, material_name, temp)

dv = (vmax - vmin)/real(nv - 1, rk)
de = (emax - emin)/real(ne - 1, rk)

! Loop over the nmaterials values of the materials
do i = 1, nmaterials

  call get_material_data(ntemp, material_name(i), temp(i,:), delta_real(i,:), delta_imag(i,:))

  write(filename, '(3a)') "graphing_", trim(material_name(i)), ".dat"
  open(unit = 10, file = trim(filename), status = "unknown")

  ! Loop over the ntemp values of the temperature
  do j = 1, ntemp

    ! Loop over the nv values of the voltage
    do k = 1, nv

      ! Compute the value of the voltage
      v(k) = vmin + real(k - 1, rk)*dv

      ! Initialize the value of the integral and weightArea, and the starting energy
      integral(i, j, k) = 0.0_rk

      e = emin
       
        ! Numerically integrate the integrand function between emin and emax
        do l = 0, ne - 1

          ! Update the value of the integral
          integral(i, j, k) = integral(i, j, k) + f(e, v(k), temp(i, j), delta_real(i, j), delta_imag(i, j), z)*de
      
          ! Update the value of the energy
          e = e + de

        end do

    end do

  end do

  do k = 1, nv

    write(10,'(11es23.15)') v(k), ((integral(i, j, k)*(1.0_rk + (z)**2))/delta_real(i, j), j = 1, ntemp)

  end do

  close(unit = 10)
  
end do

contains

  !---------------------------------------------------------------------------
  function f(e, v, temp, delta_real, delta_imag, z) result(fresult)
  !----------------------------------------
  ! The integrand function
  !----------------------------------------
  implicit none

  ! Variable declarations
  real(rk) :: e, v, temp, delta_real, delta_imag, z
  real(rk) :: fresult

    fresult = (fdirac(e - v, temp) - fdirac(e, temp))*fprob(e, delta_real, delta_imag, z)

  end function f

  !----------------------------------------------------------
  function fprob(e, delta_real, delta_imag, z) result(fresult)
  implicit none
  
  real(rk) :: e, delta_real, delta_imag, z, fresult
  real(rk) :: pa1, pa2, alfa, beta, theta, gammasq, An, Br
  complex :: usq, vsq, usqc, vsqc

! A set of dummy variables containing the real and imaginary parts of the energy gap. These are incorporated in probability currents. 

   pa1 = e**2 - delta_real**2 + delta_imag**2
   pa2 = 2_rk*delta_real*delta_imag
   theta = 0.5_rk*atan2(pa2,pa1)
   alfa = (((pa1**2 + pa2**2)**0.25)*(cos(theta)))/abs(e) !Real Part of u^2 and v^2
   beta = (((pa1**2 + pa2**2)**0.25)*(sin(theta)))/abs(e) !Imaginary Part of u^2 and v^2
   usq = (((2_rk)**0.5/2_rk)*(((1_rk + alfa)**2 + beta**2)**0.25))*cmplx(cos(0.5_rk*atan2(beta,(1_rk + alfa))),sin(0.5_rk*atan2(beta,(1_rk + alfa))))
   vsq = (((2_rk)**0.5/2_rk)*(((1_rk - alfa)**2 + beta**2)**0.25))*cmplx(cos(0.5_rk*atan2(-beta,(1_rk - alfa))),sin(0.5_rk*atan2(-beta,(1_rk - alfa))))
   usqc = (((2_rk)**0.5/2_rk)*(((1_rk + alfa)**2 + beta**2)**0.25))*cmplx(cos(0.5_rk*atan2(beta,(1_rk + alfa))),-sin(0.5_rk*atan2(beta,(1_rk + alfa))))
   vsqc = (((2_rk)**0.5/2_rk)*(((1_rk - alfa)**2 + beta**2)**0.25))*cmplx(cos(0.5_rk*atan2(-beta,(1_rk - alfa))),-sin(0.5_rk*atan2(-beta,(1_rk - alfa))))
  
    gammasq = (usq**2 + (usq**2 - vsq**2)*(z**2))*(usqc**2 + (usqc**2 - vsqc**2)*(z**2))

    An = (usq*vsq*usqc*vsqc)/gammasq !Andreev Reflection----------------

    Br = ((usqc**2 - vsqc**2)*(usq**2 - vsq**2)*(z**2 + z**4))/gammasq !Normal Reflection----------------

    fresult = 1.0_rk + An - Br

  end function fprob

  !----------------------------------------------------------------------------
  function fdirac(e, temp) result(fresult)
  !------------------------------------
  ! Fermi-Dirac distribution 
  !------------------------------------
  implicit none

  ! Variable declarations
  real(rk) :: e, temp
  real(rk) :: fresult
  
    if(abs(e) > 1.0E-03 .or. abs(temp) >= 1.0E-03) then

    fresult = 1.0_rk/(exp(e/(kc*temp)) + 1.0_rk)

  else

    fresult = 1.0_rk

  end if
  
  
  end function fdirac

  !----------------------------------------------------------
  subroutine param_input(ne, emin, emax, z, nv, vmin, vmax)
  implicit none

  integer(ik) :: ne, nv
  real(rk) :: vmin, vmax
  real(rk) :: z, emin, emax

   namelist / params / ne, emin, emax, z, nv, vmin, vmax
   open(unit = 10, file = 'param.in', status = 'old')
   read(10, NML = params) 
   close(unit = 10)

  end subroutine param_input

  !----------------------------------------------------------
  subroutine material_param_input(nmaterials, ntemp, material_name, temp)
  implicit none

  integer(ik) :: nmaterials, ntemp
  real(rk) :: temp(:,:)
  character(LEN = *) :: material_name(:)

  integer(ik) :: i, j

   open(unit = 10, file = 'materials_param.in', status = 'old')
   read(10,*) nmaterials, ntemp
   do i = 1, nmaterials
     read(10,*) material_name(i), (temp(i, j), j = 1, ntemp)
   end do

   close(unit = 10)

  end subroutine material_param_input

  !----------------------------------------------------------
  subroutine get_material_data(ntemp, material_name, temp, delta_real, delta_imag)
  implicit none

  real(rk), parameter :: epsilon = 1.0e-6_rk

  integer(ik) :: ntemp
  real(rk) :: temp(:), delta_real(:), delta_imag(:)
  character(LEN = *) :: material_name

  integer(ik) :: j, ioerr
  real(rk) :: temp_tmp, temp_min = 1.0e30_rk, deltar_tmp, deltai_tmp
  character(LEN = 32) :: filename

   write(filename, '(2a)') trim(material_name), ".dat"
   open(unit = 20, file = trim(filename), status = "unknown")

   j = 1

   do
     read(20,*,iostat = ioerr) temp_tmp, deltar_tmp, deltai_tmp
     if(ioerr > 0) then
       write(*,*) 'Check input, something is wrong.'
       exit
     else if(ioerr < 0) then
       write(*,*) 'End of material data file reached.'
       exit
     else
       if(abs(temp(j) - temp_tmp) < epsilon) then
         delta_real(j) = deltar_tmp
         delta_imag(j) = deltai_tmp        
         j = j + 1
       end if
       if(j > ntemp) then
         !write(*,*) temp(:), material_name, delta_real(:), delta_imag(:)
         write(*,'(2a)') 'All temperature values located in file: ', trim(filename)
         exit
       end if
     end if
   end do

   close(unit = 20)

   end subroutine get_material_data
!!
end program brute_force_integration
