!       driver calling nr optimisation routine
!       with reverse mode derivative
!       ilab 2017
module sizes
  real(kind=8), parameter :: missing_value = -99999._8
  integer nn, mm
end module sizes
program runopti
  use sizes
  implicit none
  character(len=*), parameter :: prog = 'runopti'
  character(len=*), parameter :: trace_fname = 'rundfprv-trace.asc'
  integer(kind=4)           :: iter, n, m, itmax, iostat
  logical :: ok
  real (kind=8)              :: f, gtol, pert
  real (kind=8), allocatable :: x0(:), x(:), sx(:)
  logical :: prior_enabled = .true. !-- whether prior term is active in cost function
  logical :: exist
  external func, dfunc

  ! handling cmdline
  integer::narg,cptArg !#of arg & counter of arg
  integer :: nconsumed
  character(len=32) :: argname !Arg name
  character(len=32) :: argval

  !-- default settings
  !-- perturbation of prior control vector for initial guess
  pert = 0._8

  !Check if any arguments are found
  narg=command_argument_count()
  !Loop over the arguments
  if(narg>0)then
     !loop across options
     nconsumed=0
     arg: do cptArg=1,narg
        if( nconsumed.gt.0 ) then
           nconsumed = nconsumed - 1
           cycle arg
        endif
        call get_command_argument(cptArg,argname)
        select case(adjustl(argname))
        case("--help","-h")
           write(*,'(/)')
           write(*,*)" This is program '"//prog//"'"
           call dump_options()
           stop 0
           exit arg
        case('--pert')
           call get_command_argument(cptArg+1, argval)
           nconsumed = 1
           read(argval, *) pert
        case('--disable_prior')
           prior_enabled = .false.
        case default
           write(*,*)"Option ",adjustl(argname),"unknown"
           write(*,'(a)') "Available options:"
           call dump_options()
           stop 0
        end select
     end do arg
  end if

  ! init configuration
  call initf_bw(n,m)
  nn = n
  mm = m
  ! init unknowns
  allocate (x0(n), x(n), sx(n))

  !-- read prior
  call initx_cd(n,x0,sx)
  !-- initial control vector
  inquire(file='x.b', exist=exist)
  if( exist ) then
     write(*, '(a)') ' INFO::'//prog//': initial control vector is read from file ***'//&
          'x.b'//'***'
     open(unit=1, file='x.b', form='unformatted')
     read(1) x
     close(1)
  else
     write(*, '(a,e10.4,a)') ' INFO::'//prog//': apply perturbation of ', pert, &
          ' to prior control vector'
     x = x0 * (1+pert)
  endif

  if( prior_enabled ) then
     write(*, '(a)') ' INFO::'//prog//': prior-term active in cost function.'
     call enable_prior()
  else
     write(*, '(a)') ' INFO::'//prog//': prior-term de-activated in cost function.'
     call disable_prior()
  endif
  
  ! intialise optimisation parameters
  iter = 0
  gtol = 1.e-8
  f = 0.
  itmax=200
  open(unit=3, file=trace_fname, form='formatted', action='write', iostat=iostat)
  if( iostat.ne.0 ) then
     write(*, '(a)') ' FATAL::'//prog//': trace file ***'//trace_fname//'*** could not be '//&
          'opened for writing!'
     stop
  endif
  ! output table
  print '(3a20)', 'f', 'norm g', 'x'
  call dfpmin(x,n,gtol,iter,f,func,dfunc)
  print*, 'dfpmin, iterations: ', iter
  call woptimum(n,x0,sx,x)
  !-- close trace file
  close(3)
  write(*, '(a)') ' INFO::'//prog//': generated trace file ***'//trace_fname//'***'
  ok = (iter.lt.itmax)
  call finishc(n,x,f,ok)

  !-- dispose resources
  deallocate(x0, x, sx)
contains

  subroutine dump_options()
    implicit none
    character(len=15) :: option
    write(*,'(a)') '=============================='
    write(*,*) " Options:"
    write(*,'(/)')
    option = "--pert"
    write(*,'(2x, a15,2x,a,e10.4,a)') option,&
         &"perturbation of prior to set initial control vector "&
         &//"(default:",pert,")"
    option = "--disable_prior"
    write(*,'(2x, a15, 2x, a)') option,&
         &"whether to deactivate prior-term in cost function."
    write(*,'(2/)')
  end subroutine dump_options
end program runopti

function func(x)
  use sizes
  implicit none
  real ( kind = 8 ) :: x(nn), func, f
  call cost_cd(nn,x,mm, f)
  func=f
  print '(e20.6,a20,6e20.6)', f, '-', x(1:min(3,nn))
  write(3,*)  f, missing_value, x(1:min(100,nn))
end function func

subroutine dfunc(x,fd)
  use sizes
  implicit none
  real ( kind = 8 ) :: x(nn), fd(nn), f, fb 
  fd = 0.
  fb = 1.
  call COST_bw(nn, x, fd, mm, f, fb)
  print '(a20,6e20.6)', '-', sqrt(sum(fd**2)), x(1:min(3,nn))  ! reverse model provides no value of f
  write(3,*)  missing_value, sqrt(sum(fd**2)), x(1:min(100,nn))
end subroutine dfunc
