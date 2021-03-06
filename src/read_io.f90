!> Subroutine to read and process 
!! input file
!! @param unitinput file unit for reading
!! @param fname filename to be read
subroutine readinputfile (unitinput,fname)

  use mod_shared

  implicit none

  integer :: unitinput
  character(len=*) :: fname
  character (len=256) :: buf, label
  integer :: ios = 0
  integer :: pos = 0

  open(unit=unitinput,file=fname)

  do while (ios.eq.0)

    read(unitinput, '(A)', iostat=ios) buf

    if (ios.eq.0) then
    
      pos = scan(buf, '    ')
      label = buf(1:pos)
      buf = buf(pos+1:)

      select case (label)
        case('SHOOTTYPE')
          select case(ADJUSTL(buf))
            case ('IN')
              stype = 2
            case ('OUT')
              stype = 1
            case ('INOUT')
              stype = 0
          end select 
        case ('RADIUSS')
          read(buf, *, iostat=ios) radiusstar
        case ('GAMMA')
          read(buf, *, iostat=ios) ratiospecificheat
        case ('MASSM')
          read(buf, *, iostat=ios) massm
        case('MASST')
          read(buf, *, iostat=ios) masst
        case('PRESSUREC')
          read(buf, *, iostat=ios) pressurec
        case('DENSITYC')
          read(buf, *, iostat=ios) densityc
        case('TEMPC')
          read(buf, *, iostat=ios) tempc
        case('HMASSF')
          read(buf, *, iostat=ios) hmassfrac
        case('LUM0')
          read(buf, *, iostat=ios) luminosity0
      end select

    end if
  end do
  
  close(unitinput)

end subroutine
