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
          if (buf.eq.'IN') then
            stype = 2
          else if (buf.eq.'OUT') then
            stype = 1
          else if (buf.eq.'INOUT') then
            stype = 0
          end if 
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
