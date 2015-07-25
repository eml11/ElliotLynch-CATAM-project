subroutine readinputfile (unitinput,filename)

  use mod_shared

  integer :: unitinput
  character(len=256) :: filename
  character (len=256) :: buf, label
  integer :: ios = 0
  integer :: pos = 0

  open(unit,file=filename)

  do while (ios.eq.0)

    read(unitinput, '(A)', iostat=ios) buf
    
    if (ios.eq.0) then
    
      pos = scan(buf, '    ')
      label = buf(1:pos)
      buf = buf(pos+1:)

      select case (label)
        case ('RADIUSS')

        case ('GAMMA')

        case ('MASSM')

        case('MASST')

        case('PRESSUREC')

        

      end select

    end if
  end do
  
  close(unitinput)

end subroutine
