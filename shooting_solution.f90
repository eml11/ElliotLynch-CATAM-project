subroutine readinputfile (unit,filename)

  integer :: unit
  character(len=256) :: filename
  character (len=256) :: buf, label
  integer :: ios = 0
  integer :: pos = 0

  do while (ios.eq.0)

    read(unit, '(A)', iostat=ios) buf
    
    if (ios.eq.0) then
    
      pos = scan(buf, '    ')
      label = buf(1:pos)
      buf = buf(pos+1:)

      select case (label)
        case ()

        case ()

      end select

    end if
  end do
  
end subroutine

