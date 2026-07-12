program tester
  use tblite_version, only : get_tblite_version
  implicit none
  character(len=:), allocatable :: version
  call get_tblite_version(string=version)
  print *, version
end program tester
