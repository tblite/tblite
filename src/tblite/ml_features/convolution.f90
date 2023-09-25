module tblite_ml_feature_convolution
    use mctc_env, only : wp
    implicit none
    private

    type, abstract, public :: convolution_type
        !> Kernel values to give scaling of property of neighbors based on atoom in question
        real(wp), allocatable :: kernel(:, :, :)
        character(len=:), allocatable :: label
    contains
        procedure :: info
    end type
contains

pure function info(self, verbosity, indent) result(str)
    !> Instance of the interaction container
    class(convolution_type), intent(in) :: self
    !> Verbosity level
    integer, intent(in) :: verbosity
    !> Indentation level
    character(len=*), intent(in) :: indent
    !> Information on the container
    character(len=:), allocatable :: str

    if (allocated(self%label)) then
        str = indent // self%label
    else
        str = "Unknown"
    end if
end function info

end module