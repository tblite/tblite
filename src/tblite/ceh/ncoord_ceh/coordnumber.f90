module tblite_ceh_coordination
    use mctc_env, only: wp
    use mctc_io, only: structure_type

    implicit none

    private

    public :: ceh_ncoord

    contains

    subroutine ceh_ncoord(cn)
        implicit none
        integer, intent(out) :: cn
        integer :: ierr

        cn = 0.0_wp
        ! call mctc_get_structure_type(ncoord, ierr)
        if (ierr /= 0) then
            write(*,*) 'ceh_ncoord: error getting structure type'
            stop 1
        end if
    end subroutine ceh_ncoord

end module tblite_ceh_coordination