module tblite_api_post_processing
    use mctc_env, only : error_type
    use tblite_post_processing_list, only : post_processing_list, post_processing_type
    use tblite_api_version, only : namespace
    use tblite_post_processing_list, only : new_post_processing
    use tblite_api_utils, only : c_f_character
    use iso_c_binding
    implicit none
    private

    public :: vp_post_processing, new_post_processing_api, delete_post_processing_api

!> Void pointer to a container instance
    type :: vp_post_processing
        !> Actual container
        type(post_processing_list) :: ptr
    end type vp_post_processing
    
    logical, parameter :: debug = .true.

contains

    function new_post_processing_api(charptr) result(vpost_proc)&
        & bind(C, name=namespace//"new_post_processing")
        type(vp_post_processing), pointer :: post_proc
        character(kind=c_char), intent(in) :: charptr(*)
        type(c_ptr) :: vpost_proc
        character(len=:), allocatable :: config_str
        type(error_type), allocatable :: error
        class(post_processing_type), allocatable :: pproc 
        if (debug) print '("[Info]", 1x, a)', "new_post_processing"
        call c_f_character(charptr, config_str)
        write(*,*) config_str
        allocate(post_proc) 
        call new_post_processing(pproc, config_str, error)
        call post_proc%ptr%push(pproc)
        write(*,*) size(post_proc%ptr%list)
        vpost_proc = c_loc(post_proc)

    end function

    subroutine delete_post_processing_api(vpost_proc) &
        & bind(C, name=namespace//"delete_post_processing")
     type(c_ptr), intent(inout) :: vpost_proc
     type(vp_post_processing), pointer :: post_proc
  
     if (debug) print '("[Info]", 1x, a)', "delete_post_processing"
  
     if (c_associated(vpost_proc)) then
        call c_f_pointer(vpost_proc, post_proc)
  
        deallocate(post_proc)
        vpost_proc = c_null_ptr
     end if
  end subroutine delete_post_processing_api

end module 
