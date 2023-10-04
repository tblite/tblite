module tblite_post_processing
   use mctc_env, only : wp, error_type, fatal_error
    use tblite_xtbml_features, only : xtbml_type, new_xtbml_features
    use tblite_param_xtbml_features, only : xtbml_features_record
    use tblite_post_processing_type, only : post_processing_type
    use tblite_param_post_processing, only : post_processing_param_list
    use tblite_toml, only : toml_error, toml_parse, toml_table, get_value
    use tblite_param, only : param_record
    use tblite_param_serde, only : serde_record
interface new_post_processing
    procedure :: new_post_processing_param
    procedure :: new_post_processing_cli
end interface
contains
subroutine new_post_processing_param(self, param)
   !> Instance of the xTB evaluator
   class(post_processing_type), allocatable, intent(inout) :: self
   class(post_processing_param_list) :: param
  do i = 1, size(param%list)
   select type(par => param%list(i)%record)
     type is (xtbml_features_record)
           block 
              type(xtbml_type), allocatable :: tmp_ml
              allocate(tmp_ml)
              call new_xtbml_features(par, tmp_ml)
              call move_alloc(tmp_ml, self)
          end block
    end select
  end do
end subroutine

subroutine new_post_processing_cli(self, config, error)
    class(post_processing_type), allocatable, intent(inout) :: self
    type(error_type), intent(inout) , allocatable:: error
    character(len=:), allocatable :: config
    class(post_processing_param_list), allocatable :: param
    allocate(param)
    select case(config)
    case("xtbml")
        block 
            type(xtbml_features_record), allocatable :: ml_param
            class(serde_record), allocatable :: cont
            allocate(ml_param)
            call populate_default_param(ml_param, .false.)
            call move_alloc(ml_param, cont)
            call param%push(cont)
        end block
    case("xtbml_xyz")
      block 
        type(xtbml_features_record), allocatable :: ml_param
        class(serde_record), allocatable :: cont 
        allocate(ml_param)
        call populate_default_param(ml_param, .true.)
        call move_alloc(ml_param, cont)
        call param%push(cont) 
      end block
    case default
        block
            type(toml_table), allocatable :: table
            integer :: io, stat
            type(xtbml_features_record), allocatable :: ml_param 
            type(toml_error), allocatable :: t_error
            type(param_record) :: record
            type(toml_table), pointer :: child

            open(file=config, newunit=io, status="old")
            
            call toml_parse(table, io, t_error)
            close(io)
            if (allocated(t_error)) then
                allocate(error)
                call move_alloc(t_error%message, error%message)
                return
             end if
            if (allocated(error)) then
                call fatal_error(error, "File name provided could not be parsed as a toml table")
            end if
            call get_value(table, "post-processing", child, requested=.false.)
            if (associated(child)) then
              call param%load(child, error)
            end if
        end block
    end select
    
    call new_post_processing(self, param) 
end subroutine

subroutine populate_default_param(param, tensor)
    type(xtbml_features_record), intent(inout) :: param    
    logical, optional :: tensor

      !> Compute geometry-based xtbml features
    param%xtbml_geometry = .true.
      !> Compute density-based xtbml features
    param%xtbml_density = .true.
      !> Return vectorial information additional to norm of the corresponding multipole moments
    if (present(tensor)) then 
      param%xtbml_tensor = tensor
    else
      param%xtbml_tensor = .false.
    end if
      !> Compute orbital energy based xtbml features
    param%xtbml_orbital_energy = .true.
      !> Compute energy based features, necessary for partitioning weights
    param%xtbml_energy = .true.
      !> Compute extended feature i.e. do CN weigthed real space convolution
    param%xtbml_convolution = .true.
      !> Scaling for logistic function, convolution over an array of values is supported
    param%xtbml_a = [1.0_wp]

end subroutine

end module