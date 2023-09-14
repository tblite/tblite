module tblite_ml_features
    use mctc_env, only : wp, error_type, fatal_error
    use tblite_ml_features_type, only : ml_features_type
    use tblite_xtbml_features, only : xtbml_type, new_xtbml_features
    use tblite_param_ml_features, only : ml_features_record
interface new_ml_features
    procedure :: new_ml_features_param
    procedure :: new_ml_features_cli
end interface
contains
subroutine new_ml_features_param(self, ml_param)
    use tblite_ml_features_methods, only : ml_features_method
    use tblite_param_ml_features, only : ml_features_record
    !> Instance of the xTB evaluator
    class(ml_features_type), allocatable, intent(inout) :: self
    type(ml_features_record) :: ml_param
    
        select case(ml_param%ml_features)
        case(ml_features_method%xtbml)
           block 
              type(xtbml_type), allocatable :: tmp_ml
              call new_xtbml_features(ml_param, tmp_ml)
              call move_alloc(tmp_ml, self)
           end block
        end select
end subroutine

subroutine new_ml_features_cli(self, config, error)
    class(ml_features_type), allocatable, intent(inout) :: self
    type(error_type), intent(inout) , allocatable:: error
    character(len=:), allocatable :: config
    
    select case(config)
    case("xtbml")
        block 
            type(xtbml_type), allocatable :: tmp_ml
            type(ml_features_record) :: ml_param
            call populate_default_param(ml_param)
            call new_xtbml_features(ml_param, tmp_ml)
            call move_alloc(tmp_ml, self)
        end block
    case default
        block
            use tblite_toml, only : toml_error, toml_parse, toml_table
            use tblite_param_ml_features, only : ml_features_record
            type(toml_table), allocatable :: table
            integer :: io
            type(ml_features_record) :: param
            type(toml_error), allocatable :: t_error
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
            call param%load_from_toml(table, error)
            call new_ml_features(self, param)
            
        end block
    end select
end subroutine

subroutine populate_default_param(param)
    type(ml_features_record) :: param    
    param%ml_features = 1

      !> Compute geometry-based xtbml features
    param%xtbml_geometry = .true.
      !> Compute density-based xtbml features
    param%xtbml_density = .true.
      !> Return vectorial information additional to norm of the corresponding multipole moments
    param%xtbml_tensor = .true.
      !> Compute orbital energy based xtbml features
    param%xtbml_orbital_energy = .false.
      !> Compute energy based features, necessary for partitioning weights
    param%xtbml_energy = .false.
      !> Compute extended feature i.e. do CN weigthed real space convolution
    param%xtbml_convolution = .true.
      !> Scaling for logistic function, convolution over an array of values is supported
    param%xtbml_a = (1.0_wp)

end subroutine

end module