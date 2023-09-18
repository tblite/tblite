module tblite_ml_features_methods
    implicit none
    private
    public :: ml_features_method
    type, public :: ml_features_enum
        !> XTBML features
        integer :: xtbml = 1
    end type

    type(ml_features_enum), parameter :: ml_features_method = ml_features_enum()
end module 