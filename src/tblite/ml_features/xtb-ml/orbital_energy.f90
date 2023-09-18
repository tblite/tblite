module tblite_xtbml_orbital_energy
   use mctc_env, only : wp
   use tblite_xtbml_feature_type, only : xtbml_feature_type
   use tblite_wavefunction_type, only : wavefunction_type
   use mctc_io, only : structure_type
   use tblite_integral_type, only : integral_type
   use tblite_basis_type, only : basis_type 
   use tblite_container, only : container_cache
   use tblite_context , only : context_type
   use tblite_double_dictionary, only : double_dictionary_type
   use tblite_xtbml_convolution, only : xtbml_convolution_type
   implicit none
   private

   type, public, extends(xtbml_feature_type) :: xtbml_density_features_type
      character(len=29) :: label = "orbital energy-based features"

    contains
      procedure :: compute_features
      procedure :: compute_extended
      procedure, private :: allocate
      procedure, private :: allocate_extended
    
    end type

end module 