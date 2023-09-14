module tblite_ml_features_type
    use mctc_env, only : wp
    use tblite_double_dictionary, only : double_dictionary_type
    use mctc_io, only : structure_type
    use tblite_basis_type, only : basis_type
    use tblite_results, only : results_type
    use tblite_integral_type, only : integral_type
    use tblite_wavefunction_type, only : wavefunction_type
    use tblite_basis_type, only : basis_type
    use tblite_container, only : container_cache
    use tblite_results, only : results_type
    use tblite_context, only : context_type
    implicit none
    private
    public :: ml_features_type

type, abstract :: ml_features_type
    
contains
    !setup container
    procedure :: compute
    procedure :: pack_res
    procedure :: print_csv
    !print toml could be possible
end type

contains

    subroutine compute(self, mol, wfn, integrals, bas, ccache, dcache, rcache, ctx, prlevel, dict)
            class(ml_features_type),intent(in) :: self
            !> Molecular structure data
            type(structure_type), intent(in) :: mol
            !> Wavefunction strcuture data
            type(wavefunction_type), intent(in) :: wfn
            !> integral container for dipole and quadrupole integrals for CAMMs
            type(integral_type) :: integrals
            type(basis_type), intent(in) :: bas
            !> Context container for writing to stdout
            type(context_type), intent(inout) :: ctx
            !> Compute cache containers
            type(container_cache),intent(inout) :: ccache, dcache, rcache
            integer :: prlevel
            type(double_dictionary_type), intent(inout) :: dict
    end subroutine

    
subroutine pack_res(self, res)
    class(ml_features_type),intent(inout) :: self
    type(results_type) :: res
end subroutine

subroutine print_csv(self, mol, dict)
    class(ml_features_type),intent(in) :: self
    type(structure_type) :: mol
    type(double_dictionary_type) :: dict
    integer :: n, i, out, j, nat
    character(len=:), allocatable :: tmp_label
    real(wp), allocatable :: tmp_array(:, :), array(:)
    integer, allocatable :: z_array(:)
    n = dict%get_n_entries()
    out = 42
    
    open (file='ml_feature_tblite.csv', newunit=out)
    
    allocate(z_array(mol%nat))
    do i=1, mol%nat
        z_array = mol%num(mol%id(i))
    end do
    write(out, '(a)', advance="no") trim("Atom")//','
    allocate(tmp_array(n, mol%nat))
    do i=1, n
      call dict%get_label(i, tmp_label)
      write(out, '(a)', advance="no") trim(tmp_label)//','
      if (i == n) write(out, '(a)', advance="no") trim(tmp_label) 
      call dict%get_entry(i, array)
      tmp_array(i, :) = array
    end do

    
    
    do j = 1, mol%nat
        write(out, '(i2,a)', advance="no") z_array(j), ','
        do i = 1, n
            write(out, '(f14.8,a)', advance="no") tmp_array(i, j), ','
            if (i == n) write(out, '(f14.8)') tmp_array(i, j) 
      end do
    end do

end subroutine
end module tblite_ml_features_type