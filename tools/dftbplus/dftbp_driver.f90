!> Milestone-1 driver: pass coordinates and settings to DFTB+ (library API)
!> and get the energy and forces back, the way qtes-dftb will use it.
!>
!> Reads, from the current directory:
!>   INxyz            qtes-dftb geometry file (6 parameter lines, then
!>                    "symbol x y z itype" per atom, in Angstrom)
!>   IN_dftbplus.hsd  DFTB+ settings in HSD format, without a Geometry block
!>
!> Usage:
!>   dftbp_driver            energy + gradient at the INxyz geometry, then
!>                           again with atom 1 moved by +0.1 bohr in z
!>   dftbp_driver fd [h]     also compare the gradient with central finite
!>                           differences of the energy (step h, default
!>                           1e-4 bohr)
!>
!> Writes:
!>   dftbp_input.hsd   the full input given to DFTB+ (incl. geometry, bohr);
!>                     usable as dftb_in.hsd for the dftb+ program
!>   dftbp_result.dat  energy and gradient at the INxyz geometry
!>   dftbp.log         DFTB+ screen output
program dftbp_driver
  use dftbplus, only : TDftbPlus, TDftbPlus_init, TDftbPlus_destruct, TDftbPlusInput, fnode,&
      & setChild, setChildValue, dumpHsd, getDftbPlusBuild, getDftbPlusApi
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  !> Angstrom per bohr, the value used by qtes-dftb (read_input_mpi.f)
  real(dp), parameter :: bohrAA = 0.52918_dp

  type(TDftbPlus) :: dftbp
  type(TDftbPlusInput) :: input
  type(fnode), pointer :: pRoot, pGeo

  character(len=2), allocatable :: symbols(:), typeNames(:)
  integer, allocatable :: species(:)
  real(dp), allocatable :: coords(:,:), coords0(:,:), grad(:,:), grad0(:,:), gradFd(:,:)
  real(dp) :: energy, energy0, ep, em, h, errMax
  integer :: nAtom, iAt, iCart, logUnit, outUnit, major, minor, patch
  character(len=32) :: arg
  character(:), allocatable :: version
  logical :: doFd

  ! ---- command line
  doFd = .false.
  h = 1.0e-4_dp
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    if (trim(arg) /= "fd") error stop "usage: dftbp_driver [fd [h]]"
    doFd = .true.
    if (command_argument_count() >= 2) then
      call get_command_argument(2, arg)
      read(arg, *) h
    end if
  end if

  ! ---- geometry from INxyz (Angstrom -> bohr)
  call readInxyz("INxyz", symbols, coords0)
  nAtom = size(symbols)
  call speciesFromSymbols(symbols, typeNames, species)
  allocate(coords(3, nAtom), grad(3, nAtom), grad0(3, nAtom), gradFd(3, nAtom))

  ! ---- DFTB+ instance; its screen output goes to dftbp.log
  open(newunit=logUnit, file="dftbp.log", action="write", status="replace")
  call getDftbPlusBuild(version)
  call getDftbPlusApi(major, minor, patch)
  write(*, "(A,A,A,I0,'.',I0,'.',I0)") "DFTB+ build: ", version, ", API ", major, minor, patch
  call TDftbPlus_init(dftbp, outputUnit=logUnit)

  ! ---- settings from IN_dftbplus.hsd, geometry added in code
  call dftbp%getInputFromFile("IN_dftbplus.hsd", input)
  call input%getRootNode(pRoot)
  call setChild(pRoot, "Geometry", pGeo)
  call setChildValue(pGeo, "Periodic", .false.)
  call setChildValue(pGeo, "TypeNames", typeNames)
  call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, nAtom]), coords0)

  open(newunit=outUnit, file="dftbp_input.hsd", action="write", status="replace")
  call dumpHsd(input%hsdTree, outUnit)
  close(outUnit)

  call dftbp%setupCalculator(input)
  if (dftbp%nrOfAtoms() /= nAtom) error stop "DFTB+ atom count differs from INxyz"

  ! ---- 1. energy and gradient at the INxyz geometry
  call calc(coords0, energy0, grad0)
  write(*, "(/,A,I0,A,*(1X,A))") "Atoms: ", nAtom, ", types:", (trim(typeNames(iAt)), iAt = 1, size(typeNames))
  call report("INxyz geometry", energy0, grad0)
  open(newunit=outUnit, file="dftbp_result.dat", action="write", status="replace")
  write(outUnit, "(A)") "# Mermin free energy (Ha), then gradient (Ha/bohr) per atom"
  write(outUnit, "(ES24.16)") energy0
  do iAt = 1, nAtom
    write(outUnit, "(A2,3ES24.16)") symbols(iAt), grad0(:, iAt)
  end do
  close(outUnit)

  ! ---- 2. same instance, new geometry: atom 1 moved by +0.1 bohr in z
  coords(:,:) = coords0
  coords(3, 1) = coords(3, 1) + 0.1_dp
  call calc(coords, energy, grad)
  call report("atom 1 moved by +0.1 bohr in z", energy, grad)
  write(*, "(A,ES14.6,A)") "Energy change: ", energy - energy0, " Ha"

  ! ---- 3. optional: gradient vs central finite differences
  if (doFd) then
    do iAt = 1, nAtom
      do iCart = 1, 3
        coords(:,:) = coords0
        coords(iCart, iAt) = coords0(iCart, iAt) + h
        call calc(coords, ep, grad)
        coords(iCart, iAt) = coords0(iCart, iAt) - h
        call calc(coords, em, grad)
        gradFd(iCart, iAt) = (ep - em) / (2.0_dp * h)
      end do
    end do
    errMax = maxval(abs(gradFd - grad0))
    write(*, "(/,A,ES10.2,A)") "Finite differences (h = ", h, " bohr):"
    write(*, "(A,ES12.4,A)") "  max |analytic - FD| = ", errMax, " Ha/bohr"
    write(*, "(A,ES12.4,A)") "  max |analytic|      = ", maxval(abs(grad0)), " Ha/bohr"
  end if

  call TDftbPlus_destruct(dftbp)
  close(logUnit)

contains

  !> One DFTB+ calculation: new coordinates in, energy and gradient out
  subroutine calc(xyz, e, g)
    real(dp), intent(in) :: xyz(:,:)
    real(dp), intent(out) :: e, g(:,:)
    call dftbp%setGeometry(xyz)
    call dftbp%getEnergy(e)
    call dftbp%getGradients(g)
  end subroutine calc


  subroutine report(label, e, g)
    character(*), intent(in) :: label
    real(dp), intent(in) :: e, g(:,:)
    integer :: i
    write(*, "(/,A,A)") "== ", label
    write(*, "(A,F20.12,A)") "Mermin free energy: ", e, " Ha"
    write(*, "(A)") "Gradient (Ha/bohr):"
    do i = 1, size(g, dim=2)
      write(*, "(I4,1X,A2,3F18.12)") i, symbols(i), g(:, i)
    end do
  end subroutine report


  !> Reads a qtes-dftb INxyz file: 6 parameter lines, then atoms until the
  !> end of the file (blank lines are skipped). Coordinates in Angstrom.
  subroutine readInxyz(fname, sym, xyz)
    character(*), intent(in) :: fname
    character(len=2), allocatable, intent(out) :: sym(:)
    real(dp), allocatable, intent(out) :: xyz(:,:)
    integer, parameter :: maxAtom = 10000
    character(len=2) :: s(maxAtom)
    real(dp) :: r(3, maxAtom)
    character(len=256) :: line
    integer :: u, ios, i, n
    open(newunit=u, file=fname, action="read", status="old")
    do i = 1, 6
      read(u, "(A)") line
    end do
    n = 0
    do
      read(u, "(A)", iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      if (n == maxAtom) error stop "too many atoms in INxyz"
      n = n + 1
      read(line, *) s(n), r(:, n)
    end do
    close(u)
    if (n == 0) error stop "no atoms found in INxyz"
    sym = s(1:n)
    xyz = r(:, 1:n) / bohrAA
  end subroutine readInxyz


  !> Type names in order of first appearance, and the species index per atom
  subroutine speciesFromSymbols(sym, names, spec)
    character(len=2), intent(in) :: sym(:)
    character(len=2), allocatable, intent(out) :: names(:)
    integer, allocatable, intent(out) :: spec(:)
    character(len=2) :: tmp(size(sym))
    integer :: i, j, n
    n = 0
    allocate(spec(size(sym)))
    do i = 1, size(sym)
      do j = 1, n
        if (tmp(j) == sym(i)) exit
      end do
      if (j > n) then
        n = n + 1
        tmp(n) = sym(i)
      end if
      spec(i) = j
    end do
    names = tmp(1:n)
  end subroutine speciesFromSymbols

end program dftbp_driver
