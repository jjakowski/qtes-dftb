!> Milestone-1 driver: pass coordinates and settings to DFTB+ (library API)
!> and get the energy and forces back, the way qtes-dftb will use it.
!>
!> Reads, from the current directory:
!>   INxyz            qtes-dftb geometry file (6 parameter lines, then
!>                    "symbol x y z itype" per atom, in Angstrom, up to
!>                    the first blank line; INdyn, if present, checks nat)
!>   IN_dftbplus.hsd  DFTB+ settings in HSD format, without a Geometry block
!>
!> Usage:
!>   dftbp_driver            energy + gradient at the INxyz geometry, then
!>                           again with atom 1 moved by +0.1 bohr in z
!>   dftbp_driver fd [h]     also compare the gradient with central finite
!>                           differences of the energy (step h, default
!>                           1e-4 bohr)
!>   dftbp_driver series FILE [file|byqt|shuffle] [fresh]
!>                           energy + gradient for every geometry in FILE, a
!>                           qtes-dftb qt.dat (2 lines per QT per large step:
!>                           positions in bohr + w, then momenta + s). Order:
!>                           file = as in FILE (one rank holding all QTs),
!>                           byqt = each QT through all steps, shuffle = fixed
!>                           pseudo-random order. fresh = new DFTB+ instance
!>                           for every geometry (SCC starts from neutral atoms)
!>
!> Writes:
!>   dftbp_input.hsd   the full input given to DFTB+ (incl. geometry, bohr);
!>                     usable as dftb_in.hsd for the dftb+ program
!>   dftbp_result.dat  energy and gradient at the INxyz geometry
!>   dftbp.log         DFTB+ screen output
!>   series mode: series_result.dat (per geometry: index, step, QT, call
!>                number, energy, time) and series_grad.dat (gradient per
!>                geometry, in file order)
program dftbp_driver
  use dftbplus, only : TDftbPlus, TDftbPlus_init, TDftbPlus_destruct, TDftbPlusInput, fnode,&
      & setChild, setChildValue, dumpHsd, getDftbPlusBuild, getDftbPlusApi
  implicit none

  integer, parameter :: dp = kind(1.0d0)

  !> Angstrom per bohr, the value used by qtes-dftb (read_input_mpi.f)
  real(dp), parameter :: bohrAA = 0.52918_dp

  type(TDftbPlus) :: dftbp

  character(len=2), allocatable :: symbols(:), typeNames(:)
  integer, allocatable :: species(:)
  real(dp), allocatable :: coords(:,:), coords0(:,:), grad(:,:), grad0(:,:), gradFd(:,:)
  real(dp) :: energy, energy0, ep, em, h, errMax
  integer :: nAtom, iAt, iCart, logUnit, outUnit, major, minor, patch
  character(len=256) :: arg, seriesFile
  character(len=16) :: order
  character(:), allocatable :: version
  logical :: doFd, doSeries, fresh

  ! ---- command line
  doFd = .false.
  doSeries = .false.
  fresh = .false.
  order = "file"
  h = 1.0e-4_dp
  if (command_argument_count() >= 1) then
    call get_command_argument(1, arg)
    select case (trim(arg))
    case ("fd")
      doFd = .true.
      if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *) h
      end if
    case ("series")
      doSeries = .true.
      if (command_argument_count() < 2) error stop "usage: dftbp_driver series FILE [order] [fresh]"
      call get_command_argument(2, seriesFile)
      do iAt = 3, command_argument_count()
        call get_command_argument(iAt, arg)
        select case (trim(arg))
        case ("file", "byqt", "shuffle")
          order = trim(arg)
        case ("fresh")
          fresh = .true.
        case default
          error stop "series options: file|byqt|shuffle, fresh"
        end select
      end do
    case default
      error stop "usage: dftbp_driver [fd [h] | series FILE [file|byqt|shuffle] [fresh]]"
    end select
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
  call setupInstance(.true.)

  if (doSeries) then
    call runSeries(trim(seriesFile), trim(order), fresh)
    call TDftbPlus_destruct(dftbp)
    close(logUnit)
    stop
  end if

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

  !> Creates the DFTB+ instance: settings from IN_dftbplus.hsd, geometry
  !> (types, INxyz coordinates) added in code. Optionally writes the full
  !> input to dftbp_input.hsd.
  subroutine setupInstance(writeInput)
    logical, intent(in) :: writeInput
    type(TDftbPlusInput), allocatable :: input
    type(fnode), pointer :: pRoot, pGeo
    integer :: u
    allocate(input)
    call TDftbPlus_init(dftbp, outputUnit=logUnit)
    call dftbp%getInputFromFile("IN_dftbplus.hsd", input)
    call input%getRootNode(pRoot)
    call setChild(pRoot, "Geometry", pGeo)
    call setChildValue(pGeo, "Periodic", .false.)
    call setChildValue(pGeo, "TypeNames", typeNames)
    call setChildValue(pGeo, "TypesAndCoordinates", reshape(species, [1, nAtom]), coords0)
    if (writeInput) then
      open(newunit=u, file="dftbp_input.hsd", action="write", status="replace")
      call dumpHsd(input%hsdTree, u)
      close(u)
    end if
    call dftbp%setupCalculator(input)
    if (dftbp%nrOfAtoms() /= nAtom) error stop "DFTB+ atom count differs from INxyz"
  end subroutine setupInstance


  !> Energy and gradient for every geometry of a qt.dat file, in the given
  !> order, on one instance (or a fresh instance per geometry). Results are
  !> written in file order, with the call number and the time per call.
  subroutine runSeries(fname, ord, newEach)
    character(*), intent(in) :: fname, ord
    logical, intent(in) :: newEach
    real(dp), allocatable :: geo(:,:,:), eSer(:), gSer(:,:,:), tSer(:)
    integer, allocatable :: perm(:), callNo(:)
    integer :: u, ios, nLine, nGeo, nQt, iGeo, iCall, j, k, tmp
    integer(8) :: c0, c1, rate
    integer(8) :: seed
    character(len=1) :: c
    real(dp) :: w
    ! count lines, then read positions (odd lines; even lines are momenta)
    open(newunit=u, file=fname, action="read", status="old")
    nLine = 0
    do
      read(u, "(A1)", iostat=ios) c
      if (ios /= 0) exit
      nLine = nLine + 1
    end do
    if (mod(nLine, 2) /= 0) error stop "series file: odd number of lines"
    nGeo = nLine / 2
    allocate(geo(3, nAtom, nGeo), eSer(nGeo), gSer(3, nAtom, nGeo), tSer(nGeo))
    allocate(perm(nGeo), callNo(nGeo))
    rewind(u)
    do iGeo = 1, nGeo
      read(u, *) geo(:,:,iGeo), w
      read(u, "(A1)") c
    end do
    close(u)
    ! QTs per step: number of geometries before the first one repeats is not
    ! known from the file, so take it from INdyn line 9 if present
    nQt = qtsPerStep(nGeo)
    ! order of the calls
    do iGeo = 1, nGeo
      perm(iGeo) = iGeo
    end do
    select case (ord)
    case ("byqt")
      iCall = 0
      do j = 1, nQt
        do k = j, nGeo, nQt
          iCall = iCall + 1
          perm(iCall) = k
        end do
      end do
    case ("shuffle")
      seed = 12345_8
      do j = nGeo, 2, -1
        seed = mod(seed * 1103515245_8 + 12345_8, 2147483648_8)
        k = 1 + int(mod(seed, int(j, 8)))
        tmp = perm(j)
        perm(j) = perm(k)
        perm(k) = tmp
      end do
    end select
    write(*, "(A,A,A,I0,A,I0,A,A,A,L1)") "Series: ", fname, ", ", nGeo, " geometries, ", nQt,&
        & " QTs per step, order ", ord, ", fresh instance per geometry: ", newEach
    call system_clock(count_rate=rate)
    do iCall = 1, nGeo
      iGeo = perm(iCall)
      callNo(iGeo) = iCall
      call system_clock(c0)
      if (newEach .and. iCall > 1) then
        call TDftbPlus_destruct(dftbp)
        call setupInstance(.false.)
      end if
      call calc(geo(:,:,iGeo), eSer(iGeo), gSer(:,:,iGeo))
      call system_clock(c1)
      tSer(iGeo) = real(c1 - c0, dp) / real(rate, dp)
    end do
    open(newunit=u, file="series_result.dat", action="write", status="replace")
    write(u, "(A)") "# index step qt call energy(Ha) time(s)   (index = order in the file)"
    do iGeo = 1, nGeo
      write(u, "(4I6,ES24.16,ES12.4)") iGeo, (iGeo - 1) / nQt + 1, mod(iGeo - 1, nQt) + 1,&
          & callNo(iGeo), eSer(iGeo), tSer(iGeo)
    end do
    close(u)
    open(newunit=u, file="series_grad.dat", action="write", status="replace")
    do iGeo = 1, nGeo
      write(u, "(*(ES24.16))") gSer(:,:,iGeo)
    end do
    close(u)
    write(*, "(A,F10.3,A,F10.5,A)") "Total time ", sum(tSer), " s, mean ", sum(tSer) / nGeo,&
        & " s per geometry"
  end subroutine runSeries


  !> Number of QTs per large step: INdyn line 9 (first value), else nGeo
  integer function qtsPerStep(nGeo)
    integer, intent(in) :: nGeo
    integer :: u, i, n
    logical :: ex
    character(len=256) :: line
    qtsPerStep = nGeo
    inquire(file="INdyn", exist=ex)
    if (.not. ex) return
    open(newunit=u, file="INdyn", action="read", status="old")
    do i = 1, 9
      read(u, "(A)") line
    end do
    close(u)
    read(line, *) n
    if (n > 0 .and. mod(nGeo, n) == 0) qtsPerStep = n
  end function qtsPerStep


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


  !> Reads a qtes-dftb INxyz file: 6 parameter lines, then the atoms
  !> (blank lines before them are skipped). The geometry ends at the first
  !> blank line or the end of the file; anything after it is ignored, as in
  !> qtes.x. If INdyn exists, its line 1 (nat) must match the atom count.
  !> Coordinates in Angstrom.
  subroutine readInxyz(fname, sym, xyz)
    character(*), intent(in) :: fname
    character(len=2), allocatable, intent(out) :: sym(:)
    real(dp), allocatable, intent(out) :: xyz(:,:)
    integer, parameter :: maxAtom = 10000
    character(len=2) :: s(maxAtom)
    real(dp) :: r(3, maxAtom)
    character(len=256) :: line
    integer :: u, ios, i, n, nat
    logical :: hasIndyn
    open(newunit=u, file=fname, action="read", status="old")
    do i = 1, 6
      read(u, "(A)") line
    end do
    n = 0
    do
      read(u, "(A)", iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) then
        if (n > 0) exit
        cycle
      end if
      if (n == maxAtom) error stop "too many atoms in INxyz"
      n = n + 1
      read(line, *) s(n), r(:, n)
    end do
    close(u)
    if (n == 0) error stop "no atoms found in INxyz"
    inquire(file="INdyn", exist=hasIndyn)
    if (hasIndyn) then
      open(newunit=u, file="INdyn", action="read", status="old")
      read(u, *) nat
      close(u)
      if (nat /= n) then
        write(*, "(A,I0,A,I0)") "INxyz has ", n, " atoms, INdyn line 1 says ", nat
        error stop "atom count mismatch"
      end if
    end if
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
