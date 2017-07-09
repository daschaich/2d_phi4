! ------------------------------------------------------------------
! Corr/twodfft.f90
! Performs two-dimensional FFT given real two-dimensional array
! Based on code by Claudio Rebbi
! Last modified 2 February 2009 by David Schaich
! ------------------------------------------------------------------



! ------------------------------------------------------------------
subroutine twodfft(lattice, pow, mom)
  implicit none
  real(8), parameter :: twopi = 6.283185307179586d0
  complex(8), parameter :: iu = (0d0, 1d0)   ! Imaginary unit

  integer :: pow
  integer :: length
  ! Momentum-space correlation functions (real and complex)
  real(8) :: mom(0:2**pow - 1)
  complex(8) :: momC(0:(2**pow / 2) - 1)

  ! Lattices, one real the other complex
  real(8) :: lattice(0:2**(2 * pow) -1)
  complex(8) :: fc(0:2**pow - 1, 0:2**pow - 1)

  integer :: arg
  complex(8) :: zSum
  integer :: i, j, k, kx, ky, z, a

  ! Copy input real lattice to two-dimensional complex lattice
  length = 2**pow
  do i = 0, length - 1
    do j = 0, length - 1
      fc(i, j) = lattice(i * length + j)
    end do
  end do

  ! Calculate the Fourier transform by the FFT algorithm
  do i = 0, length - 1
    ! Overrides array with Fourier transform
    call fftr(fc(0, i), pow, 1)
  end do
  fc = transpose(fc)
  do i = 0, length - 1
    call fftr(fc(0, i), pow, 1)
  end do
  fc = transpose(fc)

  ! Now convolute Fourier transformed lattice to calculate correlation
  ! functions in momentum space
  do k = 0, (length / 2) - 1
    momC(k) = 0
    do kx = 0, length - 1
      do ky = 0, length - 1
        zSum = (0d0, 0d0)

        do z = 0, length - 1
          do a = 0, z
            arg = (ky * (a - z) - kx * a - k * z)
            if (z == 0) then
              zSum = zSum + exp(twopi * arg * iu / length) / 4
            else
              zSum = zSum + exp(twopi * arg * iu / length) / (4 * z)
            end if
          end do
        end do

        do z = 0, length - 1
          do a = 0, z - 1
            arg = (ky * (a - z) - kx * a - k * z)
            if (z == 0) then
              zSum = zSum + exp(twopi * arg * iu / length) / 4
            else
              zSum = zSum + exp(twopi * arg * iu / length) / (4 * z)
            end if
          end do
        end do

        do z = 0, length - 1
          do a = 1, z
            arg = (ky * (a - z) - kx * a - k * z)
            if (z == 0) then
              zSum = zSum + exp(twopi * arg * iu / length) / 4
            else
              zSum = zSum + exp(twopi * arg * iu / length) / (4 * z)
            end if
          end do
        end do

        do z = 0, length - 1
          do a = 1, z - 1
            arg = (ky * (a - z) - kx * a - k * z)
            if (z == 0) then
              zSum = zSum + exp(twopi * arg * iu / length) / 4
            else
              zSum = zSum + exp(twopi * arg * iu / length) / (4 * z)
            end if
          end do
        end do

        momC(k) = momC(k) + zSum * fc(kx, ky) &
                                 * fc(length - kx, length - ky)
      end do
    end do
    momC(k) = momC(k) / (length**2 * sqrt(real(length)))
  end do

  do i = 0, length / 2 - 1
    mom(2 * i) = real(momC(i))
    mom(2 * i + 1) = aimag(momC(i))
  end do
end subroutine twodfft
! ------------------------------------------------------------------



! ------------------------------------------------------------------
subroutine fftr(a, k, sg)
  implicit none
  complex(8), parameter :: iu = (0d0, 1d0)  ! Imaginary unit

  integer :: k
  integer :: sg           ! Sign, 1 or -1
  complex(8) :: a(0:2**k - 1)
  complex(8) :: aux
  complex(8) :: ec, fc     ! Roots of 1 (exponential factors)
  real(8) :: norm
  integer :: n, i, j, ka, m1, m2, m

  ! Rearrange the elements of a in bit reversed order and normalize:
  n = ishft(1, k)
  norm = 1 / sqrt(real(n, 8))

  do i = 0, n - 1
    m1 = 1    ! "Mask"
    j = 0
    do ka = 1, k
      ! m2 is the bit of i in position ka and it is put into j
      ! in position k - ka + 1
      m2 = iand(i, m1)
      j = ior(j, ishft(m2, k - 2 * ka + 1))
      m1 = ishft(m1, 1)
    end do
    ! If the two indices are different, swap elements a(i) and a(j)
    ! At the same time, normalize all elements
    if (j == i) then
      a(i) = a(i) * norm
    else if (j > i) then
      aux  = a(i) * norm
      a(i) = a(j) * norm
      a(j) = aux
    endif
  end do

  ! Implement the FFT:
  m = 1   ! m will take values 1, 2, 4, 8, ... n/2
  fc = 1d0
  ec = 1d0

10  continue
  do i = 0, m - 1
    if (i == 0) then
      do j = 0, n - 1, 2 * m
        aux = a(j + m)
        a(j + m) = a(j) - aux
        a(j) = a(j) + aux
      end do
    else
      do j = 0, n - 1, 2 * m
        aux = a(i + j + m) * fc
        a(i + j + m) = a(i + j) - aux
        a(i + j) = a(i + j) + aux
      end do
      fc = fc * ec
    endif
  end do

  m = 2 * m           ! Double m after every iteration
  if(m == n) return   ! This is how we get out of subroutine

  ! ec equals the 2 * m primitive root of 1, or its conjugate if
  ! sg = -1: ec = exp(sg * iu * pi / m), ec**m = -1, ec**(2 * m) = 1
  if(m == 2) then
    ec = iu * sg
  else
    ec = sqrt(ec)
  endif
  fc = ec

  goto 10
end subroutine fftr
! ------------------------------------------------------------------
