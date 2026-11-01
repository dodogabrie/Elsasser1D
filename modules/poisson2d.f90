! Spectral Poisson solver on the periodic square [0,2pi)^2 (radix-2 FFT, n a power of 2).
! ponytail: the only non-compact piece; a compact/multigrid Poisson solver would replace it.
module poisson2d
   implicit none
contains
   subroutine fft(a, n, sgn)
      integer :: n, sgn, i, j, k, m, h
      complex*16 :: a(n), t, w
      real*8 :: pi
      pi = 4.d0*datan(1.d0)
      j = 1
      do i = 1, n-1
         if (i < j) then; t = a(i); a(i) = a(j); a(j) = t; endif
         m = n/2
         do while (m >= 1 .and. j > m); j = j - m; m = m/2; enddo
         j = j + m
      enddo
      h = 1
      do while (h < n)
         do k = 0, h-1
            w = exp(cmplx(0.d0, sgn*pi*k/h, kind=8))
            do i = k+1, n, 2*h
               t = w*a(i+h); a(i+h) = a(i) - t; a(i) = a(i) + t
            enddo
         enddo
         h = 2*h
      enddo
   end subroutine fft

   subroutine inv_laplacian(w, z, n)
      ! z = Laplacian^{-1} w (zero mean)
      integer :: n, i, j, kx, ky
      real*8 :: w(n,n), z(n,n)
      complex*16 :: a(n,n)
      a = cmplx(w, 0.d0, kind=8)
      do j = 1, n; call fft(a(:,j), n, -1); enddo
      do i = 1, n; call fft(a(i,:), n, -1); enddo
      do j = 1, n
         ky = j-1; if (ky > n/2) ky = ky - n
         do i = 1, n
            kx = i-1; if (kx > n/2) kx = kx - n
            if (kx == 0 .and. ky == 0) then; a(i,j) = 0.d0
            else; a(i,j) = -a(i,j)/(kx*kx + ky*ky)/(n*n); endif
         enddo
      enddo
      do j = 1, n; call fft(a(:,j), n, 1); enddo
      do i = 1, n; call fft(a(i,:), n, 1); enddo
      z = real(a)
   end subroutine inv_laplacian
end module poisson2d
