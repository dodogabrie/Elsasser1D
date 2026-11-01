! 2D periodic extension of the compact derivatives: same stencils as DFC1/DFC2
! (alpha=1/4 first derivative, 1/10 second derivative), applied direction by
! direction. The periodic tridiagonal systems reuse `tridiag` from module
! derivative, closed with a Woodbury correction for the two corner entries.
module derivative2d
   use derivative, only: tridiag, diagonal, b_mod
   implicit none
   real*8, dimension(:,:), allocatable :: q(:,:,:)
   logical :: ready(3) = .false.
   real*8 :: sm(3,2,2)
contains
   subroutine init_deriv2d(n)
      integer :: n
      allocate(diagonal(n), b_mod(n), q(n,2,3))
   end subroutine init_deriv2d

   subroutine cyc(kind, l, d, b, sol, m)
      ! Solve the periodic system (l, d, l) x = b of size m. kind = 1,2,3 caches the setup.
      integer :: kind, m
      real*8 :: l, d, b(m), sol(m), y(m), e(m), det, w1, w2, a11, a12, a21, a22
      if (.not. ready(kind)) then
         e = 0.d0; e(1) = 1.d0; call tridiag(l, d, l, e, q(1:m,1,kind), m)
         e = 0.d0; e(m) = 1.d0; call tridiag(l, d, l, e, q(1:m,2,kind), m)
         a11 = 1.d0 + l*q(m,1,kind); a12 = l*q(m,2,kind)
         a21 = l*q(1,1,kind);        a22 = 1.d0 + l*q(1,2,kind)
         det = a11*a22 - a12*a21
         sm(kind,1,1) = a22/det; sm(kind,1,2) = -a12/det
         sm(kind,2,1) = -a21/det; sm(kind,2,2) = a11/det
         ready(kind) = .true.
      endif
      call tridiag(l, d, l, b, y, m)
      w1 = l*y(m); w2 = l*y(1)
      sol = y - q(1:m,1,kind)*(sm(kind,1,1)*w1 + sm(kind,1,2)*w2) &
              - q(1:m,2,kind)*(sm(kind,2,1)*w1 + sm(kind,2,2)*w2)
   end subroutine cyc

   subroutine line1(f, der, n, dx)
      integer :: n, i, ip, im
      real*8 :: f(n), der(n), dx, r(n)
      do i = 1, n
         ip = modulo(i, n) + 1; im = modulo(i-2, n) + 1
         r(i) = 0.75d0/dx*(f(ip) - f(im))
      enddo
      call cyc(1, 0.25d0, 1.d0, r, der, n)
   end subroutine line1

   subroutine line2(f, der, n, dx)
      integer :: n, i, ip, im
      real*8 :: f(n), der(n), dx, r(n)
      do i = 1, n
         ip = modulo(i, n) + 1; im = modulo(i-2, n) + 1
         r(i) = 1.2d0/(dx*dx)*(f(ip) - 2*f(i) + f(im))
      enddo
      call cyc(2, 0.1d0, 1.d0, r, der, n)
   end subroutine line2

   subroutine line_filter(f, n)
      ! Lele filter (same coefficients as filter6, alpha=0, beta=3/10), periodic.
      ! With alpha=0 the system couples only i+-2: even and odd points decouple.
      integer :: n, i, p, h
      real*8 :: f(n), r(n), s(n), rr(n/2), ss(n/2)
      do i = 1, n
         r(i) = 0.5d0*f(i) + 0.375d0*(f(per(i+1,n)) + f(per(i-1,n))) &
              + 0.15d0*(f(per(i+2,n)) + f(per(i-2,n))) + 0.025d0*(f(per(i+3,n)) + f(per(i-3,n)))
      enddo
      h = n/2
      do p = 1, 2
         rr = r(p:n:2)
         call cyc(3, 0.3d0, 1.d0, rr, ss, h)
         s(p:n:2) = ss
      enddo
      f = s
   end subroutine line_filter

   integer function per(i, n)
      integer :: i, n
      per = modulo(i-1, n) + 1
   end function per

   subroutine d_dir(f, der, n, dx, dir, ord)
      ! dir 1 = x (first index), 2 = y; ord = 1 or 2
      integer :: n, dir, ord, i
      real*8 :: f(n,n), der(n,n), dx, t(n)
      do i = 1, n
         if (dir == 1) then
            if (ord == 1) then; call line1(f(:,i), der(:,i), n, dx); else; call line2(f(:,i), der(:,i), n, dx); endif
         else
            if (ord == 1) then; call line1(f(i,:), t, n, dx); else; call line2(f(i,:), t, n, dx); endif
            der(i,:) = t
         endif
      enddo
   end subroutine d_dir

   subroutine filter2d(f, n)
      integer :: n, i
      real*8 :: f(n,n), t(n)
      do i = 1, n
         call line_filter(f(:,i), n)
      enddo
      do i = 1, n
         t = f(i,:); call line_filter(t, n); f(i,:) = t
      enddo
   end subroutine filter2d
end module derivative2d
