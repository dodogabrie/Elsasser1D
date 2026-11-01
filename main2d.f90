! 2D incompressible MHD turbulence in Elsasser form (extension of main.f90).
! Evolved fields: Elsasser vorticities wp = w + j (for Z+) and wm = w - j (for Z-),
! with Z+- = zhat x grad(zeta+-), Laplacian(zeta+-) = w+-. Compact derivatives and
! Lele filter are applied direction by direction; same RK scheme as the 1D code.
! Initial condition: Orszag-Tang vortex. Usage: ./run2d.out  (reads input2d.dat)
module state2d
   real*8, allocatable :: wp(:,:), wm(:,:), ap(:,:), am(:,:)
   real*8 :: nu2, dt2, dx2
   integer :: n2
end module state2d

subroutine rhs2d()
   ! on exit ap, am hold d(w+-)/dt evaluated at (ap, am)
   use state2d
   use derivative2d
   use poisson2d
   implicit none
   real*8, dimension(n2,n2) :: zp, zm, phi, psi, om, cj, ax, ay, bx, by, jac, lap, t2, rp, rm
   call inv_laplacian(ap, zp, n2); call inv_laplacian(am, zm, n2)
   phi = 0.5d0*(zp+zm); psi = 0.5d0*(zp-zm); om = 0.5d0*(ap+am); cj = 0.5d0*(ap-am)
   ! -J(phi, om) + J(psi, cj),  J(a,b) = a_x b_y - a_y b_x
   call d_dir(phi,ax,n2,dx2,1,1); call d_dir(phi,ay,n2,dx2,2,1)
   call d_dir(om,bx,n2,dx2,1,1);  call d_dir(om,by,n2,dx2,2,1)
   rp = -(ax*by - ay*bx)
   call d_dir(psi,ax,n2,dx2,1,1); call d_dir(psi,ay,n2,dx2,2,1)
   call d_dir(cj,bx,n2,dx2,1,1);  call d_dir(cj,by,n2,dx2,2,1)
   rp = rp + (ax*by - ay*bx)
   rm = rp
   ! +- Laplacian J(phi, psi) from dj/dt
   call d_dir(phi,ax,n2,dx2,1,1); call d_dir(phi,ay,n2,dx2,2,1)
   call d_dir(psi,bx,n2,dx2,1,1); call d_dir(psi,by,n2,dx2,2,1)
   jac = ax*by - ay*bx
   call d_dir(jac,lap,n2,dx2,1,2); call d_dir(jac,t2,n2,dx2,2,2); lap = lap + t2
   rp = rp - lap; rm = rm + lap
   call d_dir(ap,lap,n2,dx2,1,2); call d_dir(ap,t2,n2,dx2,2,2); rp = rp + nu2*(lap+t2)
   call d_dir(am,lap,n2,dx2,1,2); call d_dir(am,t2,n2,dx2,2,2); rm = rm + nu2*(lap+t2)
   ap = rp; am = rm
end subroutine rhs2d

program main2d
   use state2d
   use derivative2d
   use derivative, only: diagonal
   implicit none
   integer :: nsteps, nout, i, j, k, s, fr
   real*8 :: pi, x, y, tend
   character(32) :: nm
   pi = 4.d0*datan(1.d0)
   open(1, file='input2d.dat', status='old')
   read(1,*) n2; read(1,*) dt2; read(1,*) tend; read(1,*) nu2; read(1,*) nout
   close(1)
   dx2 = 2*pi/n2; nsteps = nint(tend/dt2)
   allocate(wp(n2,n2), wm(n2,n2), ap(n2,n2), am(n2,n2))
   call init_deriv2d(n2)
   ! Orszag-Tang: phi = -cos x - cos y, psi = -cos y - cos(2x)/2
   do j = 1, n2; do i = 1, n2
      x = (i-1)*dx2; y = (j-1)*dx2
      wp(i,j) = (cos(x) + cos(y)) + (cos(y) + 2*cos(2*x))   ! w + j
      wm(i,j) = (cos(x) + cos(y)) - (cos(y) + 2*cos(2*x))   ! w - j
   enddo; enddo
   fr = 0; call save2d(fr)
   do s = 1, nsteps
      ! same nested RK of order 4 as RK_nl
      ap = wp; am = wm
      do k = 4, 1, -1
         call rhs2d()
         ap = wp + dt2/k*ap; am = wm + dt2/k*am
      enddo
      wp = ap; wm = am
      call filter2d(wp, n2); call filter2d(wm, n2)
      if (mod(s, nout) == 0) then
         fr = fr + 1; call save2d(fr); write(*,'(A,I6,A,I6)') 'step ', s, ' / ', nsteps
         if (maxval(abs(wp)) > 1d4) stop 'overflow'
      endif
   enddo
contains
   subroutine save2d(num)
      integer :: num
      write(nm, '(A,I5.5,A)') 'data2d/w', num, '.bin'
      open(2, file=trim(nm), form='unformatted', access='stream', status='replace')
      write(2) real(wp, 4), real(wm, 4)
      close(2)
   end subroutine save2d
end program main2d
