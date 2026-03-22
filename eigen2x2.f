      subroutine eigen2x2(a11, a12, a22, lam1, lam2, theta)
c
c  compute eigenvalues and principal axis rotation angle for a
c  2x2 symmetric positive semi-definite matrix:
c
c    A = | a11  a12 |
c        | a12  a22 |
c
c  inputs:
c    a11, a12, a22   upper triangle of the symmetric matrix
c
c  outputs:
c    lam1   larger eigenvalue  (>= 0)
c    lam2   smaller eigenvalue (>= 0)
c    theta  rotation angle (radians) such that the eigenvector for lam1
c           is (cos(theta), sin(theta)) in the original (cr,dr) frame
c
      implicit double precision (a-h,o-z)

      trace = a11 + a22
      disc  = sqrt((a11 - a22)**2 + 4.0d0 * a12 * a12)

      lam1 = (trace + disc) / 2.0d0
      lam2 = (trace - disc) / 2.0d0

c  clamp small negatives from floating point noise
c
      if (lam1 .lt. 0.0d0) lam1 = 0.0d0
      if (lam2 .lt. 0.0d0) lam2 = 0.0d0

c  rotation angle of first eigenvector
c  if matrix is already diagonal, theta = 0
c
      if (abs(a11 - a22) .lt. 1.0d-30 .and.
     *    abs(a12)        .lt. 1.0d-30) then
         theta = 0.0d0
      else
         theta = 0.5d0 * atan2(2.0d0 * a12, a11 - a22)
      endif

      return
      end
