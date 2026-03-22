      subroutine lncov(scr,sdr,shob,rcd,rch,rdh,dsig,wr0,alpha,
     *                 r95,pod,d,az,iflg,ierr)
c
c  compute probability of damage using a full 3x3 delivery error covariance.
c
c  the covariance represents errors in the weapon delivery frame:
c    axis 1: cross-range (perpendicular to flight path in horizontal plane)
c    axis 2: down-range  (along flight path in horizontal plane)
c    axis 3: hob/altitude (vertical)
c
c  the hob error is linearized: wr_eff = wr0 + alpha * delta_hob
c  where alpha = d(wr)/d(hob) is the weapon radius sensitivity to hob,
c  computed externally via two wrcalc calls.
c
c  the horizontal integral is solved via 5-point gauss-hermite quadrature
c  in the eigenbasis of the conditional horizontal covariance.  the outer
c  hob integral also uses 5-point gauss-hermite.
c
c  inputs:
c    scr    cross-range delivery 1-sigma (feet)
c    sdr    down-range delivery 1-sigma (feet)
c    shob   hob delivery 1-sigma (feet)
c    rcd    cross-range / down-range correlation  (-1 < rcd < 1)
c    rch    cross-range / hob correlation         (-1 < rch < 1)
c    rdh    down-range  / hob correlation         (-1 < rdh < 1)
c    dsig   damage sigma (vulnerability dispersion, dimensionless)
c    wr0    weapon radius at nominal hob (feet)
c    alpha  d(wr)/d(hob): weapon radius sensitivity to hob (ft/ft)
c    r95    95% target radius (nmi); circular target size descriptor
c    d      distance from dgz to target (nmi)  [in/out for iflg=6]
c    az     azimuth from dgz to target (degrees, 0=north, 90=east)
c    iflg   1: clamp pod <= 0.99
c           2: clamp pod <= 0.999
c           6: compute d such that pod equals input pod value
c    ierr   error status (in)
c
c  outputs:
c    pod    probability of damage
c    d      (iflg=6) offset distance (nmi) at which pod equals input pod
c    ierr   error status (out): 1 = desired pod unachievable
c
      implicit double precision (a-h,o-z)
      include "const.h"
      include "files.h"

      logical cross

c  5-point gauss-hermite nodes and weights for integral(f(t)*exp(-t^2)dt)
c  sum of weights = sqrt(pi); divide by pi for 2d, sqrt(pi) for each 1d axis
c
      dimension ght(5), ghw(5)
      data ght / -2.0201828704562856d0, -0.9585724646137882d0,
     *            0.0d0,
     *            0.9585724646137882d0,  2.0201828704562856d0 /
      data ghw /  0.0199532420590459d0,  0.3936193231522412d0,
     *            0.9454915028125263d0,
     *            0.3936193231522412d0,  0.0199532420590459d0 /

      double precision pvuln
      external pvuln

      ierr = 0

      if (wr0 .le. 0.001d0) then
         pod = 0.0d0
         return
      endif

c  iflg=6: iterative bisection to find d that achieves input pod value.
c  save target pod, initialise bisection, then fall through to forward calc.
c
      itch   = 0
      acc    = 0.001d0
      cross  = .false.
      dd     = 0.0d0
      if (iflg .eq. 6) then
         pod_tgt = pod
         acc     = 0.001d0
         cross   = .false.
         dd      = wr0
         d       = wr0 / cnm2ft
         itch    = 0
      endif

c  ---- forward computation of pov from current d ----
c
 5    continue

c  convert units
c
      d_ft = d * cnm2ft
      rr5  = r95 * cnm2ft

c  target position in weapon delivery frame (cross-range, down-range)
c
      az_rad = az * (3.14159265358979d0 / 180.0d0)
      tx = d_ft * sin(az_rad)
      ty = d_ft * cos(az_rad)

c  add circular target-size contribution (r95) to horizontal covariance
c  diagonal, using the same quadratic combination as lncalc:
c    adcep^2 = cep^2 + 0.231 * rr5^2
c  here 0.231 encodes the ratio of cep-to-sigma vs r95-to-sigma factors.
c
      r95v = 0.231d0 * rr5 * rr5

c  build conditional horizontal covariance via schur complement w.r.t. hob:
c    sigma_h|z = sigma_h - b*b^t / shob^2
c  where b = [cov(cr,hob), cov(dr,hob)]^t
c    = [rch*scr*shob, rdh*sdr*shob]^t
c
      if (shob .gt. 0.0d0) then
         cch = rch * scr * shob
         cdh = rdh * sdr * shob
         s2h = shob * shob
         h11 = scr*scr  - cch*cch/s2h  + r95v
         h12 = rcd*scr*sdr - cch*cdh/s2h
         h22 = sdr*sdr  - cdh*cdh/s2h  + r95v
c  conditional mean slope: e[x_cr | delta_z] = (cch/s2h)*delta_z
         slope_cr = cch / s2h
         slope_dr = cdh / s2h
         sq2hob   = sqrt(2.0d0) * shob
         nhob     = 5
      else
         h11      = scr*scr + r95v
         h12      = rcd*scr*sdr
         h22      = sdr*sdr + r95v
         slope_cr = 0.0d0
         slope_dr = 0.0d0
         sq2hob   = 0.0d0
         nhob     = 1
      endif

c  eigendecompose conditional horizontal covariance
c  lam1 >= lam2 >= 0, theta = rotation angle of first eigenvector
c
      call eigen2x2(h11, h12, h22, lam1, lam2, theta)

      sq2l1 = sqrt(2.0d0 * lam1)
      sq2l2 = sqrt(2.0d0 * lam2)
      costh = cos(theta)
      sinth = sin(theta)

c  1/pi and 1/sqrt(pi) normalisation factors for gauss-hermite
c
      piconst = 3.14159265358979d0
      piinv   = 1.0d0 / piconst
      sqpiinv = 1.0d0 / sqrt(piconst)

      psum = 0.0d0

c  outer loop: hob quadrature (or single pass when shob=0)
c
      do 100 n = 1, nhob

         if (nhob .gt. 1) then
            tn    = ght(n)
            wn    = ghw(n)
            dz_n  = sq2hob * tn
            wr_n  = wr0 + alpha * dz_n
         else
            wn    = 1.0d0
            dz_n  = 0.0d0
            wr_n  = wr0
         endif

         if (wr_n .le. 0.001d0) goto 100

c  conditional target offset (shift mean of horizontal error given hob)
c
         tx_n = tx - slope_cr * dz_n
         ty_n = ty - slope_dr * dz_n

c  rotate target offset into eigenbasis of conditional horizontal covariance
c  eigenvector 1 = (cos(theta), sin(theta))
c
         d1_n =  tx_n * costh + ty_n * sinth
         d2_n = -tx_n * sinth + ty_n * costh

c  inner 2d gauss-hermite over (u1, u2) in eigenbasis
c  u1 ~ N(0, lam1), u2 ~ N(0, lam2), independent
c
         p2d = 0.0d0

         do 90 k = 1, 5
            u1   = sq2l1 * ght(k)
            du1  = u1 - d1_n
            do 80 j = 1, 5
               u2   = sq2l2 * ght(j)
               du2  = u2 - d2_n
               dist = sqrt(du1*du1 + du2*du2)
               pv   = pvuln(dist, wr_n, dsig)
               p2d  = p2d + ghw(k) * ghw(j) * pv
 80         continue
 90      continue

c  p2d = integral(f*exp(-t1^2)*exp(-t2^2) dt1 dt2)
c  divide by pi to convert to 2d gaussian expectation
c
         psum = psum + wn * p2d * piinv

 100  continue

c  if hob was integrated, divide by sqrt(pi) for 1d gaussian normalisation
c
      if (nhob .gt. 1) then
         pov = psum * sqpiinv
      else
         pov = psum
      endif

      if (pov .lt. 0.0d0) pov = 0.0d0

c  ---- end of forward computation ----
c
      if (iflg .ne. 6) then
         if (iflg .eq. 1 .and. pov .gt. 0.99d0)  pov = 0.99d0
         if (iflg .eq. 2 .and. pov .gt. 0.999d0) pov = 0.999d0
         pod = pov
         return
      endif

c  iflg=6 bisection: iterate until d is found
c
      itch = itch + 1
      if (itch .gt. 200) then
         pod = pov
         return
      endif

      pda = abs(pod_tgt - pov)
      if (pda .lt. acc) then
         pod = pod_tgt
         return
      endif

      if (pod_tgt .gt. pov) then
         if (cross) dd = dd / 2.0d0
         d = d + dd / cnm2ft
      else
         cross = .true.
         dd = dd / 2.0d0
         d = d - dd / cnm2ft
      endif

      if (d .lt. 0.0d0) d = 0.0d0

      goto 5

      end
