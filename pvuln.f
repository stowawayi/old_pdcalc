      double precision function pvuln(dist, wr, dsig)
c
c  compute probability of damage given exact impact distance from target.
c  this is the lognormal vulnerability model extracted from lncalc.f.
c
c  inputs:
c    dist   distance from weapon impact point to target (feet)
c    wr     weapon radius at nominal hob (feet)
c    dsig   damage sigma (vulnerability dispersion)
c
c  returns:
c    pvuln  probability of damage (0 to 0.999)
c
      implicit double precision (a-h,o-z)

      if (wr .le. 0.001d0) then
         pvuln = 0.0d0
         return
      endif

      if (dist .le. 0.0d0) then
         pvuln = 0.999d0
         return
      endif

c  compute beta factor
c
      ex   = 1.0d0 - dsig * dsig
      beta = sqrt(-log(ex))

      z = (1.0d0 / beta) * log(wr * ex / dist)

      if (z .gt. 3.87d0) then
         pvuln = 0.999d0
         return
      endif

      if (z .lt. (-3.87d0)) then
         pvuln = 0.0d0
         return
      endif

      zab = abs(z)
      if (zab .lt. 5.0d-7) then
         pvuln = 0.500d0
         return
      endif

c  evaluate gaussian cdf via rational polynomial approximation of erfc
c
      c = abs(z) / 1.414213562d0

      c2 = c * c
      c3 = c2 * c
      c4 = c3 * c
      c5 = c4 * c
      c6 = c5 * c

      erfu = 1.0d0 - 1.0d0 / ((1.0d0 + 0.0705230784d0 * c  +
     *                                  0.0422820123d0 * c2 +
     *                                  0.0092705272d0 * c3 +
     *                                  0.0001520143d0 * c4 +
     *                                  0.0002765672d0 * c5 +
     *                                  0.0000430638d0 * c6)**16)

      if (z .lt. 0.0d0) then
         pvuln = 0.5d0 - 0.5d0 * erfu
      else
         pvuln = 0.5d0 + 0.5d0 * erfu
      endif

      return
      end
