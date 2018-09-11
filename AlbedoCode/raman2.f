      FUNCTION RAMALB( WAVEL )
	implicit double precision (a-h,o-z)
C
C  +-------------------------------------------------------------------+
C  !  ALTERS RAYLEIGH SCATTERING ALBEDO TO ALLOW FOR RAMAN SCATTERING  !
C  !  Returned value is modified single scattering albedo              !
C  !                                                                   !
C  !     wavel -- wavelength in MICRONS                                !
C  +-------------------------------------------------------------------+
C
      PARAMETER( H = 6.6252D-27, C = 2.9978D10, HMASS = 1.6734D-24 )
      PARAMETER( RMU = .5 * HMASS, PI = 3.141592, BOHRD = 5.2917D-9 )
      PARAMETER( SHIFTV0 = 4161. )
      PARAMETER( TEMP = 6000.   )
      PARAMETER( FACIP = H * C / ( 1.D-4 * 27.2 * 1.602D-12 ) )
      PARAMETER( FACRAY = 1.D16 * BOHRD ** 3 * 128. * PI ** 5
     \                  * BOHRD ** 3 / 9. )
      PARAMETER( FACV = 2.08 / 2.38 * FACRAY / BOHRD ** 2
     \                / ( 8. * PI * PI * RMU * C * SHIFTV0 ) * H )
C
      DIMENSION GLI(5), WLI(5), GRI(5), WRI(5)
      DIMENSION ALP(7), ARP(7)
C
      DATA ALP / 6.84, 6.96, 7.33, 8.02, 9.18, 11.1, 14.5 /
      DATA ARP / 3.66, 3.71, 3.88, 4.19, 4.70, 5.52, 6.88 /
C
      DATA GLI / 1.296, .247, .297,  .157,  .003 /
      DATA WLI /  .507, .628, .733, 1.175, 2.526 /
      DATA GRI /  .913, .239, .440,  .344,  .064 /
      DATA WRI /  .537, .639, .789, 1.304, 3.263 /
C
      SHIFT( WAVEL, WSHFT ) = 1.D4 * WAVEL / ( 1.D4 + WAVEL * WSHFT )
C
      OMEGA = FACIP / WAVEL
C
C  EXTINCTION CROSS SECTION FOR UNSHIFTED COMPONENT( i.e. Rayleigh)
C
      ALPHL = 0.
      ALPHR = 0.
C
      DO 10 I = 1, 5
         ALPHL = ALPHL + GLI(I) / ( WLI(I) ** 2 - OMEGA ** 2 )
         ALPHR = ALPHR + GRI(I) / ( WRI(I) ** 2 - OMEGA ** 2 )
   10 CONTINUE
C
      ALPHA2 = (( 2. * ALPHR + ALPHL ) / 3. ) ** 2
      GAMMA2 = ( ALPHL - ALPHR ) ** 2
      QRAY = FACRAY * ( 3. * ALPHA2 + 2./3. * GAMMA2 ) / WAVEL ** 4
C
C  EXTINCTION CROSS SECTION FOR VIBRATIONALLY SHIFTED COMPONENT
C jjf changed min1 to dmin1
c	write(*,*) 'in Raman'
      IP = DMIN1( 1.* INT( OMEGA / .05 ), 5. ) + 1
c      IP = DMIN1( INT( OMEGA / .05 ), 5 ) + 1  !there was a "real" issue
      F = OMEGA / .05 - DBLE( IP - 1 )
      ALPHPL = ( 1. - F ) * ALP(IP) + F * ALP(IP+1)
      ALPHPR = ( 1. - F ) * ARP(IP) + F * ARP(IP+1)
C
      ALPHP2 = (( 2. * ALPHPR + ALPHPL ) / 3. ) ** 2
      GAMMP2 = ( ALPHPL - ALPHPR ) ** 2
      QV = FACV / SHIFT( WAVEL, -SHIFTV0 ) ** 4
     \   * ( 3. * ALPHP2 + 2./3. * GAMMP2 )
C
C  Rayleigh-Raman single scattering albedo.
C  WARNING:  will be greater than 1 for wavelengths longer than that
C            of maximum solar intensity
C

      SF = Planck( WAVEL ,temp )
      proxy=shift(wavel, shiftv0)
      RAMALB=(QRAY+QV*planck(proxy, temp)/SF)
     \       / ( QRAY + QV )
c      write (43,*) sngl(wavel),sngl(qray),sngl(qv)
c     &   ,sngl(proxy), sngl(planck(proxy, temp)/sf), sngl(ramalb)
c     write (43,*) wavel,ramalb
C
      RETURN
      END
C
      Function Planck(wavel, temp)
      h0=6.625d-27
      c=3.d10
      xk=1.38d-16
      w=wavel*1.d-4
      planck=(2*h0*c**2)/((exp(h0*c/(w*xk*temp))-1.)*w**5)
      end

      FUNCTION SMSFLUX( WAVEL )
	implicit double precision (a-h,o-z)
C
C         +--------------------------------------------------+
C         !  SMOOTHED SOLAR SPECTRUM -- FLUX VS. WAVELENGTH  !
C         !  UNITS -- ERG / SEC / CM**2 / A AT 1AU           !
C         !                                                  !
C         !     wavel -- wavelength in MICRONS               !
C         +--------------------------------------------------+
C
      DIMENSION ISUN(383), ISUN1(190), ISUN2(190), ISUN3(3)
      EQUIVALENCE ( ISUN(1)  , ISUN1 ), ( ISUN(191), ISUN2 )
      EQUIVALENCE ( ISUN(381), ISUN3 )
C
      DATA ISUN1/
     \     425,  438,  452,  466,  480,  494,  522,  533,  555,  580,
     \     592,  620,  644,  665,  668,  712,  759,  813,  859,  878,
     \     913,  943,  969,  985, 1005, 1047, 1027, 1035, 1018, 1022,
     \    1040, 1010, 1033, 1034, 1048, 1085, 1102, 1144, 1157, 1160,
     \    1171, 1211, 1249, 1223, 1207, 1199, 1169, 1174, 1173, 1146,
     \    1149, 1176, 1200, 1263, 1342, 1410, 1479, 1525, 1582, 1641,
     \    1729, 1790, 1814, 1820, 1817, 1817, 1827, 1784, 1811, 1784,
     \    1794, 1780, 1786, 1793, 1836, 1846, 1874, 1914, 1986, 2009,
     \    2061, 2069, 2097, 2125, 2128, 2135, 2148, 2138, 2135, 2133,
     \    2127, 2124, 2130, 2136, 2127, 2081, 2074, 2080, 2070, 2055,
     \    2047, 2034, 2015, 1993, 1967, 1966, 1990, 1977, 1963, 1959,
     \    1949, 1925, 1922, 1919, 1923, 1924, 1908, 1910, 1911, 1905,
     \    1918, 1921, 1926, 1925, 1926, 1919, 1913, 1911, 1912, 1908,
     \    1905, 1893, 1880, 1879, 1874, 1859, 1856, 1857, 1858, 1856,
     \    1853, 1850, 1852, 1861, 1865, 1865, 1872, 1866, 1862, 1856,
     \    1849, 1842, 1837, 1823, 1813, 1802, 1792, 1778, 1771, 1759,
     \    1743, 1733, 1728, 1721, 1717, 1706, 1695, 1686, 1678, 1672,
     \    1667, 1669, 1661, 1653, 1647, 1642, 1641, 1638, 1639, 1640,
     \    1611, 1609, 1604, 1604, 1601, 1598, 1593, 1587, 1582, 1576 /
C
      DATA ISUN2/
     \    1569, 1590, 1584, 1577, 1565, 1556, 1547, 1539, 1531, 1524,
     \    1518, 1505, 1496, 1487, 1479, 1475, 1470, 1465, 1459, 1453,
     \    1445, 1432, 1422, 1413, 1404, 1394, 1383, 1373, 1363, 1357,
     \    1347, 1340, 1336, 1332, 1327, 1323, 1320, 1317, 1313, 1309,
     \    1302, 1297, 1269, 1253, 1241, 1227, 1220, 1213, 1206, 1201,
     \    1194, 1190, 1184, 1204, 1210, 1215, 1218, 1215, 1209, 1202,
     \    1196, 1190, 1181, 1177, 1170, 1164, 1157, 1149, 1142, 1137,
     \    1133, 1126, 1118, 1112, 1105, 1099, 1092, 1086, 1081, 1075,
     \    1069, 1063, 1059, 1054, 1048, 1043, 1037, 1022, 1015,  988,
     \     981,  976,  970,  963,  958,  937,  929,  925,  922,  923,
     \     938,  935,  929,  925,  920,  912,  923,  920,  916,  913,
     \     909,  905,  901,  897,  892,  888,  887,  885,  881,  877,
     \     873,  870,  867,  864,  862,  859,  856,  854,  851,  848,
     \     846,  843,  840,  838,  835,  833,  830,  827,  825,  822,
     \     820,  817,  815,  812,  810,  808,  805,  803,  800,  798,
     \     796,  793,  791,  788,  786,  784,  781,  779,  777,  775,
     \     772,  768,  766,  763,  761,  760,  760,  759,  758,  757,
     \     756,  756,  757,  756,  755,  754,  752,  749,  747,  744,
     \     741,  739,  736,  734,  742,  744,  748,  769,  766,  763 /
C
      DATA ISUN3 / 760, 757, 753 /
C
      DATA FS, HCKTS / 1147., 2.7326 /
      DATA FN, HCKTN / 1245., 2.7669 /
      DATA SWAVS, SWAVI, SWAVN, NSUN / 0.286, 0.002, 1.050, 383 /
C
      IF( WAVEL .LT. .296 ) THEN
         EXPF = 0.
         IF( HCKTS .LT. 25. * WAVEL )
     \      EXPF = FS / ( DEXP( HCKTS / WAVEL ) - 1. )
         SMSFLUX = EXPF / WAVEL ** 5
C
      ELSE IF( WAVEL .GT. 1.04 ) THEN
         SMSFLUX = FN / ( WAVEL ** 5 * (DEXP( HCKTN / WAVEL ) - 1. ) )
C
      ELSE
         NW = INT( ( WAVEL - SWAVS ) / SWAVI ) + 1
         X1 = SWAVS + SWAVI * DBLE( NW - 1 )
         X2 = X1 + SWAVI
         F1 = .1 * DBLE( ISUN(NW) )
         F2 = .1 * DBLE( ISUN(NW+1) )
         S1 = .5 * ( F2 - .1 * DBLE( ISUN(NW-1) ) ) / SWAVI
         S2 = .5 * ( .1 * DBLE( ISUN(NW+2) ) - F1 ) / SWAVI
         call cubic1( x1, x2,f1,
     1      f2, s1, s2,wavel, smsflux, slope )
      END IF
C
      RETURN
      END
C
      SUBROUTINE Cubic1( X1, X2, F1, F2, S1, S2, XI, YOUT, SLOPE )
      implicit double precision (a-h,o-z)
C
C  +------------------------------------------------------------------+
C  !  Cubic interpolation using function and slope at two points      !
C  !                                                                  !
C  !  x1, x2 (input) -- values of independent variable at two points  !
C  !  f1, f2 (input) -- value of function at x1 and x2                !
C  !  s1, s2 (input) -- value of slope at x1 and x2                   !
C  !  xi (input) -- value of independent variable for which           !
C  !       interpolation is wanted                                    !
C  !  yout (output) -- value of function at xi                        !
C  !  slope (output) -- value of slope at xi                          !
C  !                                                                  !
C  !              ********** WARNING **********                       !
C  !  Cubic1 requires 64 bit accuracy (more than 32 bits, anyway).    !
C  !  I can supply a double-precision version.                        !
C  +------------------------------------------------------------------+
C
      DX = X2 - X1
      A = ( ( S2 + S1 ) * DX - 2.d0 * ( F2 - F1 ) ) / DX ** 3
      B = .5d0 * ( S2 - S1 ) / DX - 1.5d0 * A * ( X2 + X1 )
      C = S1 - ( 3.d0 * A * X1 + 2.d0 * B ) * X1
      D = F1 - (( A * X1 + B ) * X1 + C ) * X1
C
      YOUT = (( ( A * XI + B ) * XI + C ) * XI + D)
      SLOPE = (( 3.d0 * A * XI + 2.d0 * B ) * XI + C)
      RETURN
      END

