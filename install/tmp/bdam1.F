*Deck BDAM1
      Block Data BDAM1
C
C     Portions of this subroutine are taken from version 5 of MOPAC,
C     which is in the public domain.
C
      Implicit Real*8(A-H,O-Z)
      COMMON /NATORB/ NATORB(107)
C**********************************************************************
C
C     COMMON BLOCKS FOR AM1
C
C**********************************************************************
      COMMON /ELEMTS/ ELEMNT(107)
     $       /ALPHA / ALP(107)
     $       /MNCORE/ CORE(107)
     $       /MULTIP/ DD(107),QQ(107),AM(107),AD(107),AQ(107)
     $       /EXPONT/ ZS(107),ZP(107),ZD(107)
     $       /ONELEC/ USS(107),UPP(107),UDD(107)
     $       /BETAS / BETAS(107),BETAP(107),BETAD(107)
     $       /TWOELE/ GSS(107),GSP(107),GPP(107),GP2(107),HSP(107),
     $                GSD(107),GPD(107),GDD(107)
     $       /ATOMIC/ EISOL(107),EHEAT(107)
     $       /VSIPS / VS(107),VP(107),VD(107)
     $       /ISTOPE/ AMS(107)
     $       /IDEAS / Guess1(107,10),Guess2(107,10),Guess3(107,10)
     $       /IDEAP / GuesP1(107,10),GuesP2(107,10),GuesP3(107,10)
C**********************************************************************
C
C     COMMON BLOCKS FOR MNDO
C
C**********************************************************************
      COMMON /MNDO/  USSM(107), UPPM(107), UDDM(107), ZSM(107),
     $ZPM(107), ZDM(107), BETASM(107), BETAPM(107), BETADM(107),
     $ALPM(107), EISOLM(107), DDM(107), QQM(107), AMM(107), ADM(107),
     $AQM(107) ,GSSM(107), GSPM(107), GPPM(107), GP2M(107), HSPM(107),
     $POLVOM(107)
C**********************************************************************
C
C     COMMON BLOCKS FOR PM3
C
C**********************************************************************
      COMMON /PM3 /  USSPM3(107), UPPPM3(107), UDDPM3(107), ZSPM3(107),
     $ZPPM3(107), ZDPM3(107), BETASP(107), BETAPP(107), BETADP(107),
     $ALPPM3(107), EISOLP(107), DDPM3(107), QQPM3(107), AMPM3(107),
     $ADPM3(107), AQPM3(107) ,GSSPM3(107), GSPPM3(107), GPPPM3(107),
     $GP2PM3(107), HSPPM3(107),POLVOP(107)
      Character*80 RefMN, RefM3, RefAM, REFPM3
      COMMON /ALPTM/ ALPTM(30), EMUDTM(30)
     $       /REFS/ REFM3(107), RefMN(107), REFAM(107), REFPM3(107)
C**********************************************************************
C
C     COMMON BLOCKS FOR MINDO/3
C
C**********************************************************************
      COMMON /ONELE3 /  USS3(18),UPP3(18)
     $       /TWOEL3 /  F03(107)
     $       /ATOMI3 /  EISOL3(18),EHEAT3(18)
     $       /BETA3  /  BETA3(153)
     $       /ALPHA3 /  ALP3(153)
     $       /EXPON3 /  ZS3(18),ZP3(18)
      CHARACTER*2 ELEMNT
C**********************************************************************
C
C     Common block with parameters availability arrays.
C
C**********************************************************************
      Parameter(NKM3=10,NKM3B=40,NKMNDO=30,NKAM1=27,NKPM3=30,
     $  LMix=(NKAM1+NKPM3))
      Common /SEKeep/KM3(NKM3),KM3B(NKM3B),KMNDO(NKMNDO),KAM1(NKAM1),
     $  KPM3(NKPM3),KMix(LMix)
      Data KM3/1,5,6,7,8,9,14,15,16,17/
      Data KM3B/1,11,15,16,20,21,22,26,27,28,29,33,34,35,36,37,41,42,43,
     $  44,45,92,97,105,106,111,120,121,126,127,128,129,136,137,142,143,
     $  144,151,152,153/
      Data KMNDO/1,3,4,5,6,7,8,9,11,13,14,15,16,17,19,24,30,32,35,50,53,
     $  80,82,90,91,102,103,104,105,106/
      Data KAM1/ 1,3,4,5,6,7,8,9,11,13,14,15,16,17,19,24,30,32,35,50,53,
     $  54,80,82,90,91,102/
      Data KPM3/1,3,4,6,7,8,9,12,13,14,15,16,17,30,31,32,33,34,35,48,49,
     $  50,51,52,53,80,81,82,83,102/
C**********************************************************************
      DATA ELEMNT/' H','He',
     $ 'Li','Be',' B',' C',' N',' O',' F','Ne',
     $ 'Na','Mg','Al','Si',' P',' S','Cl','Ar',
     $ ' K','Ca','Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu',
     $ 'Zn','Ga','Ge','As','Se','Br','Kr',
     $ 'Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag',
     $ 'Cd','In','Sn','Sb','Te',' I','Xe',
     $ 'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     $ 'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir','Pt',
     $ 'Au','Hg','Tl','Pb','Bi','Po','At','Rn',
     $ 'Fr','Ra','Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','XX',
     $ 'Fm','Md','No','++',' +','--',' -','Tv'/
C
C   NATORB IS THE NUMBER OF ATOMIC ORBITALS PER ATOM.
C
C   Lithium switched from 0 to 4 orbs to go with Thiel parameters.
      DATA NATORB/2*1, 4, 7*4, 0, 7*4, 0, 4, 9*9, 7*4,
     $2*4, 9*9, 7*4, 2*2, 14*8, 9*9, 7*4, 21*0/
C**********************************************************************
C                      VALENCE SHELLS ARE DEFINED AS                  *
C  PQN   VALENCE SHELLS                                               *
C                 P-GROUP              F-GROUP    TRANSITION METALS   *
C   1       1S                                                        *
C   2       2S 2P                                                     *
C   3       3S 3P  OR  3S 3P 3D                                       *
C   4       4S 4P                                    4S 4P 3D         *
C   5       5S 5P                                    5S 5P 4D         *
C   6       6S 6P                       6S 4F        6S 6P 5D         *
C   7  NOT ASSIGNED YET  ****DO  NOT  USE****                         *
C**********************************************************************
      DATA      POLVOM(1) /0.2287D0/
      DATA      POLVOM(6) /0.2647D0/
      DATA      POLVOM(7) /0.3584D0/
      DATA      POLVOM(8) /0.2324D0/
      DATA      POLVOM(9) /0.1982D0/
      DATA      POLVOM(17)/1.3236D0/
      DATA      POLVOM(35)/2.2583D0/
      DATA      POLVOM(53)/4.0930D0/
C
C                STANDARD ATOMIC MASSES
C
      DATA  AMS /  1.00790D0,  4.00260D0,  6.94000D0,  9.01218D0,
     $10.81000D0, 12.01100D0, 14.00670D0, 15.99940D0, 18.99840D0,
     $20.17900D0, 22.98977D0, 24.30500D0, 26.98154D0, 28.08550D0,
     $30.97376D0, 32.06000D0, 35.45300D0, 39.94800D0, 39.09830D0,
     $40.08000D0, 44.95590D0, 47.90000D0, 50.94150D0, 51.99600D0,
     $54.93800D0, 55.84700D0, 58.93320D0, 58.71000D0, 63.54600D0,
     $65.38000D0, 69.73500D0, 72.59000D0, 74.92160D0, 78.96000D0,
     $79.90400D0, 83.80000D0, 85.46780D0, 87.62000D0, 88.90590D0,
     $91.22000D0, 92.90640D0, 95.94000D0, 98.90620D0, 101.0700D0,
     $102.9055D0, 106.4000D0, 107.8680D0, 112.4100D0, 114.8200D0,
     $118.6900D0, 121.7500D0, 127.6000D0, 126.9045D0, 131.3000D0,
     $132.9054D0, 137.3300D0, 15*0.000D0, 178.4900D0, 180.9479D0,
     $183.8500D0, 186.2070D0, 190.2000D0, 192.2200D0, 195.0900D0,
     $196.9665D0, 200.5900D0, 204.3700D0, 207.2000D0, 208.9804D0,
     $24*0.000D0/
C
C   CORE IS THE CHARGE ON THE ATOM AS SEEN BY THE ELECTRONS
C
      DATA CORE/1.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     $ 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     $ 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,
     $ 11.D0,12.D0,13.D0,14.D0,15.D0,16.D0,
     $ 3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     $ 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $  16*0.D0,2.D0,1.D0,-2.D0,-1.D0,0.D0/
C
C     ENTHALPIES OF FORMATION OF GASEOUS ATOMS ARE TAKEN FROM \ANNUAL
C     REPORTS,1974,71B,P 117\  THERE ARE SOME SIGNIFICANT DIFFERENCES
C     BETWEEN THE VALUES REPORTED THERE AND THE VALUES PREVIOUSLY IN
C     THE BLOCK DATA OF THIS PROGRAM.  ONLY THE THIRD  ROW ELEMENTS
C     HAVE BEEN UPDATED.
C
C ALL THE OTHER ELEMENTS ARE TAKEN FROM CRC HANDBOOK 1981-1982
      DATA EHEAT( 1) / 52.102D0/
      DATA EHEAT( 2) /  0.000D0/
C
      DATA EHEAT( 3) / 38.410D0/
      DATA EHEAT( 4) / 76.960D0/
      DATA EHEAT( 5) /135.700D0/
      DATA EHEAT( 6) /170.890D0/
      DATA EHEAT( 7) /113.000D0/
      DATA EHEAT( 8) / 59.559D0/
      DATA EHEAT( 9) / 18.890D0/
      DATA EHEAT(10) /  0.000D0/
C
      DATA EHEAT(11) / 25.850D0/
      DATA EHEAT(12) / 35.000D0/
      DATA EHEAT(13) / 79.490D0/
      DATA EHEAT(14) /108.390D0/
      DATA EHEAT(15) / 75.570D0/
      DATA EHEAT(16) / 66.400D0/
      DATA EHEAT(17) / 28.990D0/
      DATA EHEAT(18) /  0.000D0/
C
      DATA EHEAT(19) / 21.420D0/
      DATA EHEAT(20) / 42.600D0/
      DATA EHEAT(21) / 90.300D0/
      DATA EHEAT(22) /112.300D0/
      DATA EHEAT(23) /122.900D0/
      DATA EHEAT(24) / 95.000D0/
      DATA EHEAT(25) / 67.700D0/
      DATA EHEAT(26) / 99.300D0/
      DATA EHEAT(27) /102.400D0/
      DATA EHEAT(28) /102.800D0/
      DATA EHEAT(29) / 80.700D0/
      DATA EHEAT(30) / 31.170D0/
      DATA EHEAT(31) / 65.400D0/
      DATA EHEAT(32) / 89.500D0/
      DATA EHEAT(33) / 72.300D0/
      DATA EHEAT(34) / 54.300D0/
      DATA EHEAT(35) / 26.740D0/
      DATA EHEAT(36) /  0.000D0/
C
      DATA EHEAT(37) / 19.600D0/
      DATA EHEAT(38) / 39.100D0/
      DATA EHEAT(39) /101.500D0/
      DATA EHEAT(40) /145.500D0/
      DATA EHEAT(41) /172.400D0/
      DATA EHEAT(42) /157.300D0/
      DATA EHEAT(44) /155.500D0/
      DATA EHEAT(45) /133.000D0/
      DATA EHEAT(46) / 90.000D0/
      DATA EHEAT(47) / 68.100D0/
      DATA EHEAT(48) / 26.720D0/
      DATA EHEAT(49) / 58.000D0/
      DATA EHEAT(50) / 72.200D0/
      DATA EHEAT(51) / 63.200D0/
      DATA EHEAT(52) / 47.000D0/
      DATA EHEAT(53) / 25.517D0/
      DATA EHEAT(54) /  0.000D0/
C
      DATA EHEAT(55) / 18.700D0/
      DATA EHEAT(56) / 42.500D0/
      DATA EHEAT(58) /101.300D0/
      DATA EHEAT(62) / 49.400D0/
      DATA EHEAT(68) / 75.800D0/
      DATA EHEAT(70) / 36.350D0/
      DATA EHEAT(72) /148.000D0/
      DATA EHEAT(73) /186.900D0/
      DATA EHEAT(74) /203.100D0/
      DATA EHEAT(75) /185.000D0/
      DATA EHEAT(76) /188.000D0/
      DATA EHEAT(77) /160.000D0/
      DATA EHEAT(78) /135.200D0/
      DATA EHEAT(79) / 88.000D0/
      DATA EHEAT(80) / 14.690D0/
      DATA EHEAT(81) / 43.550D0/
      DATA EHEAT(82) / 46.620D0/
      DATA EHEAT(83) / 50.100D0/
      DATA EHEAT(86) /  0.000D0/
C
      DATA VS( 1) /-13.605D0/
      DATA VS( 5) /-15.16D00/
      DATA VS( 6) /-21.34D00/
      DATA VS( 7) /-27.51D00/
      DATA VS( 8) /-35.30D00/
      DATA VS( 9) /-43.70D00/
      DATA VS(14) /-17.82D00/
      DATA VS(15) /-21.10D00/
      DATA VS(16) /-23.84D00/
      DATA VS(17) /-25.26D00/
      DATA VP( 1) /  0.0D00 /
      DATA VP( 5) / -8.52D00/
      DATA VP( 6) /-11.54D00/
      DATA VP( 7) /-14.34D00/
      DATA VP( 8) /-17.91D00/
      DATA VP( 9) /-20.89D00/
      DATA VP(14) / -8.51D00/
      DATA VP(15) /-10.29D00/
      DATA VP(16) /-12.41D00/
      DATA VP(17) /-15.09D00/
C      DATA NPQ/1,1, 2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3, 4,4,4,4,4,4,4,4,
C     +4,4,4,4,4,4,4,4,4,4, 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/
C
C *** ONE CENTER REPULSION INTEGRALS
C     GSS ::= (SS,SS)
C     GPP ::= (PP,PP)
C     GSP ::= (SS,PP)
C     GP2 ::= (PP,P*P*)
C     HSP ::= (SP,SP)
C***********************************************************************
      DATA GSSM( 1) /12.848D00 /
      DATA GSSM( 4) / 9.00D00/
      DATA GSSM( 5) /10.59D00/
      DATA GSSM( 6) /12.23D00 /
      DATA GSSM( 7) /13.59D00/
      DATA GSSM( 8) /15.42D00/
      DATA GSSM( 9) /16.92D00/
      DATA GSSM(11) / 6.5D00/
      DATA GSSM(13) / 8.09D00/
      DATA GSSM(14) / 9.82D00/
      DATA GSSM(15) /11.56D00/
      DATA GSSM(16) /12.88D00/
      DATA GSSM(17) /15.03D00/
      DATA GSSM(19) /15.03D00/
      DATA GSSM(35) /15.03643948D0/
      DATA GSSM(53) /15.04044855D0/
      DATA GPPM( 4) / 6.97D00/
      DATA GPPM( 5) / 8.86D00/
      DATA GPPM( 6) /11.08D00 /
      DATA GPPM( 7) /12.98D00/
      DATA GPPM( 8) /14.52D00/
      DATA GPPM( 9) /16.71D00/
      DATA GPPM(11) / 5.0D00/
      DATA GPPM(13) / 5.98D00/
      DATA GPPM(14) / 7.31D00/
      DATA GPPM(15) / 8.64D00/
      DATA GPPM(16) / 9.90D00/
      DATA GPPM(17) /11.30D00/
      DATA GPPM(19) /11.30D00/
      DATA GPPM(35) /11.27632539D0/
      DATA GPPM(53) /11.14778369D0/
      DATA GSPM( 4) / 7.43D00/
      DATA GSPM( 5) / 9.56D00/
      DATA GSPM( 6) /11.47D00 /
      DATA GSPM( 7) /12.66D00/
      DATA GSPM( 8) /14.48D00/
      DATA GSPM( 9) /17.25D00/
      DATA GSPM(11) / 5.D00/
      DATA GSPM(13) / 6.63D00/
      DATA GSPM(14) / 8.36D00/
      DATA GSPM(15) /10.08D00/
      DATA GSPM(16) /11.26D00/
      DATA GSPM(17) /13.16D00/
      DATA GSPM(19) /13.16D00/
      DATA GSPM(35) /13.03468242D0/
      DATA GSPM(53) /13.05655798D0/
      DATA GP2M( 4) / 6.22D00/
      DATA GP2M( 5) / 7.86D00/
      DATA GP2M( 6) / 9.84D00 /
      DATA GP2M( 7) /11.59D00/
      DATA GP2M( 8) /12.98D00/
      DATA GP2M( 9) /14.91D00/
      DATA GP2M(11) / 4.70D00/
      DATA GP2M(13) / 5.40D00/
      DATA GP2M(14) / 6.54D00/
      DATA GP2M(15) / 7.68D00/
      DATA GP2M(16) / 8.83D00/
      DATA GP2M(17) / 9.97D00/
      DATA GP2M(19) / 9.97D00/
      DATA GP2M(35) / 9.85442552D0/
      DATA GP2M(53) / 9.91409071D0/
      DATA HSPM( 4) / 1.28D00/
      DATA HSPM( 5) / 1.81D00/
      DATA HSPM( 6) / 2.43D00 /
      DATA HSPM( 7) / 3.14D00/
      DATA HSPM( 8) / 3.94D00/
      DATA HSPM( 9) / 4.83D00/
      DATA HSPM(11) / 0.70D00/
      DATA HSPM(13) / 0.70D00/
      DATA HSPM(14) / 1.32D00/
      DATA HSPM(15) / 1.92D00/
      DATA HSPM(16) / 2.26D00/
      DATA HSPM(17) / 2.42D00/
      DATA HSPM(19) / 2.42D00/
      DATA HSPM(35) / 2.45586832D0/
      DATA HSPM(53) / 2.45638202D0/
C
C     THE MONOCENTRIC INTEGRALS HSP AND GSP FOR ALUMINIUM ARE ONLY
C     ESTIMATES. A VALUE OF G1 FOR AL IS NEEDED TO RESOLVE OLEARIS
C     INTEGRALS.
C
C     OPTIMIZED MNDO PARAMETERS FOR H, BE, B, C, N, O, F
C                                                     CL
C     ESTIMATED MNDO PARAMETERS FOR       AL,SI, P, S
C
C     ELEMENTS H, C, N, O WERE PARAMETERIZED BY WALTER THIEL
C     ELEMENTS B,SI,P,S   WERE      ..          MICHAEL MCKEE
C     ELEMENTS BE,F,AL,CL WERE      ..          HENRY RZEPA
C
C**********************************************************************
C
C    START OF MINDO/3 PARAMETERS
C
C**********************************************************************
      Data RefM3( 1)( 1:40)/' H: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3( 1)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3( 5)( 1:40)/' B: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3( 5)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3( 6)( 1:40)/' C: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3( 6)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3( 7)( 1:40)/' N: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3( 7)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3( 8)( 1:40)/' O: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3( 8)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3( 9)( 1:40)/' F: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3( 9)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3(14)( 1:40)/'SI: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3(14)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3(15)( 1:40)/' P: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3(15)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3(16)( 1:40)/' S: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3(16)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
      Data RefM3(17)( 1:40)/'CL: (MINDO/3): R.C.BINGHAM ET.AL., J.AM.'/
      Data RefM3(17)(41:80)/'CHEM.SOC. 97,1285,1294,1302,1307 (1975)'/
C *** F03 IS THE ONE CENTER AVERAGED REPULSION INTEGRAL FOR USE IN THE
C        TWO CENTER ELECTRONIC REPULSION INTEGRAL EVALUATION.
      Data F03              /  12.848D0, 10.0D0, 10.0D0, 0.0D0,
     $  8.958D0, 10.833D0, 12.377D0, 13.985D0, 16.250D0,
     $         10.000D0, 10.000D0, 0.000D0, 0.000D0,7.57D0 ,  9.00D0 ,
     $         10.20D0 , 11.73d0,10.0D0,35*0.D0,10.D0,53*10.D0/
C *** USS AND UPP ARE THE ONE-CENTER CORE ELECTRON ATTRACTION AND KINETI
C     ENERGY INTEGRALS FOR S AND P ELECTRONS RESPECTIVELY IN E.V.
      Data USS3             / -12.505D0, 0.000D0, 0.000D0, 0.000D0,
     $                       -33.61D0, -51.79D0, -66.06D0, -91.73D0 ,
     $                       -129.86D0,
     $                        0.0000D0 , 0.000 D0 ,0.000D0 , 0.000D0 ,
     $          -39.82D0 , -56.23D0 , -73.39D0 , -98.99D0 ,.0D0/
      Data UPP3             /   0.0D0, 0.0D0, 0.0D0, 0.0D0,
     $     -25.11D0 , -39.18D0 , -56.40D0 , -78.80D0 , -105.93D0 ,
     $                        0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,
     $         -29.15D0 , -42.31D0 , -57.25D0 , -76.43D0 ,.0D0/
C *** EISOL3 AND EHEAT3 ARE THE GS ELECTRONIC ENERGY OF THE NEUTRAL ATOM
C     (IN E.V.) AND THE HEAT OF FORMATION IF THE FREE ATOM (IN KCAL/MOL)
      Data EISOL3             /-12.505D0 , 0.0D0 , 0.0D0 ,0.0D0 ,
     $        -61.70D0 ,-119.47D0 , -187.51D0 , -307.07D0 , -475.00D0 ,
     $                         0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     $          -90.98D0 , -150.81D0 , -229.15D0 , -345.93D0 , 0.0D0/
      Data EHEAT3             / 52.102D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     $     135.7 D0 , 170.89D0 ,  113.0 D0 ,  59.559D0 ,  18.86D0 ,
     $                         0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     $     106.0D0 ,   79.8D0 ,  65.65D0 ,  28.95D0 , 0.0D0 /
C *** BETA3 AND ALP3 ARE THE BOND PARAMETERS USED IN THE
C     RESONANCE INTEGRAL AND THE CORE CORE REPULSION INTEGRAL RESPECTIVE
C     THAT IS ACCORDING TO THE FOLLOWING CONVENTION
C
C     HERE IS THE
C     BOND TYPE DESIGNATION
C
C
C         H   B   C   N   O   F  SI   P   S  CL
C       -----------------------------------------
C      H  1  11  16  22  29  37  92 106 121 137
C      B     15  20  26  33  41
C      C         21  27  34  42  97 111 126 142
C      N             28  35  43         127 143
C      O                 36  44     113 128 144
C      F                     45     114 129
C     SI                        105
C      P                            120     151
C      S                                136 152
C     CL                                    153
      Data BETA3(  1),ALP3(  1) /  0.244770D0 ,  1.489450D0 /
      Data BETA3( 11),ALP3( 11) /  0.185347D0 ,  2.090352D0 /
      Data BETA3( 15),ALP3( 15) /  0.151324D0 ,  2.280544D0 /
      Data BETA3( 16),ALP3( 16) /  0.315011D0 ,  1.475836D0 /
      Data BETA3( 20),ALP3( 20) /  0.250031D0 ,  2.138291D0 /
      Data BETA3( 21),ALP3( 21) /  0.419907D0 ,  1.371208D0 /
      Data BETA3( 22),ALP3( 22) /  0.360776D0 ,  0.589380D0 /
      Data BETA3( 26),ALP3( 26) /  0.310959D0 ,  1.909763D0 /
      Data BETA3( 27),ALP3( 27) /  0.410886D0 ,  1.635259D0 /
      Data BETA3( 28),ALP3( 28) /  0.377342D0 ,  2.029618D0 /
      Data BETA3( 29),ALP3( 29) /  0.417759D0 ,  0.478901D0 /
      Data BETA3( 33),ALP3( 33) /  0.349745D0 ,  2.484827D0 /
      Data BETA3( 34),ALP3( 34) /  0.464514D0 ,  1.820975D0 /
      Data BETA3( 35),ALP3( 35) /  0.458110D0 ,  1.873859D0 /
      Data BETA3( 36),ALP3( 36) /  0.659407D0 ,  1.537190D0 /
      Data BETA3( 37),ALP3( 37) /  0.195242D0 ,  3.771362D0 /
      Data BETA3( 41),ALP3( 41) /  0.219591D0 ,  2.862183D0 /
      Data BETA3( 42),ALP3( 42) /  0.247494D0 ,  2.725913D0 /
      Data BETA3( 43),ALP3( 43) /  0.205347D0 ,  2.861667D0 /
      Data BETA3( 44),ALP3( 44) /  0.334044D0 ,  2.266949D0 /
      Data BETA3( 45),ALP3( 45) /  0.197464D0 ,  3.864997D0 /
      Data BETA3( 92),ALP3( 92) /  0.289647D0 ,  0.940789D0 /
      Data BETA3( 97),ALP3( 97) /  0.411377D0 ,  1.101382D0 /
      Data BETA3(105),ALP3(105) /  0.291703D0 ,  0.918432D0 /
      Data BETA3(106),ALP3(106) /  0.320118D0 ,  0.923170D0 /
      Data BETA3(111),ALP3(111) /  0.457816D0 ,  1.029693D0 /
      Data BETA3(113),ALP3(113) /  0.470000D0 ,  1.662500D0 /
      Data BETA3(114),ALP3(114) /  0.300000D0 ,  1.750000D0 /
      Data BETA3(120),ALP3(120) /  0.311790D0 ,  1.186652D0 /
      Data BETA3(121),ALP3(121) /  0.220654D0 ,  1.700698D0 /
      Data BETA3(126),ALP3(126) /  0.284620D0 ,  1.761370D0 /
      Data BETA3(127),ALP3(127) /  0.313170D0 ,  1.878176D0 /
      Data BETA3(128),ALP3(128) /  0.422890D0 ,  2.077240D0 /
      Data BETA3(129),ALP3(129) /  0.000000D0 ,  0.000000D0 /
      Data BETA3(136),ALP3(136) /  0.202489D0 ,  1.751617D0 /
      Data BETA3(137),ALP3(137) /  0.231653D0 ,  2.089404D0 /
      Data BETA3(142),ALP3(142) /  0.315480D0 ,  1.676222D0 /
      Data BETA3(143),ALP3(143) /  0.302298D0 ,  1.817064D0 /
      Data BETA3(144),ALP3(144) /  0.000000D0 ,  0.000000D0 /
      Data BETA3(151),ALP3(151) /  0.277322D0 ,  1.543720D0 /
      Data BETA3(152),ALP3(152) /  0.221764D0 ,  1.950318D0 /
      Data BETA3(153),ALP3(153) /  0.258969D0 ,  1.792125D0 /
C *** HERE COMES THE OPTIMIZED SLATER_S EXPONENTS FOR THE EVALUATION
C     OF THE OVERLAP INTEGRALS AND MOLECULAR DIPOLE MOMENTS.
      Data ZS3( 1),ZP3( 1)    /  1.3D0      ,  0.0D0      /
      Data ZS3( 5),ZP3( 5)    /  1.211156D0 ,  0.972826D0 /
      Data ZS3( 6),ZP3( 6)    /  1.739391D0 ,  1.709645D0 /
      Data ZS3( 7),ZP3( 7)    /  2.704546D0 ,  1.870839D0 /
      Data ZS3( 8),ZP3( 8)    /  3.640575D0 ,  2.168448D0 /
      Data ZS3( 9),ZP3( 9)    /  3.111270D0 ,  1.41986D0  /
      Data ZS3(14),ZP3(14)    /  1.629173D0 ,  1.381721D0 /
      Data ZS3(15),ZP3(15)    /  1.926108D0 ,  1.590665D0 /
      Data ZS3(16),ZP3(16)    /  1.719480D0 ,  1.403205D0 /
      Data ZS3(17),ZP3(17)    /  3.430887D0 ,  1.627017D0 /
C************************************************************
C                                                           *
C               Data FOR THE SPARKLES                       *
C                                                           *
C************************************************************
C                               Data FOR THE " ++ " SPARKLE
      Data EHEAT (103) / 0.0D0/
      Data VS    (103) /10.0D0/
      Data ALP   (103) / 1.5D0/
      Data EISOL (103) / 0.0D0/
      Data AM    (103) / 0.5D0/
      Data ALPM  (103) / 1.5D0/
      Data EISOLM(103) / 0.0D0/
      Data AMM   (103) / 0.5D0/
      Data ALPPM3(103) / 1.5D0/
      Data EISOLP(103) / 0.0D0/
      Data AMPM3 (103) / 0.5D0/
      Data DD    (103) / 0.0d0/
      Data DDM   (103) / 0.0d0/
      Data DDPM3 (103) / 0.0d0/
C                               Data FOR THE " + " SPARKLE
      Data EHEAT (104) / 0.0D0/
      Data VS    (104) /10.0D0/
      Data ALP   (104) / 1.5D0/
      Data EISOL (104) / 0.0D0/
      Data AM    (104) / 0.5D0/
      Data ALPM  (104) / 1.5D0/
      Data EISOLM(104) / 0.0D0/
      Data AMM   (104) / 0.5D0/
      Data ALPPM3(104) / 1.5D0/
      Data EISOLP(104) / 0.0D0/
      Data AMPM3 (104) / 0.5D0/
      Data DD    (104) / 0.0d0/
      Data DDM   (104) / 0.0d0/
      Data DDPM3 (104) / 0.0d0/
C                               Data FOR THE " -- " SPARKLE
      Data EHEAT (105) / 0.0D0/
      Data VS    (105) /10.0D0/
      Data ALP   (105) / 1.5D0/
      Data EISOL (105) / 0.0D0/
      Data AM    (105) / 0.5D0/
      Data ALPM  (105) / 1.5D0/
      Data EISOLM(105) / 0.0D0/
      Data AMM   (105) / 0.5D0/
      Data ALPPM3(105) / 1.5D0/
      Data EISOLP(105) / 0.0D0/
      Data AMPM3 (105) / 0.5D0/
      Data DD    (105) / 0.0d0/
      Data DDM   (105) / 0.0d0/
      Data DDPM3 (105) / 0.0d0/
C                               Data FOR THE " - " SPARKLE
      Data EHEAT (106) / 0.0D0/
      Data VS    (106) /10.0D0/
      Data ALP   (106) / 1.5D0/
      Data EISOL (106) / 0.0D0/
      Data AM    (106) / 0.5D0/
      Data ALPM  (106) / 1.5D0/
      Data EISOLM(106) / 0.0D0/
      Data AMM   (106) / 0.5D0/
      Data ALPPM3(106) / 1.5D0/
      Data EISOLP(106) / 0.0D0/
      Data AMPM3 (106) / 0.5D0/
      Data DD    (106) / 0.0d0/
      Data DDM   (106) / 0.0d0/
      Data DDPM3 (106) / 0.0d0/
C**********************************************************************
C
C    START OF MNDO PARAMETERS
C
C**********************************************************************
C                    Data FOR ELEMENT  1        HYDROGEN
      Data RefMN( 1)( 1:40)/' H: (MNDO):  M.J.S. DEWAR. THIEL, J. AM.'/
      Data RefMN( 1)(41:80)/' CHEM. SOC., 99, 4907, (1977)'/
      Data USSM   ( 1)/     -11.9062760D0/
      Data BETASM ( 1)/      -6.9890640D0/
      Data ZSM    ( 1)/       1.3319670D0/
      Data ALPM   ( 1)/       2.5441341D0/
      Data EISOLM ( 1)/     -11.9062760D0/
      Data AMM    ( 1)/       0.4721793D0/
      Data ADM    ( 1)/       0.4721793D0/
      Data AQM    ( 1)/       0.4721793D0/
C
      Data RefMN( 3)( 1:40)/'Li: (MNDO):  TAKEN FROM MNDOC BY W.THIEL'/
      Data RefMN( 3)(41:80)/',QCPE NO.438, V. 2, P.63, (1982).'/
      Data USSM   ( 3)/      -5.1280000D0/
      Data UPPM   ( 3)/      -2.7212000D0/
      Data BETASM ( 3)/      -1.3500400D0/
      Data BETAPM ( 3)/      -1.3500400D0/
      Data ZSM    ( 3)/       0.7023800D0/
      Data ZPM    ( 3)/       0.7023800D0/
      Data ALPM   ( 3)/       1.2501400D0/
      Data EISOLM ( 3)/      -5.1280000D0/
      Data GSSM   ( 3)/       7.3000000D0/
      Data GSPM   ( 3)/       5.4200000D0/
      Data GPPM   ( 3)/       5.0000000D0/
      Data GP2M   ( 3)/       4.5200000D0/
      Data HSPM   ( 3)/       0.8300000D0/
      Data DDM    ( 3)/       2.0549783D0/
      Data QQM    ( 3)/       1.7437069D0/
      Data AMM    ( 3)/       0.2682837D0/
      Data ADM    ( 3)/       0.2269793D0/
      Data AQM    ( 3)/       0.2614581D0/
C                    Data FOR ELEMENT  4        BERYLLIUM
      Data RefMN( 4)( 1:40)/'BE: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J'/
      Data RefMN( 4)(41:80)/'. AM. CHEM. SOC., 100, 777, (1978)'/
      Data USSM   ( 4)/     -16.6023780D0/
      Data UPPM   ( 4)/     -10.7037710D0/
      Data BETASM ( 4)/      -4.0170960D0/
      Data BETAPM ( 4)/      -4.0170960D0/
      Data ZSM    ( 4)/       1.0042100D0/
      Data ZPM    ( 4)/       1.0042100D0/
      Data ALPM   ( 4)/       1.6694340D0/
      Data EISOLM ( 4)/     -24.2047560D0/
      Data DDM    ( 4)/       1.4373245D0/
      Data QQM    ( 4)/       1.2196103D0/
      Data AMM    ( 4)/       0.3307607D0/
      Data ADM    ( 4)/       0.3356142D0/
      Data AQM    ( 4)/       0.3846373D0/
C                    Data FOR ELEMENT  5        BORON
      Data RefMN( 5)( 1:40)/' B: (MNDO):  M.J.S. DEWAR, M.L. MCKEE, J'/
      Data RefMN( 5)(41:80)/'. AM. CHEM. SOC., 99, 5231, (1977)     '/
      Data USSM   ( 5)/     -34.5471300D0/
      Data UPPM   ( 5)/     -23.1216900D0/
      Data BETASM ( 5)/      -8.2520540D0/
      Data BETAPM ( 5)/      -8.2520540D0/
      Data ZSM    ( 5)/       1.5068010D0/
      Data ZPM    ( 5)/       1.5068010D0/
      Data ALPM   ( 5)/       2.1349930D0/
      Data EISOLM ( 5)/     -64.3159500D0/
      Data DDM    ( 5)/       0.9579073D0/
      Data QQM    ( 5)/       0.8128113D0/
      Data AMM    ( 5)/       0.3891951D0/
      Data ADM    ( 5)/       0.4904730D0/
      Data AQM    ( 5)/       0.5556979D0/
C                    Data FOR ELEMENT  6        CARBON
      Data RefMN( 6)( 1:40)/' C: (MNDO):  M.J.S. DEWAR, W. THIEL, J. '/
      Data RefMN( 6)(41:80)/'AM. CHEM. SOC., 99, 4907, (1977)       '/
      Data USSM   ( 6)/     -52.2797450D0/
      Data UPPM   ( 6)/     -39.2055580D0/
      Data BETASM ( 6)/     -18.9850440D0/
      Data BETAPM ( 6)/      -7.9341220D0/
      Data ZSM    ( 6)/       1.7875370D0/
      Data ZPM    ( 6)/       1.7875370D0/
      Data ALPM   ( 6)/       2.5463800D0/
      Data EISOLM ( 6)/    -120.5006060D0/
      Data DDM    ( 6)/       0.8074662D0/
      Data QQM    ( 6)/       0.6851578D0/
      Data AMM    ( 6)/       0.4494671D0/
      Data ADM    ( 6)/       0.6149474D0/
      Data AQM    ( 6)/       0.6685897D0/
C                    Data FOR ELEMENT  7        NITROGEN
      Data RefMN( 7)( 1:40)/' N: (MNDO):  M.J.S. DEWAR, W. THIEL, J. '/
      Data RefMN( 7)(41:80)/'AM. CHEM. SOC., 99, 4907, (1977)       '/
      Data USSM   ( 7)/     -71.9321220D0/
      Data UPPM   ( 7)/     -57.1723190D0/
      Data BETASM ( 7)/     -20.4957580D0/
      Data BETAPM ( 7)/     -20.4957580D0/
      Data ZSM    ( 7)/       2.2556140D0/
      Data ZPM    ( 7)/       2.2556140D0/
      Data ALPM   ( 7)/       2.8613420D0/
      Data EISOLM ( 7)/    -202.5662010D0/
      Data DDM    ( 7)/       0.6399037D0/
      Data QQM    ( 7)/       0.5429763D0/
      Data AMM    ( 7)/       0.4994487D0/
      Data ADM    ( 7)/       0.7843643D0/
      Data AQM    ( 7)/       0.81264450D0/
C                    Data FOR ELEMENT  8        OXYGEN
      Data RefMN( 8)( 1:40)/' O: (MNDO):  M.J.S. DEWAR, W. THIEL, J. '/
      Data RefMN( 8)(41:80)/'AM. CHEM. SOC., 99, 4907, (1977)       '/
      Data USSM   ( 8)/     -99.6443090D0/
      Data UPPM   ( 8)/     -77.7974720D0/
      Data BETASM ( 8)/     -32.6880820D0/
      Data BETAPM ( 8)/     -32.6880820D0/
      Data ZSM    ( 8)/       2.6999050D0/
      Data ZPM    ( 8)/       2.6999050D0/
      Data ALPM   ( 8)/       3.1606040D0/
      Data EISOLM ( 8)/    -317.8685060D0/
      Data DDM    ( 8)/       0.5346024D0/
      Data QQM    ( 8)/       0.4536252D0/
      Data AMM    ( 8)/       0.5667034D0/
      Data ADM    ( 8)/       0.9592562D0/
      Data AQM    ( 8)/       0.9495934D0/
C                    Data FOR ELEMENT  9        FLUORINE
      Data RefMN( 9)( 1:40)/' F: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J'/
      Data RefMN( 9)(41:80)/'. AM. CHEM. SOC., 100, 777, (1978)     '/
      Data USSM   ( 9)/    -131.0715480D0/
      Data UPPM   ( 9)/    -105.7821370D0/
      Data BETASM ( 9)/     -48.2904660D0/
      Data BETAPM ( 9)/     -36.5085400D0/
      Data ZSM    ( 9)/       2.8484870D0/
      Data ZPM    ( 9)/       2.8484870D0/
      Data ALPM   ( 9)/       3.4196606D0/
      Data EISOLM ( 9)/    -476.6837810D0/
      Data DDM    ( 9)/       0.5067166D0/
      Data QQM    ( 9)/       0.4299633D0/
      Data AMM    ( 9)/       0.6218302D0/
      Data ADM    ( 9)/       1.0850301D0/
      Data AQM    ( 9)/       1.0343643D0/
C                               Data FOR THE SODIUM-LIKE SPARKLE
      Data VS     (11) /10.0D0/
      Data ALP    (11) / 1.668D0/
      Data EISOL  (11) / 0.0D0/
      Data DD     (11) / 0.0d0/
      Data AM     (11) / 0.5D0/
      Data ALPM   (11) / 1.660D0/
      Data EISOLM (11) / 0.0D0/
      Data DDM    (11) / 0.0d0/
      Data AMM    (11) / 0.5D0/
C                    Data FOR ELEMENT 13        ALUMINUM
      Data RefMN(13)( 1:40)/'AL: (MNDO):  L.P. DAVIS, ET.AL.  J. COMP'/
      Data RefMN(13)(41:80)/'. CHEM., 2, 433, (1981) SEE MANUAL.    '/
      Data USSM   (13)/     -23.8070970D0/
      Data UPPM   (13)/     -17.5198780D0/
      Data BETASM (13)/      -2.6702840D0/
      Data BETAPM (13)/      -2.6702840D0/
      Data ZSM    (13)/       1.4441610D0/
      Data ZPM    (13)/       1.4441610D0/
      Data ZDM    (13)/       1.0000000D0/
      Data ALPM   (13)/       1.8688394D0/
      Data EISOLM (13)/     -44.4840720D0/
      Data DDM    (13)/       1.3992387D0/
      Data QQM    (13)/       1.1586797D0/
      Data AMM    (13)/       0.2973172D0/
      Data ADM    (13)/       0.2635574D0/
      Data AQM    (13)/       0.3673560D0/
C                    Data FOR ELEMENT 14
C    PARAMETERS FOR SILICON STARTING WITH -37.0375 ARE GILL'S PARAMETERS
      Data RefMN(14)( 1:40)/'SI: (MNDO): M.J.S.DEWAR, ET. AL. ORGANOM'/
      Data RefMN(14)(41:80)/'ETALLICS  5, 375 (1986)                '/
      Data USSM   (14)/     -37.0375330D0/
      Data UPPM   (14)/     -27.7696780D0/
      Data BETASM (14)/      -9.0868040D0/
      Data BETAPM (14)/      -1.0758270D0/
      Data ZSM    (14)/       1.3159860D0/
      Data ZPM    (14)/       1.7099430D0/
      Data ZDM    (14)/       1.0000000D0/
      Data ALPM   (14)/       2.2053160D0/
      Data EISOLM (14)/     -82.8394220D0/
      Data DDM    (14)/       1.2580349D0/
      Data QQM    (14)/       0.9785824D0/
      Data AMM    (14)/       0.3608967D0/
      Data ADM    (14)/       0.3664244D0/
      Data AQM    (14)/       0.4506740D0/
C                    Data FOR ELEMENT 15        PHOSPHORUS
      Data RefMN(15)( 1:40)/' P: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.'/
      Data RefMN(15)(41:80)/'RZEPA,J. AM. CHEM. SOC., 100 3607 1978'/
      Data USSM   (15)/     -56.1433600D0/
      Data UPPM   (15)/     -42.8510800D0/
      Data BETASM (15)/      -6.7916000D0/
      Data BETAPM (15)/      -6.7916000D0/
      Data ZSM    (15)/       2.1087200D0/
      Data ZPM    (15)/       1.7858100D0/
      Data ZDM    (15)/       1.0000000D0/
      Data ALPM   (15)/       2.4152800D0/
      Data EISOLM (15)/    -152.9599600D0/
      Data DDM    (15)/       1.0129699D0/
      Data QQM    (15)/       0.9370090D0/
      Data AMM    (15)/       0.4248438D0/
      Data ADM    (15)/       0.4882420D0/
      Data AQM    (15)/       0.4979406D0/
C
C                  Data for element 16  MNDO: Sulfer
      Data RefMN(16)( 1:40)/' S: (MNDO): M.J.S.DEWAR, C.H. REYNOLDS, '/
      Data RefMN(16)(41:80)/'J. COMP. CHEM. 7, 140-143 (1986)      '/
      Data USSM   (16)/     -72.2422810D0/
      Data UPPM   (16)/     -56.9732070D0/
      Data BETASM (16)/     -10.7616700D0/
      Data BETAPM (16)/     -10.1084330D0/
      Data ZSM    (16)/       2.3129620D0/
      Data ZPM    (16)/       2.0091460D0/
      Data ZDM    (16)/       1.0000000D0/
      Data ALPM   (16)/       2.4780260D0/
      Data EISOLM (16)/    -226.0123900D0/
      Data DDM    (16)/       0.9189935D0/
      Data QQM    (16)/       0.8328514D0/
      Data AMM    (16)/       0.4733554D0/
      Data ADM    (16)/       0.5544502D0/
      Data AQM    (16)/       0.5585244D0/
C                    Data FOR ELEMENT 17        CHLORINE
      Data RefMN(17)( 1:40)/'CL: (MNDO): M.J.S.DEWAR, H.S.RZEPA, J. C'/
      Data RefMN(17)(41:80)/'OMP. CHEM., 4, 158, (1983)             '/
      Data USSM   (17)/    -100.2271660D0/
      Data UPPM   (17)/     -77.3786670D0/
      Data BETASM (17)/     -14.2623200D0/
      Data BETAPM (17)/     -14.2623200D0/
      Data ZSM    (17)/       3.7846450D0/
      Data ZPM    (17)/       2.0362630D0/
      Data ZDM    (17)/       1.0000000D0/
      Data ALPM   (17)/       2.5422010D0/
      Data EISOLM (17)/    -353.1176670D0/
      Data DDM    (17)/       0.4986870D0/
      Data QQM    (17)/       0.8217603D0/
      Data AMM    (17)/       0.5523705D0/
      Data ADM    (17)/       0.8061220D0/
      Data AQM    (17)/       0.6053435D0/
C                               Data FOR THE POTASSIUM-LIKE SPARKLE
      Data RefMN(19)( 1:40)/' K: (MNDO):  POTASSIUM-LIKE SPARKLE.  US'/
      Data RefMN(19)(41:80)/'E WITHCARE.'/
      Data VS     (19) /10.0D0/
      Data ALP    (19) / 1.405D0/
      Data EISOL  (19) / 0.0D0/
      Data DD     (19) / 0.0d0/
      Data DDM    (19) / 0.0d0/
      Data AM     (19) / 0.5D0/
      Data ALPM   (19) / 1.396D0/
      Data EISOLM (19) / 0.0D0/
      Data AMM    (19) / 0.5D0/
C                    Data FOR ELEMENT 24  CHROMIUM
      Data RefMN(24)( 1:40)/'CR: (MNDO):  M.J.S. DEWAR, E.F. HEALY, J'/
      Data RefMN(24)(41:80)/'.J.P.STEWART (IN PREPARATION)      '/
      Data USSM   (24)/     -17.5170270D0/
      Data UPPM   (24)/     -12.5337290D0/
      Data UDDM   (24)/     -44.1249280D0/
      Data BETASM (24)/      -0.1000000D0/
      Data BETAPM (24)/      -0.1000000D0/
      Data BETADM (24)/      -8.7766360D0/
      Data ZSM    (24)/       1.5000000D0/
      Data ZPM    (24)/       1.5000000D0/
      Data ZDM    (24)/       2.8845490D0/
      Data ALPM   (24)/       3.0683070D0/
      Data EISOLM (24)/    -134.8187920D0/
      Data GSSM   (24)/       6.0000000D0/
      Data GSPM   (24)/       4.1500000D0/
      Data GPPM   (24)/       5.0000000D0/
      Data GP2M   (24)/       3.5000000D0/
      Data HSPM   (24)/       1.0000000D0/
      Data GSD    (24)/       2.8746410D0/
      Data GPD    (24)/       3.0000000D0/
      Data GDD    (24)/       8.8949670D0/
      Data DDM    (24)/       1.7320508D0/
      Data QQM    (24)/       1.4142136D0/
      Data AMM    (24)/       0.2205072D0/
      Data ADM    (24)/       0.2711332D0/
      Data AQM    (24)/       0.4464656D0/
C                    Data FOR ELEMENT 30       ZINC
      Data RefMN(30)( 1:40)/'Zn: (MNDO):  M.J.S. DEWAR, K.M. MERZ, OR'/
      Data RefMN(30)(41:80)/'GANOMETALLICS, 5, 1494-1496 (1986)     '/
      Data USSM   (30)/     -20.8397160D0/
      Data UPPM   (30)/     -19.6252240D0/
      Data BETASM (30)/      -1.0000000D0/
      Data BETAPM (30)/      -2.0000000D0/
      Data ZSM    (30)/       2.0473590D0/
      Data ZPM    (30)/       1.4609460D0/
      Data ZDM    (30)/       1.0000000D0/
      Data ALPM   (30)/       1.5064570D0/
      Data EISOLM (30)/     -29.8794320D0/
      Data GSSM   (30)/      11.8000000D0/
      Data GSPM   (30)/      11.1820180D0/
      Data GPPM   (30)/      13.3000000D0/
      Data GP2M   (30)/      12.9305200D0/
      Data HSPM   (30)/       0.4846060D0/
      Data DDM    (30)/       1.3037826D0/
      Data QQM    (30)/       1.4520183D0/
      Data AMM    (30)/       0.4336641D0/
      Data ADM    (30)/       0.2375912D0/
      Data AQM    (30)/       0.2738858D0/
C                    Data FOR ELEMENT 32        GERMANIUM
      Data RefMN(32)( 1:40)/'Ge: (MNDO): M.J.S.DEWAR, G.L.GRADY, E.F.'/
      Data RefMN(32)(41:80)/'HEALY,ORGANOMETALLICS 6 186-189, (1987)'/
      Data USSM   (32)/     -33.9493670D0/
      Data UPPM   (32)/     -27.4251050D0/
      Data BETASM (32)/      -4.5164790D0/
      Data BETAPM (32)/      -1.7555170D0/
      Data ZSM    (32)/       1.2931800D0/
      Data ZPM    (32)/       2.0205640D0/
      Data ALPM   (32)/       1.9784980D0/
      Data EISOLM (32)/     -76.2489440D0/
      Data GSSM   (32)/       9.8000000D0/
      Data GSPM   (32)/       8.3000000D0/
      Data GPPM   (32)/       7.3000000D0/
      Data GP2M   (32)/       6.5000000D0/
      Data HSPM   (32)/       1.3000000D0/
      Data DDM    (32)/       1.2556091D0/
      Data QQM    (32)/       1.0498655D0/
      Data AMM    (32)/       0.3601617D0/
      Data ADM    (32)/       0.3643722D0/
      Data AQM    (32)/       0.4347337D0/
C                    Data FOR ELEMENT 35        BROMINE
      Data RefMN(35)( 1:40)/'BR: (MNDO): M.J.S.DEWAR, E.F. HEALY, J. '/
      Data RefMN(35)(41:80)/'COMP. CHEM., 4, 542, (1983)            '/
      Data USSM   (35)/     -99.9864405D0/
      Data UPPM   (35)/     -75.6713075D0/
      Data BETASM (35)/      -8.9171070D0/
      Data BETAPM (35)/      -9.9437400D0/
      Data ZSM    (35)/       3.8543019D0/
      Data ZPM    (35)/       2.1992091D0/
      Data ZDM    (35)/       1.0000000D0/
      Data ALPM   (35)/       2.4457051D0/
      Data EISOLM (35)/    -346.6812500D0/
      Data DDM    (35)/       0.6051074D0/
      Data QQM    (35)/       0.9645873D0/
      Data AMM    (35)/       0.5526068D0/
      Data ADM    (35)/       0.7258330D0/
      Data AQM    (35)/       0.5574589D0/
C                    Data FOR ELEMENT 50        TIN
      Data RefMN(50)( 1:40)/'SN: (MNDO): M.J.S.DEWAR,G.L.GRADY,J.J.P.'/
      Data RefMN(50)(41:80)/'STEWART, J.AM.CHEM.SOC.,106 6771 (1984)'/
      Data USSM   (50)/     -40.8518020D0/
      Data UPPM   (50)/     -28.5602490D0/
      Data BETASM (50)/      -3.2351470D0/
      Data BETAPM (50)/      -4.2904160D0/
      Data ZSM    (50)/       2.0803800D0/
      Data ZPM    (50)/       1.9371060D0/
      Data ALPM   (50)/       1.8008140D0/
      Data EISOLM (50)/     -92.3241020D0/
      Data GSSM   (50)/       9.8000000D0/
      Data GSPM   (50)/       8.3000000D0/
      Data GPPM   (50)/       7.3000000D0/
      Data GP2M   (50)/       6.5000000D0/
      Data HSPM   (50)/       1.3000000D0/
      Data DDM    (50)/       1.5697766D0/
      Data QQM    (50)/       1.3262292D0/
      Data AMM    (50)/       0.3601617D0/
      Data ADM    (50)/       0.3219998D0/
      Data AQM    (50)/       0.3713827D0/
C                    Data FOR ELEMENT 53        IODINE
      Data RefMN(53)( 1:40)/' I: (MNDO): M.J.S.DEWAR, E.F. HEALY, J.J'/
      Data RefMN(53)(41:80)/'.P. STEWART, J.COMP.CHEM., 5,358,(1984)'/
      Data USSM   (53)/    -100.0030538D0/
      Data UPPM   (53)/     -74.6114692D0/
      Data BETASM (53)/      -7.4144510D0/
      Data BETAPM (53)/      -6.1967810D0/
      Data ZSM    (53)/       2.2729610D0/
      Data ZPM    (53)/       2.1694980D0/
      Data ZDM    (53)/       1.0000000D0/
      Data ALPM   (53)/       2.2073200D0/
      Data EISOLM (53)/    -340.5983600D0/
      Data DDM    (53)/       1.4253233D0/
      Data QQM    (53)/       1.1841707D0/
      Data AMM    (53)/       0.5527541D0/
      Data ADM    (53)/       0.4593451D0/
      Data AQM    (53)/       0.4585376D0/
C                    Data FOR ELEMENT 80        MERCURY
      Data RefMN(80)( 1:40)/'HG: (MNDO): M.J.S.DEWAR,  ET. AL. ORGANO'/
      Data RefMN(80)(41:80)/'METALLICS 4, 1964, (1985) SEE MANUAL   '/
      Data USSM   (80)/     -19.8095740D0/
      Data UPPM   (80)/     -13.1025300D0/
      Data BETASM (80)/      -0.4045250D0/
      Data BETAPM (80)/      -6.2066830D0/
      Data ZSM    (80)/       2.2181840D0/
      Data ZPM    (80)/       2.0650380D0/
      Data ALPM   (80)/       1.3356410D0/
      Data EISOLM (80)/     -28.8191480D0/
      Data GSSM   (80)/      10.8000000D0/
      Data GSPM   (80)/       9.3000000D0/
      Data GPPM   (80)/      14.3000000D0/
      Data GP2M   (80)/      13.5000000D0/
      Data HSPM   (80)/       1.3000000D0/
      Data DDM    (80)/       1.7378048D0/
      Data QQM    (80)/       1.4608064D0/
      Data AMM    (80)/       0.3969129D0/
      Data ADM    (80)/       0.3047694D0/
      Data AQM    (80)/       0.3483102D0/
C                    Data FOR ELEMENT 82        LEAD
      Data RefMN(82)( 1:40)/'PB: (MNDO): M.J.S.DEWAR, ET.AL ORGANOMET'/
      Data RefMN(82)(41:80)/'ALLICS4 1973-1980 (1985)              '/
      Data USSM   (82)/     -47.3196920D0/
      Data UPPM   (82)/     -28.8475600D0/
      Data BETASM (82)/      -8.0423870D0/
      Data BETAPM (82)/      -3.0000000D0/
      Data ZSM    (82)/       2.4982860D0/
      Data ZPM    (82)/       2.0820710D0/
      Data ALPM   (82)/       1.7283330D0/
      Data EISOLM (82)/    -105.8345040D0/
      Data GSSM   (82)/       9.8000000D0/
      Data GSPM   (82)/       8.3000000D0/
      Data GPPM   (82)/       7.3000000D0/
      Data GP2M   (82)/       6.5000000D0/
      Data HSPM   (82)/       1.3000000D0/
      Data DDM    (82)/       1.5526624D0/
      Data QQM    (82)/       1.4488558D0/
      Data AMM    (82)/       0.3601617D0/
      Data ADM    (82)/       0.3239309D0/
      Data AQM    (82)/       0.3502057D0/
C
C     START OF "OLD" ELEMENTS: THESE ARE OLD PARAMETERS WHICH
C     CAN BE USED, IF DESIRED, BY SPECIFYING "<CHEMICAL SYMBOL>YEAR"
C     AS IN SI1978 OR  S1983.
C
C                    Data FOR ELEMENT 90        SILICON
      Data RefMN(90)( 1:40)/'SI: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.'/
      Data RefMN(90)(41:80)/'RZEPA,J. AM. CHEM. SOC., 100 3607 1978'/
      Data USSM   (90)/     -40.5682920D0/
      Data UPPM   (90)/     -28.0891870D0/
      Data BETASM (90)/      -4.2562180D0/
      Data BETAPM (90)/      -4.2562180D0/
      Data ZSM    (90)/       1.4353060D0/
      Data ZPM    (90)/       1.4353060D0/
      Data ZDM    (90)/       1.0000000D0/
      Data ALPM   (90)/       2.1961078D0/
      Data EISOLM (90)/     -90.5399580D0/
      Data DDM    (90)/       1.4078712D0/
      Data QQM    (90)/       1.1658281D0/
      Data AMM    (90)/       0.3608967D0/
      Data ADM    (90)/       0.3441817D0/
      Data AQM    (90)/       0.3999442D0/
      Data HSPM   (90)/       1.32D00/
      Data GP2M   (90)/       6.54D00/
      Data GPPM   (90)/       7.31D00/
      Data GSPM   (90)/       8.36D00/
      Data GSSM   (90)/       9.82D00/
      Data RefMN(91)( 1:40)/' S: (MNDO): M.J.S.DEWAR, H.S. RZEPA, M.L'/
      Data RefMN(91)(41:80)/'.MCKEE, J.AM.CHEM.SOC.100, 3607 (1978).'/
      Data USSM   (91)/     -75.2391520D0/
      Data UPPM   (91)/     -57.8320130D0/
      Data BETASM (91)/     -11.1422310D0/
      Data BETAPM (91)/     -11.1422310D0/
      Data ZSM    (91)/       2.6135910D0/
      Data ZPM    (91)/       2.0343930D0/
      Data ZDM    (91)/       1.0000000D0/
      Data ALPM   (91)/       2.4916445D0/
      Data EISOLM (91)/    -235.4413560D0/
      Data DDM    (91)/       0.8231596D0/
      Data QQM    (91)/       0.8225156D0/
      Data AMM    (91)/       0.4733554D0/
      Data ADM    (91)/       0.5889395D0/
      Data AQM    (91)/       0.5632724D0/
      Data RefMN(102)( 1:40)/'Cb: (MNDO):  Capped Bond  (Hydrogen-like'/
      Data RefMN(102)(41:80)/', takes on a  zero charge.)            '/
      Data USSM  (102)/     -11.9062760D0/
      Data BETASM(102)/-9999999.0000000D0/
      Data ZSM   (102)/       4.0000000D0/
      Data ZPM   (102)/       0.3000000D0/
      Data ZDM   (102)/       0.3000000D0/
      Data ALPM  (102)/       2.5441341D0/
      Data EISOLM(102)/       4.0000000D0/
      Data GSSM  (102)/      12.8480000D0/
      Data HSPM  (102)/       0.1000000D0/
      Data DDM   (102)/       0.0684105D0/
      Data QQM   (102)/       1.0540926D0/
      Data AMM   (102)/       0.4721793D0/
      Data ADM   (102)/       0.9262742D0/
      Data AQM   (102)/       0.2909059D0/
C**********************************************************************
C
C    START OF AM1 PARAMETERS
C
C**********************************************************************
C                    Data FOR ELEMENT  1
      Data RefAM( 1)( 1:40)/' H: (AM1): M.J.S. DEWAR ET AL, J. AM. CH'/
      Data RefAM( 1)(41:80)/'EM. SOC. 107 3902-3909 (1985)'/
      Data USS   ( 1)/     -11.3964270D0/
      Data BETAS ( 1)/      -6.1737870D0/
      Data ZS    ( 1)/       1.1880780D0/
      Data ALP   ( 1)/       2.8823240D0/
      Data EISOL ( 1)/     -11.3964270D0/
      Data GSS   ( 1)/      12.8480000D0/
      Data AM    ( 1)/       0.4721793D0/
      Data AD    ( 1)/       0.4721793D0/
      Data AQ    ( 1)/       0.4721793D0/
      Data Guess1( 1,1)/     0.1227960D0/
      Data Guess2( 1,1)/     5.0000000D0/
      Data Guess3( 1,1)/     1.2000000D0/
      Data Guess1( 1,2)/     0.0050900D0/
      Data Guess2( 1,2)/     5.0000000D0/
      Data Guess3( 1,2)/     1.8000000D0/
      Data Guess1( 1,3)/    -0.0183360D0/
      Data Guess2( 1,3)/     2.0000000D0/
      Data Guess3( 1,3)/     2.1000000D0/
      Data (Guess1(1,J),J=4,10)/7*0.0d0/
      Data (Guess2(1,J),J=4,10)/7*0.0d0/
      Data (Guess3(1,J),J=4,10)/7*0.0d0/
C                    Data for element 3
      Data RefAM( 3)/'LI: (MNDO):  Lithium from work of Thiel?!'/
      Data GSS   ( 3)/       7.300D00/
      Data GPP   ( 3)/       5.0d0/
      Data GSP   ( 3)/       5.43D00/
      Data GP2   ( 3)/       4.52d0/
      Data HSP   ( 3)/       0.83d0/
      Data USS   ( 3)/      -5.128d0/
      Data UPP   ( 3)/      -2.7217d0/
      Data BETAS ( 3)/      -1.3500d0/
      Data BETAP ( 3)/      -1.3500d0/
      Data ZS    ( 3)/       0.70238d0/
      Data ZP    ( 3)/       0.70238d0/
      Data ALP   ( 3)/       1.25014d0/
      Data EISOL ( 3)/      -5.128d0/
      Data DD    ( 3)/       2.05497823d0/
      Data QQ    ( 3)/       1.7437069d0/
      Data AM    ( 3)/       0.2682837d0/
      Data AD    ( 3)/       0.22697935d0/
      Data AQ    ( 3)/       0.26145812d0/
      Data (Guess1(3,J),J=1,10)/10*0.0d0/
      Data (Guess2(3,J),J=1,10)/10*0.0d0/
      Data (Guess3(3,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT  4       AM1:   BERYLLIUM  *
      Data RefAM( 4)( 1:40)/'BE: (MNDO):  M.J.S. DEWAR, H.S. RZEPA, J'/
      Data RefAM( 4)(41:80)/'. AM. CHEM. SOC., 100, 777, (1978)    '/
      Data USS   ( 4)/     -16.6023780D0/
      Data UPP   ( 4)/     -10.7037710D0/
      Data BETAS ( 4)/      -4.0170960D0/
      Data BETAP ( 4)/      -4.0170960D0/
      Data ZS    ( 4)/       1.0042100D0/
      Data ZP    ( 4)/       1.0042100D0/
      Data ALP   ( 4)/       1.6694340D0/
      Data EISOL ( 4)/     -24.2047560D0/
      Data GSS   ( 4)/       9.0000000D0/
      Data GSP   ( 4)/       7.4300000D0/
      Data GPP   ( 4)/       6.9700000D0/
      Data GP2   ( 4)/       6.2200000D0/
      Data HSP   ( 4)/       1.2800000D0/
      Data DD    ( 4)/       1.4373245D0/
      Data QQ    ( 4)/       1.2196103D0/
      Data AM    ( 4)/       0.3307607D0/
      Data AD    ( 4)/       0.3356142D0/
      Data AQ    ( 4)/       0.3846373D0/
      Data (Guess1(4,J),J=1,10)/10*0.0d0/
      Data (Guess2(4,J),J=1,10)/10*0.0d0/
      Data (Guess3(4,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT  5       AM1:   BORON  *
      Data RefAM( 5)( 1:40)/' B: (AM1):  M.J.S. DEWAR, C. JIE, E. G. '/
      Data RefAM( 5)(41:80)/'ZOEBISCH ORGANOMETALLICS 7, 513 (1988)'/
C                    Data FOR ELEMENT  5
      Data USS   ( 5)/     -34.4928700D0/
      Data UPP   ( 5)/     -22.6315250D0/
      Data BETAS ( 5)/      -9.5991140D0/
      Data BetaP ( 5)/      -6.2737570D0/
      Data ZS    ( 5)/       1.6117090D0/
      Data ZP    ( 5)/       1.5553850D0/
      Data ALP   ( 5)/       2.4469090D0/
      Data EISOL ( 5)/     -63.7172650D0/
      Data GSS   ( 5)/      10.5900000D0/
      Data GSP   ( 5)/       9.5600000D0/
      Data GPP   ( 5)/       8.8600000D0/
      Data GP2   ( 5)/       7.8600000D0/
      Data HSP   ( 5)/       1.8100000D0/
      Data DD    ( 5)/       0.9107622D0/
      Data QQ    ( 5)/       0.7874223D0/
      Data AM    ( 5)/       0.3891951D0/
      Data AD    ( 5)/       0.5045152D0/
      Data AQ    ( 5)/       0.5678856D0/
C                    Data FOR ELEMENT  6       AM1:   CARBON
      Data RefAM( 6)( 1:40)/' C: (AM1): M.J.S. DEWAR ET AL, J. AM. CH'/
      Data RefAM( 6)(41:80)/'EM. SOC. 107 3902-3909 (1985)'/
      Data USS   ( 6)/     -52.0286580D0/
      Data UPP   ( 6)/     -39.6142390D0/
      Data BETAS ( 6)/     -15.7157830D0/
      Data BETAP ( 6)/      -7.7192830D0/
      Data ZS    ( 6)/       1.8086650D0/
      Data ZP    ( 6)/       1.6851160D0/
      Data ALP   ( 6)/       2.6482740D0/
      Data EISOL ( 6)/    -120.8157940D0/
      Data GSS   ( 6)/      12.2300000D0/
      Data GSP   ( 6)/      11.4700000D0/
      Data GPP   ( 6)/      11.0800000D0/
      Data GP2   ( 6)/       9.8400000D0/
      Data HSP   ( 6)/       2.4300000D0/
      Data DD    ( 6)/       0.8236736D0/
      Data QQ    ( 6)/       0.7268015D0/
      Data AM    ( 6)/       0.4494671D0/
      Data AD    ( 6)/       0.6082946D0/
      Data AQ    ( 6)/       0.6423492D0/
      Data Guess1( 6,1)/     0.0113550D0/
      Data Guess2( 6,1)/     5.0000000D0/
      Data Guess3( 6,1)/     1.6000000D0/
      Data Guess1( 6,2)/     0.0459240D0/
      Data Guess2( 6,2)/     5.0000000D0/
      Data Guess3( 6,2)/     1.8500000D0/
      Data Guess1( 6,3)/    -0.0200610D0/
      Data Guess2( 6,3)/     5.0000000D0/
      Data Guess3( 6,3)/     2.0500000D0/
      Data Guess1( 6,4)/    -0.0012600D0/
      Data Guess2( 6,4)/     5.0000000D0/
      Data Guess3( 6,4)/     2.6500000D0/
      Data (Guess1(6,J),J=5,10)/6*0.0d0/
      Data (Guess2(6,J),J=5,10)/6*0.0d0/
      Data (Guess3(6,J),J=5,10)/6*0.0d0/
C                    Data FOR ELEMENT  7       AM1:   NITROGEN
      Data RefAM( 7)( 1:40)/' N: (AM1): M.J.S. DEWAR ET AL, J. AM. CH'/
      Data RefAM( 7)(41:80)/'EM. SOC. 107 3902-3909 (1985)'/
      Data USS   ( 7)/     -71.8600000D0/
      Data UPP   ( 7)/     -57.1675810D0/
      Data BETAS ( 7)/     -20.2991100D0/
      Data BETAP ( 7)/     -18.2386660D0/
      Data ZS    ( 7)/       2.3154100D0/
      Data ZP    ( 7)/       2.1579400D0/
      Data ALP   ( 7)/       2.9472860D0/
      Data EISOL ( 7)/    -202.4077430D0/
      Data GSS   ( 7)/      13.5900000D0/
      Data GSP   ( 7)/      12.6600000D0/
      Data GPP   ( 7)/      12.9800000D0/
      Data GP2   ( 7)/      11.5900000D0/
      Data HSP   ( 7)/       3.1400000D0/
      Data DD    ( 7)/       0.6433247D0/
      Data QQ    ( 7)/       0.5675528D0/
      Data AM    ( 7)/       0.4994487D0/
      Data AD    ( 7)/       0.7820840D0/
      Data AQ    ( 7)/       0.7883498D0/
      Data Guess1( 7,1)/     0.0252510D0/
      Data Guess2( 7,1)/     5.0000000D0/
      Data Guess3( 7,1)/     1.5000000D0/
      Data Guess1( 7,2)/     0.0289530D0/
      Data Guess2( 7,2)/     5.0000000D0/
      Data Guess3( 7,2)/     2.1000000D0/
      Data Guess1( 7,3)/    -0.0058060D0/
      Data Guess2( 7,3)/     2.0000000D0/
      Data Guess3( 7,3)/     2.4000000D0/
      Data (Guess1(7,J),J=4,10)/7*0.0d0/
      Data (Guess2(7,J),J=4,10)/7*0.0d0/
      Data (Guess3(7,J),J=4,10)/7*0.0d0/
C                    Data FOR ELEMENT  8       AM1:   OXYGEN
      Data RefAM( 8)( 1:40)/' O: (AM1): M.J.S. DEWAR ET AL, J. AM. CH'/
      Data RefAM( 8)(41:80)/'EM. SOC. 107 3902-3909 (1985)'/
      Data USS   ( 8)/     -97.8300000D0/
      Data UPP   ( 8)/     -78.2623800D0/
      Data BETAS ( 8)/     -29.2727730D0/
      Data BETAP ( 8)/     -29.2727730D0/
      Data ZS    ( 8)/       3.1080320D0/
      Data ZP    ( 8)/       2.5240390D0/
      Data ALP   ( 8)/       4.4553710D0/
      Data EISOL ( 8)/    -316.0995200D0/
      Data GSS   ( 8)/      15.4200000D0/
      Data GSP   ( 8)/      14.4800000D0/
      Data GPP   ( 8)/      14.5200000D0/
      Data GP2   ( 8)/      12.9800000D0/
      Data HSP   ( 8)/       3.9400000D0/
      Data DD    ( 8)/       0.4988896D0/
      Data QQ    ( 8)/       0.4852322D0/
      Data AM    ( 8)/       0.5667034D0/
      Data AD    ( 8)/       0.9961066D0/
      Data AQ    ( 8)/       0.9065223D0/
      Data Guess1( 8,1)/     0.2809620D0/
      Data Guess2( 8,1)/     5.0000000D0/
      Data Guess3( 8,1)/     0.8479180D0/
      Data Guess1( 8,2)/     0.0814300D0/
      Data Guess2( 8,2)/     7.0000000D0/
      Data Guess3( 8,2)/     1.4450710D0/
      Data (Guess1(8,J),J=3,10)/8*0.0d0/
      Data (Guess2(8,J),J=3,10)/8*0.0d0/
      Data (Guess3(8,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  9       AM1:   FLUORINE
      Data RefAM( 9)( 1:40)/' F: (AM1): M.J.S. DEWAR AND E. G. ZOEBIS'/
      Data RefAM( 9)(41:80)/'CH, THEOCHEM, 180, 1 (1988).          '/
      Data USS   ( 9)/    -136.1055790D0/
      Data UPP   ( 9)/    -104.8898850D0/
      Data BETAS ( 9)/     -69.5902770D0/
      Data BETAP ( 9)/     -27.9223600D0/
      Data ZS    ( 9)/       3.7700820D0/
      Data ZP    ( 9)/       2.4946700D0/
      Data ALP   ( 9)/       5.5178000D0/
      Data EISOL ( 9)/    -482.2905830D0/
      Data GSS   ( 9)/      16.9200000D0/
      Data GSP   ( 9)/      17.2500000D0/
      Data GPP   ( 9)/      16.7100000D0/
      Data GP2   ( 9)/      14.9100000D0/
      Data HSP   ( 9)/       4.8300000D0/
      Data DD    ( 9)/       0.4145203D0/
      Data QQ    ( 9)/       0.4909446D0/
      Data AM    ( 9)/       0.6218302D0/
      Data AD    ( 9)/       1.2088792D0/
      Data AQ    ( 9)/       0.9449355D0/
      Data Guess1( 9,1)/     0.2420790D0/
      Data Guess2( 9,1)/     4.8000000D0/
      Data Guess3( 9,1)/     0.9300000D0/
      Data Guess1( 9,2)/     0.0036070D0/
      Data Guess2( 9,2)/     4.6000000D0/
      Data Guess3( 9,2)/     1.6600000D0/
      Data (Guess1(9,J),J=3,10)/8*0.0d0/
      Data (Guess2(9,J),J=3,10)/8*0.0d0/
      Data (Guess3(9,J),J=3,10)/8*0.0d0/
C                               Data FOR THE SODIUM-LIKE SPARKLE
      Data RefAM(11)( 1:40)/'NA: (MNDO):  SODIUM-LIKE SPARKLE.   USE'/
      Data RefAM(11)(41:80)/' WITH CARE.                             '/
      Data (Guess1(11,J),J=1,10)/10*0.0d0/
      Data (Guess2(11,J),J=1,10)/10*0.0d0/
      Data (Guess3(11,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 12       AM1:   Magnesium *
      Data RefAM(12)( 1:40)/'MG: (AM1): M. C. Hutter, ET.AL. JPC. B10'/
      Data RefAM(12)(41:80)/'2 8080 (1998).'/
      Data USS    (12)/    -14.96959313D0/
      Data UPP    (12)/    -11.56229248D0/
      Data BETAS  (12)/     -1.25974355D0/
      Data BETAP  (12)/     -0.77836604D0/
      Data ZS     (12)/      1.22339270D0/
      Data ZP     (12)/      1.02030798D0/
      Data ALP    (12)/      1.67049799D0/
      Data EISOL  (12)/    -22.43786349D0/
      Data GSS    (12)/      7.50132277D0/
      Data GSP    (12)/      6.34591536D0/
      Data GPP    (12)/      4.77534467D0/
      Data GP2    (12)/      4.34017279D0/
      Data HSP    (12)/      0.48930466D0/
      Data DD     (12)/      1.75012115D0/
      Data QQ     (12)/      1.64001467D0/
      Data AM     (12)/      0.27568257D0/
      Data AD     (12)/      0.19957236D0/
      Data AQ     (12)/      0.26440152D0/
      Data GUESS1 (12,1)/    2.55017735D0/
      Data GUESS2 (12,1)/    4.29397225D0/
      Data GUESS3 (12,1)/    0.79989601D0/
      Data GUESS1 (12,2)/   -0.00565806D0/
      Data GUESS2 (12,2)/    2.96053910D0/
      Data GUESS3 (12,2)/    1.47499983D0/
      Data GUESS1 (12,3)/   -0.00610286D0/
      Data GUESS2 (12,3)/    2.61416919D0/
      Data GUESS3 (12,3)/    2.42604040D0/
      Data (Guess1(12,J),J=4,10)/7*0.0d0/
      Data (Guess2(12,J),J=4,10)/7*0.0d0/
      Data (Guess3(12,J),J=4,10)/7*0.0d0/
C                    Data FOR ELEMENT 13       AM1:   ALUMINUM  *
      Data RefAM(13)( 1:40)/'Al: (AM1):  M. J. S. Dewar, A. J. Holder'/
      Data RefAM(13)(41:80)/', Organometallics, 9, 508-511 (1990).'/
      Data USS    (13)/     -24.3535850D0/
      Data UPP    (13)/     -18.3636450D0/
      Data BETAS  (13)/      -3.8668220D0/
      Data BetaP  (13)/      -2.3171460D0/
      Data ZS     (13)/       1.5165930D0/
      Data ZP     (13)/       1.3063470D0/
      Data ZD     (13)/       1.0000000D0/
      Data ALP    (13)/       1.9765860D0/
      Data EISOL  (13)/     -46.4208150D0/
      Data GSS    (13)/       8.0900000D0/
      Data GSP    (13)/       6.6300000D0/
      Data GPP    (13)/       5.9800000D0/
      Data GP2    (13)/       5.4000000D0/
      Data HSP    (13)/       0.7000000D0/
      Data DD     (13)/       1.4040443D0/
      Data QQ     (13)/       1.2809154D0/
      Data AM     (13)/       0.2973172D0/
      Data AD     (13)/       0.2630229D0/
      Data AQ     (13)/       0.3427832D0/
      Data Guess1 (13,1)/     0.0900000D0/
      Data Guess2 (13,1)/    12.3924430D0/
      Data Guess3 (13,1)/     2.0503940D0/
      Data (Guess1(13,J),J=2,10)/9*0.0d0/
      Data (Guess2(13,J),J=2,10)/9*0.0d0/
      Data (Guess3(13,J),J=2,10)/9*0.0d0/
C                    Data FOR ELEMENT 14       AM1:   SILICON  *
      Data RefAM(14)( 1:40)/'Si: (AM1): M.J.S.DEWAR, C. JIE, ORGANOME'/
      Data RefAM(14)(41:80)/'TALLICS, 6, 1486-1490 (1987).'/
      Data USS    (14)/     -33.9536220D0/
      Data UPP    (14)/     -28.9347490D0/
      Data BETAS  (14)/      -3.784852D0/
      Data BETAP  (14)/      -1.968123D0/
      Data ZS     (14)/       1.830697D0/
      Data ZP     (14)/       1.2849530D0/
      Data ZD     (14)/       1.0000000D0/
      Data ALP    (14)/       2.257816D0/
      Data EISOL  (14)/     -79.0017420D0/
      Data GSS    (14)/       9.8200000D0/
      Data GSP    (14)/       8.3600000D0/
      Data GPP    (14)/       7.3100000D0/
      Data GP2    (14)/       6.5400000D0/
      Data HSP    (14)/       1.3200000D0/
      Data DD     (14)/       1.1631107D0/
      Data QQ     (14)/       1.3022422D0/
      Data AM     (14)/       0.3608967D0/
      Data AD     (14)/       0.3829813D0/
      Data AQ     (14)/       0.3712106D0/
      Data Guess1 (14,1)/     0.25D0/
      Data Guess2 (14,1)/     9.000D0/
      Data Guess3 (14,1)/     0.911453D0/
      Data Guess1 (14,2)/     0.061513D0/
      Data Guess2 (14,2)/     5.00D0/
      Data Guess3 (14,2)/     1.995569D0/
      Data Guess1 (14,3)/     0.0207890D0/
      Data Guess2 (14,3)/     5.00D0/
      Data Guess3 (14,3)/     2.990610D0/
      Data (Guess1(14,J),J=4,10)/7*0.0d0/
      Data (Guess2(14,J),J=4,10)/7*0.0d0/
      Data (Guess3(14,J),J=4,10)/7*0.0d0/
C                    Data FOR ELEMENT 15       AM1:   PHOSPHORUS  *
      Data RefAM(15)( 1:40)/' P: (AM1): M.J.S.DEWAR, JIE, C, THEOCHEM'/
      Data RefAM(15)(41:80)/', 187,1 (1989)                        '/
      Data USS    (15)/     -42.0298630D0/
      Data UPP    (15)/     -34.0307090D0/
      Data BETAS  (15)/      -6.3537640D0/
      Data BetaP  (15)/      -6.5907090D0/
      Data ZS     (15)/       1.9812800D0/
      Data ZP     (15)/       1.8751500D0/
      Data ZD     (15)/       1.0000000D0/
      Data ALP    (15)/       2.4553220D0/
      Data EISOL  (15)/    -124.4368355D0/
      Data GSS    (15)/      11.5600050D0/
      Data GSP    (15)/       5.2374490D0/
      Data GPP    (15)/       7.8775890D0/
      Data GP2    (15)/       7.3076480D0/
      Data HSP    (15)/       0.7792380D0/
      Data DD     (15)/       1.0452022D0/
      Data QQ     (15)/       0.8923660D0/
      Data AM     (15)/       0.4248440D0/
      Data AD     (15)/       0.3275319D0/
      Data AQ     (15)/       0.4386854D0/
      Data Guess1 (15,1)/    -0.0318270D0/
      Data Guess2 (15,1)/     6.0000000D0/
      Data Guess3 (15,1)/     1.4743230D0/
      Data Guess1 (15,2)/     0.0184700D0/
      Data Guess2 (15,2)/     7.0000000D0/
      Data Guess3 (15,2)/     1.7793540D0/
      Data Guess1 (15,3)/     0.0332900D0/
      Data Guess2 (15,3)/     9.0000000D0/
      Data Guess3 (15,3)/     3.0065760D0/
      Data (Guess1(15,J),J=4,10)/7*0.0d0/
      Data (Guess2(15,J),J=4,10)/7*0.0d0/
      Data (Guess3(15,J),J=4,10)/7*0.0d0/
C                    Data FOR ELEMENT 16       AM1:   SULFUR  *
      Data RefAM(16)( 1:40)/' S: (AM1): M.J.S.DEWAR and Y-C YUAN, Ino'/
      Data RefAM(16)(41:80)/'rg. Chem. 1990?'/
      Data USS    (16)/     -56.6940560D0/
      Data UPP    (16)/     -48.7170490D0/
      Data BETAS  (16)/      -3.9205660D0/
      Data BetaP  (16)/      -7.9052780D0/
      Data ZS     (16)/       2.3665150D0/
      Data ZP     (16)/       1.6672630D0/
      Data ZD     (16)/       1.0000000D0/
      Data ALP    (16)/       2.4616480D0/
      Data EISOL  (16)/    -191.7321930D0/
      Data GSS    (16)/      11.7863290D0/
      Data GSP    (16)/       8.6631270D0/
      Data GPP    (16)/      10.0393080D0/
      Data GP2    (16)/       7.7816880D0/
      Data HSP    (16)/       2.5321370D0/
      Data DD     (16)/       0.9004265D0/
      Data QQ     (16)/       1.0036329D0/
      Data AM     (16)/       0.4331617D0/
      Data AD     (16)/       0.5907115D0/
      Data AQ     (16)/       0.6454943D0/
      Data Guess1 (16,1)/    -0.5091950D0/
      Data Guess2 (16,1)/     4.5936910D0/
      Data Guess3 (16,1)/     0.7706650D0/
      Data Guess1 (16,2)/    -0.0118630D0/
      Data Guess2 (16,2)/     5.8657310D0/
      Data Guess3 (16,2)/     1.5033130D0/
      Data Guess1 (16,3)/     0.0123340D0/
      Data Guess2 (16,3)/    13.5573360D0/
      Data Guess3 (16,3)/     2.0091730D0/
      Data (Guess1(16,J),J=4,10)/7*0.0d0/
      Data (Guess2(16,J),J=4,10)/7*0.0d0/
      Data (Guess3(16,J),J=4,10)/7*0.0d0/
C                    Data FOR ELEMENT 17       AM1:   CHLORINE
      Data RefAM(17)( 1:40)/'Cl: (AM1): M.J.S. DEWAR AND E. G. ZOEBIS'/
      Data RefAM(17)(41:80)/'CH, THEOCHEM, 180, 1 (1988).           '/
      Data USS    (17)/    -111.6139480D0/
      Data UPP    (17)/     -76.6401070D0/
      Data BETAS  (17)/     -24.5946700D0/
      Data BETAP  (17)/     -14.6372160D0/
      Data ZS     (17)/       3.6313760D0/
      Data ZP     (17)/       2.0767990D0/
      Data ZD     (17)/       1.0000000D0/
      Data ALP    (17)/       2.9193680D0/
      Data EISOL  (17)/    -372.1984310D0/
      Data GSS    (17)/      15.0300000D0/
      Data GSP    (17)/      13.1600000D0/
      Data GPP    (17)/      11.3000000D0/
      Data GP2    (17)/       9.9700000D0/
      Data HSP    (17)/       2.4200000D0/
      Data DD     (17)/       0.5406286D0/
      Data QQ     (17)/       0.8057208D0/
      Data AM     (17)/       0.5523705D0/
      Data AD     (17)/       0.7693200D0/
      Data AQ     (17)/       0.6133369D0/
      Data Guess1 (17,1)/     0.0942430D0/
      Data Guess2 (17,1)/     4.0000000D0/
      Data Guess3 (17,1)/     1.3000000D0/
      Data Guess1 (17,2)/     0.0271680D0/
      Data Guess2 (17,2)/     4.0000000D0/
      Data Guess3 (17,2)/     2.1000000D0/
      Data (Guess1(17,J),J=3,10)/8*0.0d0/
      Data (Guess2(17,J),J=3,10)/8*0.0d0/
      Data (Guess3(17,J),J=3,10)/8*0.0d0/
C                               Data FOR THE POTASSIUM-LIKE SPARKLE
      Data RefAM(19)( 1:40)/'K:  (AM1):  POTASSIUM-LIKE SPARKLE.   US'/
      Data RefAM(19)(41:80)/'E WITH CARE.                          '/
      Data (Guess1(19,J),J=1,10)/10*0.0d0/
      Data (Guess2(19,J),J=1,10)/10*0.0d0/
      Data (Guess3(19,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 24  CHROMIUM         *
      Data RefAM(24)( 1:40)/'CR: (MNDO):  M.J.S. DEWAR, E.F. HEALY, J'/
      Data RefAM(24)(41:80)/'.J.P.STEWART (IN PREPARATION)      '/
      Data USS    (24)/     -17.5170270D0/
      Data UPP    (24)/     -12.5337290D0/
      Data UDD    (24)/     -44.1249280D0/
      Data BETAS  (24)/      -0.1000000D0/
      Data BETAP  (24)/      -0.1000000D0/
      Data BETAD  (24)/      -8.7766360D0/
      Data ZS     (24)/       1.5000000D0/
      Data ZP     (24)/       1.5000000D0/
      Data ZD     (24)/       2.8845490D0/
      Data ALP    (24)/       3.0683070D0/
      Data EISOL  (24)/    -134.8187920D0/
      Data GSS    (24)/       6.0000000D0/
      Data GSP    (24)/       4.1500000D0/
      Data GPP    (24)/       5.0000000D0/
      Data GP2    (24)/       3.5000000D0/
      Data HSP    (24)/       1.0000000D0/
      Data DD     (24)/       1.7320508D0/
      Data QQ     (24)/       1.4142136D0/
      Data AM     (24)/       0.2205072D0/
      Data AD     (24)/       0.2711332D0/
      Data AQ     (24)/       0.4464656D0/
C     Data ALPTM  (24)/       2.9499070D0/
C     Data EMUDTM (24)/       2.2085050D0/
      Data ALPTM  (24)/       0.0d0/
      Data EMUDTM (24)/       0.0d0/
      Data (Guess1(24,J),J=1,10)/10*0.0d0/
      Data (Guess2(24,J),J=1,10)/10*0.0d0/
      Data (Guess3(24,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 30
      Data USS    (30)/     -21.0400080D0/
      Data UPP    (30)/     -17.6555740D0/
      Data BETAS  (30)/      -1.9974290D0/
      Data BETAP  (30)/      -4.7581190D0/
      Data ZS     (30)/       1.9542990D0/
      Data ZP     (30)/       1.3723650D0/
      Data ZD     (30)/       1.0000000D0/
      Data ALP    (30)/       1.4845630D0/
      Data EISOL  (30)/     -30.2800160D0/
      Data GSS    (30)/      11.8000000D0/
      Data GSP    (30)/      11.1820180D0/
      Data GPP    (30)/      13.3000000D0/
      Data GP2    (30)/      12.9305200D0/
      Data HSP    (30)/       0.4846060D0/
      Data DD     (30)/       1.3581113D0/
      Data QQ     (30)/       1.5457406D0/
      Data AM     (30)/       0.4336641D0/
      Data AD     (30)/       0.2317423D0/
      Data AQ     (30)/       0.2621165D0/
      Data (Guess1(30,J),J=1,10)/10*0.0d0/
      Data (Guess2(30,J),J=1,10)/10*0.0d0/
      Data (Guess3(30,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 32        GERMANIUM *
      Data RefAM(32)( 1:40)/'Ge: (AM1): M.J.S.Dewar and C.Jie, Organo'/
      Data RefAM(32)(41:80)/'metallics, 8, 1544, (1989)             '/
      Data USS    (32)/     -34.1838890D0/
      Data UPP    (32)/     -28.6408110D0/
      Data BETAS  (32)/      -4.3566070D0/
      Data BETAP  (32)/      -0.9910910D0/
      Data ZS     (32)/       1.2196310D0/
      Data ZP     (32)/       1.9827940D0/
      Data ZD     (32)/       1.0000000D0/
      Data ALP    (32)/       2.1364050D0/
      Data EISOL  (32)/     -78.7084810D0/
      Data GSS    (32)/      10.1686050D0/
      Data GSP    (32)/       8.1444730D0/
      Data GPP    (32)/       6.6719020D0/
      Data GP2    (32)/       6.2697060D0/
      Data HSP    (32)/       0.9370930D0/
      Data DD     (32)/       1.2472095D0/
      Data QQ     (32)/       1.0698642D0/
      Data AM     (32)/       0.3737084D0/
      Data AD     (32)/       0.3180310D0/
      Data AQ     (32)/       0.3485612D0/
      Data (Guess1(32,J),J=1,10)/10*0.0d0/
      Data (Guess2(32,J),J=1,10)/10*0.0d0/
      Data (Guess3(32,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 35       AM1:   BROMINE
      Data RefAM(35)( 1:40)/'Br: (AM1): M.J.S. DEWAR AND E. G. ZOEBIS'/
      Data RefAM(35)(41:80)/'CH, THEOCHEM, 180, 1 (1988).          '/
      Data USS    (35)/    -104.6560630D0/
      Data UPP    (35)/     -74.9300520D0/
      Data BETAS  (35)/     -19.3998800D0/
      Data BETAP  (35)/      -8.9571950D0/
      Data ZS     (35)/       3.0641330D0/
      Data ZP     (35)/       2.0383330D0/
      Data ZD     (35)/       1.0000000D0/
      Data ALP    (35)/       2.5765460D0/
      Data EISOL  (35)/    -352.3142087D0/
      Data GSS    (35)/      15.0364395D0/
      Data GSP    (35)/      13.0346824D0/
      Data GPP    (35)/      11.2763254D0/
      Data GP2    (35)/       9.8544255D0/
      Data HSP    (35)/       2.4558683D0/
      Data DD     (35)/       0.8458104D0/
      Data QQ     (35)/       1.0407133D0/
      Data AM     (35)/       0.5526071D0/
      Data AD     (35)/       0.6024598D0/
      Data AQ     (35)/       0.5307555D0/
      Data Guess1 (35,1)/     0.0666850D0/
      Data Guess2 (35,1)/     4.0000000D0/
      Data Guess3 (35,1)/     1.5000000D0/
      Data Guess1 (35,2)/     0.0255680D0/
      Data Guess2 (35,2)/     4.0000000D0/
      Data Guess3 (35,2)/     2.3000000D0/
      Data (Guess1(35,J),J=3,10)/8*0.0d0/
      Data (Guess2(35,J),J=3,10)/8*0.0d0/
      Data (Guess3(35,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 50        TIN  *
      Data RefAM(50)/'SN: (AM1): taken from AMPAC 2.1'/
      Data USS    (50)/     -35.4967410D0/
      Data UPP    (50)/     -28.0976360D0/
      Data BETAS  (50)/      -3.2350000D0/
      Data BETAP  (50)/      -2.5778900D0/
      Data ZS     (50)/       2.5993760D0/
      Data ZP     (50)/       1.6959620D0/
      Data ZD     (50)/       1.0000000D0/
      Data ALP    (50)/       1.8369360D0/
      Data EISOL  (50)/     -80.6887540D0/
      Data GSS    (50)/       9.8000000D0/
      Data GSP    (50)/       8.3000000D0/
      Data GPP    (50)/       7.3000000D0/
      Data GP2    (50)/       6.5000000D0/
      Data HSP    (50)/       1.3000000D0/
      Data DD     (50)/       1.1528229D0/
      Data QQ     (50)/       1.5148019D0/
      Data AM     (50)/       0.3601617D0/
      Data AD     (50)/       0.3823775D0/
      Data AQ     (50)/       0.3400755D0/
      Data (Guess1(50,J),J=1,10)/10*0.0d0/
      Data (Guess2(50,J),J=1,10)/10*0.0d0/
      Data (Guess3(50,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 53       AM1:   IODINE
      Data RefAM(53)/' I: (AM1): DEWAR AND ZOEBISCH, UNPUBLISHED WORK'/
      Data USS    (53)/    -103.5896630D0/
      Data UPP    (53)/     -74.4299970D0/
      Data BETAS  (53)/      -8.4433270D0/
      Data BETAP  (53)/      -6.3234050D0/
      Data ZS     (53)/       2.1028580D0/
      Data ZP     (53)/       2.1611530D0/
      Data ZD     (53)/       1.0000000D0/
      Data ALP    (53)/       2.2994240D0/
      Data EISOL  (53)/    -346.8642857D0/
      Data GSS    (53)/      15.0404486D0/
      Data GSP    (53)/      13.0565580D0/
      Data GPP    (53)/      11.1477837D0/
      Data GP2    (53)/       9.9140907D0/
      Data HSP    (53)/       2.4563820D0/
      Data DD     (53)/       1.4878778D0/
      Data QQ     (53)/       1.1887388D0/
      Data AM     (53)/       0.5527544D0/
      Data AD     (53)/       0.4497523D0/
      Data AQ     (53)/       0.4631775D0/
      Data Guess1 (53,1)/     0.0043610D0/
      Data Guess2 (53,1)/     2.3000000D0/
      Data Guess3 (53,1)/     1.8000000D0/
      Data Guess1 (53,2)/     0.0157060D0/
      Data Guess2 (53,2)/     3.0000000D0/
      Data Guess3 (53,2)/     2.2400000D0/
      Data (Guess1(53,J),J=3,10)/8*0.0d0/
      Data (Guess2(53,J),J=3,10)/8*0.0d0/
      Data (Guess3(53,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 54
      Data RefAM(54)/' Xe: (MNDO?) source unknown'/
      Data USS    (54)/    -100.0030538D0/
      Data UPP    (54)/     -70.0000000D0/
      Data BETAS  (54)/      -7.4144510D0/
      Data BETAP  (54)/      -6.1967810D0/
      Data ZS     (54)/       2.2729610D0/
      Data ZP     (54)/       1.1694980D0/
      Data ALP    (54)/       2.2073200D0/
      Data EISOL  (54)/    -191.0061076D0/
      Data GSS    (54)/       9.0000000D0/
      Data GSP    (54)/      11.2000000D0/
      Data GPP    (54)/      11.0000000D0/
      Data GP2    (54)/       9.2300000D0/
      Data HSP    (54)/       2.4000000D0/
      Data DD     (54)/       1.0775919D0/
      Data QQ     (54)/       2.5794150D0/
      Data AM     (54)/       0.3307607D0/
      Data AD     (54)/       0.5236997D0/
      Data AQ     (54)/       0.3397123D0/
      Data (Guess1(54,J),J=1,10)/10*0.0d0/
      Data (Guess2(54,J),J=1,10)/10*0.0d0/
      Data (Guess3(54,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 80        MERCURY *
      Data RefAM(80)( 1:40)/'Hg: (AM1): M.J.S.Dewar and C.Jie, Organo'/
      Data RefAM(80)(41:80)/'metallics 8, 1547, (1989)              '/
      Data USS    (80)/     -19.9415780D0/
      Data UPP    (80)/     -11.1108700D0/
      Data BETAS  (80)/      -0.9086570D0/
      Data BETAP  (80)/      -4.9093840D0/
      Data ZS     (80)/       2.0364130D0/
      Data ZP     (80)/       1.9557660D0/
      Data ZD     (80)/       1.0000000D0/
      Data ALP    (80)/       1.4847340D0/
      Data EISOL  (80)/     -29.0831560D0/
      Data GSS    (80)/      10.8000000D0/
      Data GSP    (80)/       9.3000000D0/
      Data GPP    (80)/      14.3000000D0/
      Data GP2    (80)/      13.5000000D0/
      Data HSP    (80)/       1.3000000D0/
      Data DD     (80)/       1.8750829D0/
      Data QQ     (80)/       1.5424241D0/
      Data AM     (80)/       0.3969129D0/
      Data AD     (80)/       0.2926605D0/
      Data AQ     (80)/       0.3360599D0/
      Data (Guess1(80,J),J=1,10)/10*0.0d0/
      Data (Guess2(80,J),J=1,10)/10*0.0d0/
      Data (Guess3(80,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 82        LEAD  *
      Data RefAM(82)/' PB: (MNDO): M.J.S.DEWAR, ET.AL (SUBMITTED)'/
      Data USS    (82)/     -47.3196920D0/
      Data UPP    (82)/     -28.8475600D0/
      Data BETAS  (82)/      -8.0423870D0/
      Data BETAP  (82)/      -3.0000000D0/
      Data ZS     (82)/       2.4982860D0/
      Data ZP     (82)/       2.0820710D0/
      Data ALP    (82)/       1.7283330D0/
      Data EISOL  (82)/    -105.8345040D0/
      Data GSS    (82)/       9.8000000D0/
      Data GSP    (82)/       8.3000000D0/
      Data GPP    (82)/       7.3000000D0/
      Data GP2    (82)/       6.5000000D0/
      Data HSP    (82)/       1.3000000D0/
      Data DD     (82)/       1.5526624D0/
      Data QQ     (82)/       1.4488558D0/
      Data AM     (82)/       0.3601617D0/
      Data AD     (82)/       0.3239309D0/
      Data AQ     (82)/       0.3502057D0/
      Data (Guess1(82,J),J=1,10)/10*0.0d0/
      Data (Guess2(82,J),J=1,10)/10*0.0d0/
      Data (Guess3(82,J),J=1,10)/10*0.0d0/
C
C     START OF "OLD" ELEMENTS: THESE ARE OLD PARAMETERS WHICH
C     CAN BE USED, IF DESIRED, BY SPECIFYING "<CHEMICAL SYMBOL>YEAR"
C     AS IN SI1978 OR  S1983.
C
C                    Data FOR ELEMENT 90        SILICON
      Data RefAM(90)( 1:40)/'SI: (MNDO): M.J.S.DEWAR, M.L.MCKEE, H.S.'/
      Data RefAM(90)(41:80)/'RZEPA,J. AM. CHEM. SOC., 100 3607 1978'/
      Data USS    (90)/     -40.5682920D0/
      Data UPP    (90)/     -28.0891870D0/
      Data BETAS  (90)/      -4.2562180D0/
      Data BETAP  (90)/      -4.2562180D0/
      Data ZS     (90)/       1.4353060D0/
      Data ZP     (90)/       1.4353060D0/
      Data ZD     (90)/       1.0000000D0/
      Data ALP    (90)/       2.1961078D0/
      Data EISOL  (90)/     -90.5399580D0/
      Data DD     (90)/       1.4078712D0/
      Data QQ     (90)/       1.1658281D0/
      Data AM     (90)/       0.3608967D0/
      Data AD     (90)/       0.3441817D0/
      Data AQ     (90)/       0.3999442D0/
      Data HSP    (90)/       1.32D00/
      Data GP2    (90)/       6.54D00/
      Data GPP    (90)/       7.31D00/
      Data GSP    (90)/       8.36D00/
      Data GSS    (90)/       9.82D00/
      Data (Guess1(90,J),J=1,10)/10*0.0d0/
      Data (Guess2(90,J),J=1,10)/10*0.0d0/
      Data (Guess3(90,J),J=1,10)/10*0.0d0/
      Data RefAM(91)( 1:40)/' S: (MNDO): M.J.S.DEWAR, H.S. RZEPA, M.L'/
      Data RefAM(91)(41:80)/'.MCKEE, J.AM.CHEM.SOC.100, 3607 (1978).'/
      Data USS    (91)/     -75.2391520D0/
      Data UPP    (91)/     -57.8320130D0/
      Data BETAS  (91)/     -11.1422310D0/
      Data BETAP  (91)/     -11.1422310D0/
      Data ZS     (91)/       2.6135910D0/
      Data ZP     (91)/       2.0343930D0/
      Data ZD     (91)/       1.0000000D0/
      Data ALP    (91)/       2.4916445D0/
      Data EISOL  (91)/    -235.4413560D0/
      Data GSS    (91)/      12.8800000D0/
      Data GSP    (91)/      11.2600000D0/
      Data GPP    (91)/       9.9000000D0/
      Data GP2    (91)/       8.8300000D0/
      Data HSP    (91)/       2.2600000D0/
      Data DD     (91)/       0.8231596D0/
      Data QQ     (91)/       0.8225156D0/
      Data AM     (91)/       0.4733554D0/
      Data AD     (91)/       0.5889395D0/
      Data AQ     (91)/       0.5632724D0/
      Data (Guess1(91,J),J=1,10)/10*0.0d0/
      Data (Guess2(91,J),J=1,10)/10*0.0d0/
      Data (Guess3(91,J),J=1,10)/10*0.0d0/
      Data RefAM(102)( 1:40)/'Cb: (AM1):  Capped Bond  (Hydrogen-like,'/
      Data RefAM(102)(41:80)/' takes on zero charge.)'/
C                    Data FOR ELEMENT102
      Data USS    (102)/     -11.9062760D0/
      Data BETAS  (102)/-9999999.0000000D0/
      Data ZS     (102)/       4.0000000D0/
      Data ZP     (102)/       0.3000000D0/
      Data ZD     (102)/       0.3000000D0/
      Data ALP    (102)/       2.5441341D0/
      Data EISOL  (102)/       4.0000000D0/
      Data GSS    (102)/      12.8480000D0/
      Data HSP    (102)/       0.1000000D0/
      Data DD     (102)/       0.0684105D0/
      Data QQ     (102)/       1.0540926D0/
      Data AM     (102)/       0.4721793D0/
      Data AD     (102)/       0.9262742D0/
      Data AQ     (102)/       0.2909059D0/
      Data (Guess1(102,J),J=1,10)/10*0.0d0/
      Data (Guess2(102,J),J=1,10)/10*0.0d0/
      Data (Guess3(102,J),J=1,10)/10*0.0d0/
C
C    NOTE "*" MEANS THAT AM1 PARAMETERS ARE NOT AVAILABLE YET, AND THAT
C             MNDO PARAMETERS ARE TO BE USED MEANWHILE, IF THOUGHT
C             SUITABLE.
C
C     START OF MNDO-PM3 PARAMETER SET
C
C                    Data FOR ELEMENT  1        HYDROGEN
      Data REFPM3( 1)( 1:40)/' H: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3( 1)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 ( 1)/     -13.0733210D0/
      Data BETASP ( 1)/      -5.6265120D0/
      Data ZSPM3  ( 1)/       0.9678070D0/
      Data ALPPM3 ( 1)/       3.3563860D0/
      Data EISOLP ( 1)/     -13.0733210D0/
      Data GSSPM3 ( 1)/      14.7942080D0/
      Data AMPM3  ( 1)/       0.5437048D0/
      Data ADPM3  ( 1)/       0.5437048D0/
      Data AQPM3  ( 1)/       0.5437048D0/
      Data GuesP1 ( 1,1)/     1.1287500D0/
      Data GuesP2 ( 1,1)/     5.0962820D0/
      Data GuesP3 ( 1,1)/     1.5374650D0/
      Data GuesP1 ( 1,2)/    -1.0603290D0/
      Data GuesP2 ( 1,2)/     6.0037880D0/
      Data GuesP3 ( 1,2)/     1.5701890D0/
      Data (GuesP1(1,J),J=3,10)/8*0.0d0/
      Data (GuesP2(1,J),J=3,10)/8*0.0d0/
      Data (GuesP3(1,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  3        LITHIUM
      Data REFPM3( 3)( 1:40)/'Li: (PM3): ANDERS, KOCH, FREUNSCHT,J. CO'/
      Data REFPM3( 3)(41:80)/'MP. CHEM. 14, 1301 (1993).'/
      Data USSPM3 ( 3)/      -5.3000000D0/
      Data UPPPM3 ( 3)/      -3.4000000D0/
      Data BETASP ( 3)/      -0.5500000D0/
      Data BETAPP ( 3)/      -1.5000000D0/
      Data ZSPM3  ( 3)/       0.6500000D0/
      Data ZPPM3  ( 3)/       0.7500000D0/
      Data ALPPM3 ( 3)/       1.2550000D0/
      Data GSSPM3 ( 3)/       4.5000000D0/
      Data GSPPM3 ( 3)/       3.0000000D0/
      Data GPPPM3 ( 3)/       5.2500000D0/
      Data GP2PM3 ( 3)/       4.5000000D0/
      Data HSPPM3 ( 3)/       0.1500000D0/
      Data DDPM3  ( 3)/       2.0357652D0/
      Data QQPM3  ( 3)/       1.6329932D0/
      Data AMPM3  ( 3)/       0.1653804D0/
      Data ADPM3  ( 3)/       0.1156738D0/
      Data AQPM3  ( 3)/       0.3166080D0/
      Data GUESP1 ( 3,1)/    -0.4500000D0/
      Data GUESP2 ( 3,1)/     5.0000000D0/
      Data GUESP3 ( 3,1)/     1.0000000D0/
      Data GUESP1 ( 3,2)/     0.8000000D0/
      Data GUESP2 ( 3,2)/     6.5000000D0/
      Data GUESP3 ( 3,2)/     1.0000000D0/
C                    Old Data FOR ELEMENT  3
C     Data USSPM3 ( 3)/      -4.9414130D0/
C     Data UPPPM3 ( 3)/      -4.8996300D0/
C     Data BETASP ( 3)/      -0.2783390D0/
C     Data BETAPP ( 3)/      -2.8900310D0/
C     Data ZSPM3  ( 3)/       1.5000000D0/
C     Data ZPPM3  ( 3)/       0.7000000D0/
C     Data ALPPM3 ( 3)/       1.2303780D0/
C     Data EISOLP ( 3)/      -4.9414130D0/
C     Data GSSPM3 ( 3)/       4.2382110D0/
C     Data GSPPM3 ( 3)/      16.2028110D0/
C     Data GPPPM3 ( 3)/       5.8425930D0/
C     Data GP2PM3 ( 3)/      13.2881150D0/
C     Data HSPPM3 ( 3)/       0.0000001D0/
C     Data DDPM3  ( 3)/       0.9204422D0/
C     Data QQPM3  ( 3)/       1.7496355D0/
C     Data AMPM3  ( 3)/       0.1557593D0/
C     Data ADPM3  ( 3)/       0.0016309D0/
C     Data AQPM3  ( 3)/       0.2010255D0/
C     Data GuesP1 ( 3,1)/     1.7757820D0/
C     Data GuesP2 ( 3,1)/     5.7070980D0/
C     Data GuesP3 ( 3,1)/     1.1333310D0/
C     Data GuesP1 ( 3,2)/    -1.8605030D0/
C     Data GuesP2 ( 3,2)/     3.8427580D0/
C     Data GuesP3 ( 3,2)/     1.0267740D0/
C     Data (GuesP1(3,J),J=3,10)/8*0.0d0/
C     Data (GuesP2(3,J),J=3,10)/8*0.0d0/
C     Data (GuesP3(3,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  4        BERYLLIUM
      Data REFPM3( 4)( 1:40)/'Be: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3( 4)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 ( 4)/     -17.2647520D0/
      Data UPPPM3 ( 4)/     -11.3042430D0/
      Data BETASP ( 4)/      -3.9620530D0/
      Data BETAPP ( 4)/      -2.7806840D0/
      Data ZSPM3  ( 4)/       0.8774390D0/
      Data ZPPM3  ( 4)/       1.5087550D0/
      Data ALPPM3 ( 4)/       1.5935360D0/
      Data EISOLP ( 4)/     -25.5166530D0/
      Data GSSPM3 ( 4)/       9.0128510D0/
      Data GSPPM3 ( 4)/       6.5761990D0/
      Data GPPPM3 ( 4)/       6.0571820D0/
      Data GP2PM3 ( 4)/       9.0052190D0/
      Data HSPPM3 ( 4)/       0.5446790D0/
      Data DDPM3  ( 4)/       1.0090531D0/
      Data QQPM3  ( 4)/       0.8117586D0/
      Data AMPM3  ( 4)/       0.3312330D0/
      Data ADPM3  ( 4)/       0.2908996D0/
      Data AQPM3  ( 4)/       0.3530008D0/
      Data GuesP1 ( 4,1)/     1.6315720D0/
      Data GuesP2 ( 4,1)/     2.6729620D0/
      Data GuesP3 ( 4,1)/     1.7916860D0/
      Data GuesP1 ( 4,2)/    -2.1109590D0/
      Data GuesP2 ( 4,2)/     1.9685940D0/
      Data GuesP3 ( 4,2)/     1.7558710D0/
      Data (GuesP1( 4,J),J=3,10)/8*0.0d0/
      Data (GuesP2( 4,J),J=3,10)/8*0.0d0/
      Data (GuesP3( 4,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  5
      Data USSPM3 ( 5)/     -28.3539160D0/
      Data UPPPM3 ( 5)/     -23.9738780D0/
      Data BETASP ( 5)/      -5.5215820D0/
      Data BETAPP ( 5)/      -5.9611770D0/
      Data ZSPM3  ( 5)/       0.9885410D0/
      Data ZPPM3  ( 5)/       1.9094410D0/
      Data ALPPM3 ( 5)/       2.4981070D0/
      Data EISOLP ( 5)/     -52.1098330D0/
      Data GSSPM3 ( 5)/       8.8104530D0/
      Data GSPPM3 ( 5)/      10.2992280D0/
      Data GPPPM3 ( 5)/       4.9891680D0/
      Data GP2PM3 ( 5)/       7.2899770D0/
      Data HSPPM3 ( 5)/       0.8370320D0/
      Data DDPM3  ( 5)/       0.7633744D0/
      Data QQPM3  ( 5)/       0.6414154D0/
      Data AMPM3  ( 5)/       0.3237947D0/
      Data ADPM3  ( 5)/       0.4075660D0/
      Data AQPM3  ( 5)/       0.4208854D0/
      Data GuesP1 ( 5,1)/     3.0046640D0/
      Data GuesP2 ( 5,1)/     6.0288690D0/
      Data GuesP3 ( 5,1)/     0.3472230D0/
      Data GuesP1 ( 5,2)/     0.0201950D0/
      Data GuesP2 ( 5,2)/     6.0106580D0/
      Data GuesP3 ( 5,2)/     1.2140670D0/
      Data (GuesP1( 5,J),J=3,10)/8*0.0d0/
      Data (GuesP2( 5,J),J=3,10)/8*0.0d0/
      Data (GuesP3( 5,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  6      CARBON
      Data REFPM3( 6)( 1:40)/' C: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3( 6)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 ( 6)/     -47.2703200D0/
      Data UPPPM3 ( 6)/     -36.2669180D0/
      Data BETASP ( 6)/     -11.9100150D0/
      Data BETAPP ( 6)/      -9.8027550D0/
      Data ZSPM3  ( 6)/       1.5650850D0/
      Data ZPPM3  ( 6)/       1.8423450D0/
      Data ALPPM3 ( 6)/       2.7078070D0/
      Data EISOLP ( 6)/    -111.2299170D0/
      Data GSSPM3 ( 6)/      11.2007080D0/
      Data GSPPM3 ( 6)/      10.2650270D0/
      Data GPPPM3 ( 6)/      10.7962920D0/
      Data GP2PM3 ( 6)/       9.0425660D0/
      Data HSPPM3 ( 6)/       2.2909800D0/
      Data DDPM3  ( 6)/       0.8332396D0/
      Data QQPM3  ( 6)/       0.6647750D0/
      Data AMPM3  ( 6)/       0.4116394D0/
      Data ADPM3  ( 6)/       0.5885862D0/
      Data AQPM3  ( 6)/       0.7647667D0/
      Data GuesP1 ( 6,1)/     0.0501070D0/
      Data GuesP2 ( 6,1)/     6.0031650D0/
      Data GuesP3 ( 6,1)/     1.6422140D0/
      Data GuesP1 ( 6,2)/     0.0507330D0/
      Data GuesP2 ( 6,2)/     6.0029790D0/
      Data GuesP3 ( 6,2)/     0.8924880D0/
      Data (GuesP1( 6,J),J=3,10)/8*0.0d0/
      Data (GuesP2( 6,J),J=3,10)/8*0.0d0/
      Data (GuesP3( 6,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  7      NITROGEN
      Data REFPM3( 7)( 1:40)/' N: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3( 7)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 ( 7)/     -49.3356720D0/
      Data UPPPM3 ( 7)/     -47.5097360D0/
      Data BETASP ( 7)/     -14.0625210D0/
      Data BETAPP ( 7)/     -20.0438480D0/
      Data ZSPM3  ( 7)/       2.0280940D0/
      Data ZPPM3  ( 7)/       2.3137280D0/
      Data ALPPM3 ( 7)/       2.8305450D0/
      Data EISOLP ( 7)/    -157.6137755D0/
      Data GSSPM3 ( 7)/      11.9047870D0/
      Data GSPPM3 ( 7)/       7.3485650D0/
      Data GPPPM3 ( 7)/      11.7546720D0/
      Data GP2PM3 ( 7)/      10.8072770D0/
      Data HSPPM3 ( 7)/       1.1367130D0/
      Data DDPM3  ( 7)/       0.6577006D0/
      Data QQPM3  ( 7)/       0.5293383D0/
      Data AMPM3  ( 7)/       0.4375151D0/
      Data ADPM3  ( 7)/       0.5030995D0/
      Data AQPM3  ( 7)/       0.7364933D0/
      Data GuesP1 ( 7,1)/     1.5016740D0/
      Data GuesP2 ( 7,1)/     5.9011480D0/
      Data GuesP3 ( 7,1)/     1.7107400D0/
      Data GuesP1 ( 7,2)/    -1.5057720D0/
      Data GuesP2 ( 7,2)/     6.0046580D0/
      Data GuesP3 ( 7,2)/     1.7161490D0/
      Data (GuesP1( 7,J),J=3,10)/8*0.0d0/
      Data (GuesP2( 7,J),J=3,10)/8*0.0d0/
      Data (GuesP3( 7,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  8      OXYGEN
      Data REFPM3( 8)( 1:40)/' O: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3( 8)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 ( 8)/     -86.9930020D0/
      Data UPPPM3 ( 8)/     -71.8795800D0/
      Data BETASP ( 8)/     -45.2026510D0/
      Data BETAPP ( 8)/     -24.7525150D0/
      Data ZSPM3  ( 8)/       3.7965440D0/
      Data ZPPM3  ( 8)/       2.3894020D0/
      Data ALPPM3 ( 8)/       3.2171020D0/
      Data EISOLP ( 8)/    -289.3422065D0/
      Data GSSPM3 ( 8)/      15.7557600D0/
      Data GSPPM3 ( 8)/      10.6211600D0/
      Data GPPPM3 ( 8)/      13.6540160D0/
      Data GP2PM3 ( 8)/      12.4060950D0/
      Data HSPPM3 ( 8)/       0.5938830D0/
      Data DDPM3  ( 8)/       0.4086173D0/
      Data QQPM3  ( 8)/       0.5125738D0/
      Data AMPM3  ( 8)/       0.5790430D0/
      Data ADPM3  ( 8)/       0.5299630D0/
      Data AQPM3  ( 8)/       0.8179630D0/
      Data GuesP1 ( 8,1)/    -1.1311280D0/
      Data GuesP2 ( 8,1)/     6.0024770D0/
      Data GuesP3 ( 8,1)/     1.6073110D0/
      Data GuesP1 ( 8,2)/     1.1378910D0/
      Data GuesP2 ( 8,2)/     5.9505120D0/
      Data GuesP3 ( 8,2)/     1.5983950D0/
      Data (GuesP1( 8,J),J=3,10)/8*0.0d0/
      Data (GuesP2( 8,J),J=3,10)/8*0.0d0/
      Data (GuesP3( 8,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT  9      FLUORINE
      Data REFPM3( 9)( 1:40)/' F: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3( 9)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 ( 9)/    -110.4353030D0/
      Data UPPPM3 ( 9)/    -105.6850470D0/
      Data BETASP ( 9)/     -48.4059390D0/
      Data BETAPP ( 9)/     -27.7446600D0/
      Data ZSPM3  ( 9)/       4.7085550D0/
      Data ZPPM3  ( 9)/       2.4911780D0/
      Data ALPPM3 ( 9)/       3.3589210D0/
      Data EISOLP ( 9)/    -437.5171690D0/
      Data GSSPM3 ( 9)/      10.4966670D0/
      Data GSPPM3 ( 9)/      16.0736890D0/
      Data GPPPM3 ( 9)/      14.8172560D0/
      Data GP2PM3 ( 9)/      14.4183930D0/
      Data HSPPM3 ( 9)/       0.7277630D0/
      Data DDPM3  ( 9)/       0.3125302D0/
      Data QQPM3  ( 9)/       0.4916328D0/
      Data AMPM3  ( 9)/       0.3857650D0/
      Data ADPM3  ( 9)/       0.6768503D0/
      Data AQPM3  ( 9)/       0.6120047D0/
      Data GuesP1 ( 9,1)/    -0.0121660D0/
      Data GuesP2 ( 9,1)/     6.0235740D0/
      Data GuesP3 ( 9,1)/     1.8568590D0/
      Data GuesP1 ( 9,2)/    -0.0028520D0/
      Data GuesP2 ( 9,2)/     6.0037170D0/
      Data GuesP3 ( 9,2)/     2.6361580D0/
      Data (GuesP1( 9,J),J=3,10)/8*0.0d0/
      Data (GuesP2( 9,J),J=3,10)/8*0.0d0/
      Data (GuesP3( 9,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 11
      Data USSPM3 (11)/      -4.7683450D0/
      Data UPPPM3 (11)/      -4.7703090D0/
      Data BETASP (11)/      -0.5354170D0/
      Data BETAPP (11)/       0.0548760D0/
      Data ZSPM3  (11)/       0.6767940D0/
      Data ZPPM3  (11)/       0.9876530D0/
      Data ALPPM3 (11)/       1.681D0/
      Data EISOLP (11)/      -4.7683450D0/
      Data GSSPM3 (11)/       5.1531420D0/
      Data GSPPM3 (11)/       3.1778570D0/
      Data GPPPM3 (11)/       3.7335330D0/
      Data GP2PM3 (11)/      11.5386740D0/
      Data HSPPM3 (11)/       0.5455510D0/
      Data DDPM3  (11)/       2.1443749D0/
      Data QQPM3  (11)/       1.6942388D0/
      Data AMPM3  (11)/       0.1893841D0/
      Data ADPM3  (11)/       0.1850465D0/
      Data AQPM3  (11)/       0.2057464D0/
      Data GuesP1 (11,1)/     0.4175960D0/
      Data GuesP2 (11,1)/     8.5969780D0/
      Data GuesP3 (11,1)/     1.6731780D0/
      Data GuesP1 (11,2)/     0.4432600D0/
      Data GuesP2 (11,2)/     3.3594960D0/
      Data GuesP3 (11,2)/     2.2315110D0/
      Data (GuesP1(11,J),J=3,10)/8*0.0d0/
      Data (GuesP2(11,J),J=3,10)/8*0.0d0/
      Data (GuesP3(11,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 12        MAGNESIUM
      Data REFPM3(12)( 1:40)/'Mg: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(12)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (12)/     -14.6236880D0/
      Data UPPPM3 (12)/     -14.1734600D0/
      Data BETASP (12)/      -2.0716910D0/
      Data BETAPP (12)/      -0.5695810D0/
      Data ZSPM3  (12)/       0.6985520D0/
      Data ZPPM3  (12)/       1.4834530D0/
      Data ALPPM3 (12)/       1.3291470D0/
      Data EISOLP (12)/     -22.5530760D0/
      Data GSSPM3 (12)/       6.6943000D0/
      Data GSPPM3 (12)/       6.7939950D0/
      Data GPPPM3 (12)/       6.9104460D0/
      Data GP2PM3 (12)/       7.0908230D0/
      Data HSPPM3 (12)/       0.5433000D0/
      Data DDPM3  (12)/       1.1403950D0/
      Data QQPM3  (12)/       1.1279899D0/
      Data AMPM3  (12)/       0.2460235D0/
      Data ADPM3  (12)/       0.2695751D0/
      Data AQPM3  (12)/       0.2767522D0/
      Data GuesP1 (12,1)/     2.1170500D0/
      Data GuesP2 (12,1)/     6.0094770D0/
      Data GuesP3 (12,1)/     2.0844060D0/
      Data GuesP1 (12,2)/    -2.5477670D0/
      Data GuesP2 (12,2)/     4.3953700D0/
      Data GuesP3 (12,2)/     2.0636740D0/
      Data (GuesP1(12,J),J=3,10)/8*0.0d0/
      Data (GuesP2(12,J),J=3,10)/8*0.0d0/
      Data (GuesP3(12,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 13      ALUMINUM
      Data REFPM3(13)( 1:40)/'Al: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(13)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (13)/     -24.8454040D0/
      Data UPPPM3 (13)/     -22.2641590D0/
      Data BETASP (13)/      -0.5943010D0/
      Data BETAPP (13)/      -0.9565500D0/
      Data ZSPM3  (13)/       1.7028880D0/
      Data ZPPM3  (13)/       1.0736290D0/
      Data ZDPM3  (13)/       1.0000000D0/
      Data ALPPM3 (13)/       1.5217030D0/
      Data EISOLP (13)/     -46.8647630D0/
      Data GSSPM3 (13)/       5.7767370D0/
      Data GSPPM3 (13)/      11.6598560D0/
      Data GPPPM3 (13)/       6.3477900D0/
      Data GP2PM3 (13)/       6.1210770D0/
      Data HSPPM3 (13)/       4.0062450D0/
      Data DDPM3  (13)/       1.2102799D0/
      Data QQPM3  (13)/       1.5585645D0/
      Data AMPM3  (13)/       0.2123020D0/
      Data ADPM3  (13)/       0.6418584D0/
      Data AQPM3  (13)/       0.2262838D0/
      Data GuesP1 (13,1)/    -0.4730900D0/
      Data GuesP2 (13,1)/     1.9158250D0/
      Data GuesP3 (13,1)/     1.4517280D0/
      Data GuesP1 (13,2)/    -0.1540510D0/
      Data GuesP2 (13,2)/     6.0050860D0/
      Data GuesP3 (13,2)/     2.5199970D0/
      Data (GuesP1(13,J),J=3,10)/8*0.0d0/
      Data (GuesP2(13,J),J=3,10)/8*0.0d0/
      Data (GuesP3(13,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 14      SILICON
      Data REFPM3(14)( 1:40)/'Si: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(14)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (14)/     -26.7634830D0/
      Data UPPPM3 (14)/     -22.8136350D0/
      Data BETASP (14)/      -2.8621450D0/
      Data BETAPP (14)/      -3.9331480D0/
      Data ZSPM3  (14)/       1.6350750D0/
      Data ZPPM3  (14)/       1.3130880D0/
      Data ZDPM3  (14)/       1.0000000D0/
      Data ALPPM3 (14)/       2.1358090D0/
      Data EISOLP (14)/     -67.7882140D0/
      Data GSSPM3 (14)/       5.0471960D0/
      Data GSPPM3 (14)/       5.9490570D0/
      Data GPPPM3 (14)/       6.7593670D0/
      Data GP2PM3 (14)/       5.1612970D0/
      Data HSPPM3 (14)/       0.9198320D0/
      Data DDPM3  (14)/       1.3144550D0/
      Data QQPM3  (14)/       1.2743396D0/
      Data AMPM3  (14)/       0.1854905D0/
      Data ADPM3  (14)/       0.3060715D0/
      Data AQPM3  (14)/       0.4877432D0/
      Data GuesP1 (14,1)/    -0.3906000D0/
      Data GuesP2 (14,1)/     6.0000540D0/
      Data GuesP3 (14,1)/     0.6322620D0/
      Data GuesP1 (14,2)/     0.0572590D0/
      Data GuesP2 (14,2)/     6.0071830D0/
      Data GuesP3 (14,2)/     2.0199870D0/
      Data (GuesP1(14,J),J=3,10)/8*0.0d0/
      Data (GuesP2(14,J),J=3,10)/8*0.0d0/
      Data (GuesP3(14,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 15      PHOSPHORUS
      Data REFPM3(15)( 1:40)/' P: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(15)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (15)/     -40.4130960D0/
      Data UPPPM3 (15)/     -29.5930520D0/
      Data BETASP (15)/     -12.6158790D0/
      Data BETAPP (15)/      -4.1600400D0/
      Data ZSPM3  (15)/       2.0175630D0/
      Data ZPPM3  (15)/       1.5047320D0/
      Data ZDPM3  (15)/       1.0000000D0/
      Data ALPPM3 (15)/       1.9405340D0/
      Data EISOLP (15)/    -117.9591740D0/
      Data GSSPM3 (15)/       7.8016150D0/
      Data GSPPM3 (15)/       5.1869490D0/
      Data GPPPM3 (15)/       6.6184780D0/
      Data GP2PM3 (15)/       6.0620020D0/
      Data HSPPM3 (15)/       1.5428090D0/
      Data DDPM3  (15)/       1.0644947D0/
      Data QQPM3  (15)/       1.1120386D0/
      Data AMPM3  (15)/       0.2867187D0/
      Data ADPM3  (15)/       0.4309446D0/
      Data AQPM3  (15)/       0.3732517D0/
      Data GuesP1 (15,1)/    -0.6114210D0/
      Data GuesP2 (15,1)/     1.9972720D0/
      Data GuesP3 (15,1)/     0.7946240D0/
      Data GuesP1 (15,2)/    -0.0939350D0/
      Data GuesP2 (15,2)/     1.9983600D0/
      Data GuesP3 (15,2)/     1.9106770D0/
      Data (GuesP1(15,J),J=3,10)/8*0.0d0/
      Data (GuesP2(15,J),J=3,10)/8*0.0d0/
      Data (GuesP3(15,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 16      SULFUR
      Data REFPM3(16)( 1:40)/' S: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(16)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (16)/     -49.8953710D0/
      Data UPPPM3 (16)/     -44.3925830D0/
      Data BETASP (16)/      -8.8274650D0/
      Data BETAPP (16)/      -8.0914150D0/
      Data ZSPM3  (16)/       1.8911850D0/
      Data ZPPM3  (16)/       1.6589720D0/
      Data ZDPM3  (16)/       1.0000000D0/
      Data ALPPM3 (16)/       2.2697060D0/
      Data EISOLP (16)/    -183.4537395D0/
      Data GSSPM3 (16)/       8.9646670D0/
      Data GSPPM3 (16)/       6.7859360D0/
      Data GPPPM3 (16)/       9.9681640D0/
      Data GP2PM3 (16)/       7.9702470D0/
      Data HSPPM3 (16)/       4.0418360D0/
      Data DDPM3  (16)/       1.1214313D0/
      Data QQPM3  (16)/       1.0086488D0/
      Data AMPM3  (16)/       0.3294622D0/
      Data ADPM3  (16)/       0.6679118D0/
      Data AQPM3  (16)/       0.6137472D0/
      Data GuesP1 (16,1)/    -0.3991910D0/
      Data GuesP2 (16,1)/     6.0006690D0/
      Data GuesP3 (16,1)/     0.9621230D0/
      Data GuesP1 (16,2)/    -0.0548990D0/
      Data GuesP2 (16,2)/     6.0018450D0/
      Data GuesP3 (16,2)/     1.5799440D0/
      Data (GuesP1(16,J),J=3,10)/8*0.0d0/
      Data (GuesP2(16,J),J=3,10)/8*0.0d0/
      Data (GuesP3(16,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 17      CHLORINE
      Data REFPM3(17)( 1:40)/'Cl: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(17)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (17)/    -100.6267470D0/
      Data UPPPM3 (17)/     -53.6143960D0/
      Data BETASP (17)/     -27.5285600D0/
      Data BETAPP (17)/     -11.5939220D0/
      Data ZSPM3  (17)/       2.2462100D0/
      Data ZPPM3  (17)/       2.1510100D0/
      Data ZDPM3  (17)/       1.0000000D0/
      Data ALPPM3 (17)/       2.5172960D0/
      Data EISOLP (17)/    -315.1949480D0/
      Data GSSPM3 (17)/      16.0136010D0/
      Data GSPPM3 (17)/       8.0481150D0/
      Data GPPPM3 (17)/       7.5222150D0/
      Data GP2PM3 (17)/       7.5041540D0/
      Data HSPPM3 (17)/       3.4811530D0/
      Data DDPM3  (17)/       0.9175856D0/
      Data QQPM3  (17)/       0.7779230D0/
      Data AMPM3  (17)/       0.5885190D0/
      Data ADPM3  (17)/       0.6814522D0/
      Data AQPM3  (17)/       0.3643694D0/
      Data GuesP1 (17,1)/    -0.1715910D0/
      Data GuesP2 (17,1)/     6.0008020D0/
      Data GuesP3 (17,1)/     1.0875020D0/
      Data GuesP1 (17,2)/    -0.0134580D0/
      Data GuesP2 (17,2)/     1.9666180D0/
      Data GuesP3 (17,2)/     2.2928910D0/
      Data (GuesP1(17,J),J=3,10)/8*0.0d0/
      Data (GuesP2(17,J),J=3,10)/8*0.0d0/
      Data (GuesP3(17,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 19
      Data USSPM3 (19)/      -4.1079950D0/
      Data UPPPM3 (19)/      -1.9126300D0/
      Data BETASP (19)/      -0.7235020D0/
      Data BETAPP (19)/       0.0109380D0/
      Data ZSPM3  (19)/       0.9158730D0/
      Data ZPPM3  (19)/       1.3430990D0/
      Data ALPPM3 (19)/       4.4452300D0/
      Data EISOLP (19)/      -4.1079950D0/
      Data GSSPM3 (19)/       4.0502660D0/
      Data GSPPM3 (19)/       2.9879670D0/
      Data GPPPM3 (19)/       2.9983050D0/
      Data GP2PM3 (19)/       4.9936140D0/
      Data HSPPM3 (19)/       0.5653070D0/
      Data DDPM3  (19)/       1.9524869D0/
      Data QQPM3  (19)/       1.5794222D0/
      Data AMPM3  (19)/       0.1488521D0/
      Data ADPM3  (19)/       0.1983782D0/
      Data AQPM3  (19)/       0.2164587D0/
      Data GuesP1 (19,1)/     6.0409170D0/
      Data GuesP2 (19,1)/    10.0465540D0/
      Data GuesP3 (19,1)/     1.6005120D0/
      Data GuesP1 (19,2)/     0.7242240D0/
      Data GuesP2 (19,2)/     3.0082750D0/
      Data GuesP3 (19,2)/     2.4643520D0/
      Data (GuesP1(19,J),J=3,10)/8*0.0d0/
      Data (GuesP2(19,J),J=3,10)/8*0.0d0/
      Data (GuesP3(19,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 20
      Data USSPM3 (20)/     -11.6466770D0/
      Data UPPPM3 (20)/     -11.7003040D0/
      Data BETASP (20)/      -1.3663020D0/
      Data BETAPP (20)/      -0.6800840D0/
      Data ZSPM3  (20)/       0.8224500D0/
      Data ZPPM3  (20)/       1.2791170D0/
      Data ALPPM3 (20)/       3.2251430D0/
      Data EISOLP (20)/     -17.7272050D0/
      Data GSSPM3 (20)/       5.5661490D0/
      Data GSPPM3 (20)/      15.8559960D0/
      Data GPPPM3 (20)/       6.7057840D0/
      Data GP2PM3 (20)/       4.8642700D0/
      Data HSPPM3 (20)/       1.0053200D0/
      Data DDPM3  (20)/       1.9888793D0/
      Data QQPM3  (20)/       1.6584256D0/
      Data AMPM3  (20)/       0.2045626D0/
      Data ADPM3  (20)/       0.2518009D0/
      Data AQPM3  (20)/       0.4412688D0/
      Data GuesP1 (20,1)/     0.7940250D0/
      Data GuesP2 (20,1)/     3.8548560D0/
      Data GuesP3 (20,1)/     1.5284410D0/
      Data GuesP1 (20,2)/     0.3400220D0/
      Data GuesP2 (20,2)/     5.7455560D0/
      Data GuesP3 (20,2)/     2.3571110D0/
      Data (GuesP1(20,J),J=3,10)/8*0.0d0/
      Data (GuesP2(20,J),J=3,10)/8*0.0d0/
      Data (GuesP3(20,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 24
      Data USSPM3 (24)/     -17.5170270D0/
      Data UPPPM3 (24)/     -12.5337290D0/
      Data UDDPM3 (24)/     -44.1249280D0/
      Data BETASP (24)/      -0.1000000D0/
      Data BETAPP (24)/      -0.1000000D0/
      Data BETADP (24)/      -8.7766360D0/
      Data ZSPM3  (24)/       1.5000000D0/
      Data ZPPM3  (24)/       1.5000000D0/
      Data ZDPM3  (24)/       2.8845490D0/
      Data ALPPM3 (24)/       3.0683070D0/
      Data EISOLP (24)/    -134.8187920D0/
      Data GSSPM3 (24)/       6.0000000D0/
      Data GSPPM3 (24)/       4.1500000D0/
      Data GPPPM3 (24)/       5.0000000D0/
      Data GP2PM3 (24)/       3.5000000D0/
      Data HSPPM3 (24)/       1.0000000D0/
      Data DDPM3  (24)/       1.7320508D0/
      Data QQPM3  (24)/       1.4142136D0/
      Data AMPM3  (24)/       0.2205072D0/
      Data ADPM3  (24)/       0.2711332D0/
      Data AQPM3  (24)/       0.4464656D0/
      Data (GuesP1(24,J),J=1,10)/10*0.0d0/
      Data (GuesP2(24,J),J=1,10)/10*0.0d0/
      Data (GuesP3(24,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 30        ZINC
      Data REFPM3(30)( 1:40)/'Zn: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(30)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (30)/     -18.5321980D0/
      Data UPPPM3 (30)/     -11.0474090D0/
      Data BETASP (30)/      -0.7155780D0/
      Data BETAPP (30)/      -6.3518640D0/
      Data ZSPM3  (30)/       1.8199890D0/
      Data ZPPM3  (30)/       1.5069220D0/
      Data ZDPM3  (30)/       1.0000000D0/
      Data ALPPM3 (30)/       1.3501260D0/
      Data EISOLP (30)/     -27.3872000D0/
      Data GSSPM3 (30)/       9.6771960D0/
      Data GSPPM3 (30)/       7.7362040D0/
      Data GPPPM3 (30)/       4.9801740D0/
      Data GP2PM3 (30)/       4.6696560D0/
      Data HSPPM3 (30)/       0.6004130D0/
      Data DDPM3  (30)/       1.5005758D0/
      Data QQPM3  (30)/       1.4077174D0/
      Data AMPM3  (30)/       0.3556485D0/
      Data ADPM3  (30)/       0.2375689D0/
      Data AQPM3  (30)/       0.2661069D0/
      Data GUESP1 (30,1)/    -0.1112340D0/
      Data GUESP2 (30,1)/     6.0014780D0/
      Data GUESP3 (30,1)/     1.5160320D0/
      Data GUESP1 (30,2)/    -0.1323700D0/
      Data GUESP2 (30,2)/     1.9958390D0/
      Data GUESP3 (30,2)/     2.5196420D0/
C                    Data FOR ELEMENT 31        GALLIUM
      Data REFPM3(31)( 1:40)/'Ga: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(31)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (31)/     -29.8555930D0/
      Data UPPPM3 (31)/     -21.8753710D0/
      Data BETASP (31)/      -4.9456180D0/
      Data BETAPP (31)/      -0.4070530D0/
      Data ZSPM3  (31)/       1.8470400D0/
      Data ZPPM3  (31)/       0.8394110D0/
      Data ALPPM3 (31)/       1.6051150D0/
      Data EISOLP (31)/     -57.3280250D0/
      Data GSSPM3 (31)/       8.4585540D0/
      Data GSPPM3 (31)/       8.9256190D0/
      Data GPPPM3 (31)/       5.0868550D0/
      Data GP2PM3 (31)/       4.9830450D0/
      Data HSPPM3 (31)/       2.0512600D0/
      Data DDPM3  (31)/       0.9776692D0/
      Data QQPM3  (31)/       2.5271534D0/
      Data AMPM3  (31)/       0.3108620D0/
      Data ADPM3  (31)/       0.5129360D0/
      Data AQPM3  (31)/       0.1546208D0/
      Data GuesP1 (31,1)/    -0.5601790D0/
      Data GuesP2 (31,1)/     5.6232730D0/
      Data GuesP3 (31,1)/     1.5317800D0/
      Data GuesP1 (31,2)/    -0.2727310D0/
      Data GuesP2 (31,2)/     1.9918430D0/
      Data GuesP3 (31,2)/     2.1838640D0/
      Data (GuesP1(31,J),J=3,10)/8*0.0d0/
      Data (GuesP2(31,J),J=3,10)/8*0.0d0/
      Data (GuesP3(31,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 32        GERMANIUM
      Data REFPM3(32)( 1:40)/'Ge: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(32)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (32)/     -35.4671955D0/
      Data UPPPM3 (32)/     -31.5863583D0/
      Data BETASP (32)/      -5.3250024D0/
      Data BETAPP (32)/      -2.2501567D0/
      Data ZSPM3  (32)/       2.2373526D0/
      Data ZPPM3  (32)/       1.5924319D0/
      Data ALPPM3 (32)/       1.9723370D0/
      Data EISOLP (32)/     -84.0156006D0/
      Data GSSPM3 (32)/       5.3769635D0/
      Data GSPPM3 (32)/      10.2095293D0/
      Data GPPPM3 (32)/       7.6718647D0/
      Data GP2PM3 (32)/       6.9242663D0/
      Data HSPPM3 (32)/       1.3370204D0/
      Data DDPM3  (32)/       1.1920304D0/
      Data QQPM3  (32)/       1.3321263D0/
      Data AMPM3  (32)/       0.1976098D0/
      Data ADPM3  (32)/       0.3798182D0/
      Data AQPM3  (32)/       0.3620669D0/
      Data GuesP1 (32,1)/     0.9631726D0/
      Data GuesP2 (32,1)/     6.0120134D0/
      Data GuesP3 (32,1)/     2.1633655D0/
      Data GuesP1 (32,2)/    -0.9593891D0/
      Data GuesP2 (32,2)/     5.7491802D0/
      Data GuesP3 (32,2)/     2.1693724D0/
      Data (GuesP1(32,J),J=3,10)/8*0.0d0/
      Data (GuesP2(32,J),J=3,10)/8*0.0d0/
      Data (GuesP3(32,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 33        ARSENIC
      Data REFPM3(33)( 1:40)/'As: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(33)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (33)/     -38.5074240D0/
      Data UPPPM3 (33)/     -35.1524150D0/
      Data BETASP (33)/      -8.2321650D0/
      Data BETAPP (33)/      -5.0173860D0/
      Data ZSPM3  (33)/       2.6361770D0/
      Data ZPPM3  (33)/       1.7038890D0/
      Data ALPPM3 (33)/       1.7944770D0/
      Data EISOLP (33)/    -122.6326140D0/
      Data GSSPM3 (33)/       8.7890010D0/
      Data GSPPM3 (33)/       5.3979830D0/
      Data GPPPM3 (33)/       8.2872500D0/
      Data GP2PM3 (33)/       8.2103460D0/
      Data HSPPM3 (33)/       1.9510340D0/
      Data DDPM3  (33)/       0.9679655D0/
      Data QQPM3  (33)/       1.2449874D0/
      Data AMPM3  (33)/       0.3230063D0/
      Data ADPM3  (33)/       0.5042239D0/
      Data AQPM3  (33)/       0.2574219D0/
      Data GuesP1 (33,1)/    -0.4600950D0/
      Data GuesP2 (33,1)/     1.9831150D0/
      Data GuesP3 (33,1)/     1.0867930D0/
      Data GuesP1 (33,2)/    -0.0889960D0/
      Data GuesP2 (33,2)/     1.9929440D0/
      Data GuesP3 (33,2)/     2.1400580D0/
      Data (GuesP1(33,J),J=3,10)/8*0.0d0/
      Data (GuesP2(33,J),J=3,10)/8*0.0d0/
      Data (GuesP3(33,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 34        SELENIUM
      Data REFPM3(34)( 1:40)/'Se: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(34)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (34)/     -55.3781350D0/
      Data UPPPM3 (34)/     -49.8230760D0/
      Data BETASP (34)/      -6.1578220D0/
      Data BETAPP (34)/      -5.4930390D0/
      Data ZSPM3  (34)/       2.8280510D0/
      Data ZPPM3  (34)/       1.7325360D0/
      Data ALPPM3 (34)/       3.0439570D0/
      Data EISOLP (34)/    -192.7748115D0/
      Data GSSPM3 (34)/       7.4325910D0/
      Data GSPPM3 (34)/      10.0604610D0/
      Data GPPPM3 (34)/       9.5683260D0/
      Data GP2PM3 (34)/       7.7242890D0/
      Data HSPPM3 (34)/       4.0165580D0/
      Data DDPM3  (34)/       0.8719813D0/
      Data QQPM3  (34)/       1.2244019D0/
      Data AMPM3  (34)/       0.2731566D0/
      Data ADPM3  (34)/       0.7509697D0/
      Data AQPM3  (34)/       0.5283737D0/
      Data GuesP1 (34,1)/     0.0478730D0/
      Data GuesP2 (34,1)/     6.0074000D0/
      Data GuesP3 (34,1)/     2.0817170D0/
      Data GuesP1 (34,2)/     0.1147200D0/
      Data GuesP2 (34,2)/     6.0086720D0/
      Data GuesP3 (34,2)/     1.5164230D0/
      Data (GuesP1(34,J),J=3,10)/8*0.0d0/
      Data (GuesP2(34,J),J=3,10)/8*0.0d0/
      Data (GuesP3(34,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 35      BROMINE
      Data REFPM3(35)( 1:40)/'Br: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(35)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (35)/    -116.6193110D0/
      Data UPPPM3 (35)/     -74.2271290D0/
      Data BETASP (35)/     -31.1713420D0/
      Data BETAPP (35)/      -6.8140130D0/
      Data ZSPM3  (35)/       5.3484570D0/
      Data ZPPM3  (35)/       2.1275900D0/
      Data ZDPM3  (35)/       1.0000000D0/
      Data ALPPM3 (35)/       2.5118420D0/
      Data EISOLP (35)/    -352.5398970D0/
      Data GSSPM3 (35)/      15.9434250D0/
      Data GSPPM3 (35)/      16.0616800D0/
      Data GPPPM3 (35)/       8.2827630D0/
      Data GP2PM3 (35)/       7.8168490D0/
      Data HSPPM3 (35)/       0.5788690D0/
      Data DDPM3  (35)/       0.2759025D0/
      Data QQPM3  (35)/       0.9970532D0/
      Data AMPM3  (35)/       0.5859399D0/
      Data ADPM3  (35)/       0.6755383D0/
      Data AQPM3  (35)/       0.3823719D0/
      Data GuesP1 (35,1)/     0.9604580D0/
      Data GuesP2 (35,1)/     5.9765080D0/
      Data GuesP3 (35,1)/     2.3216540D0/
      Data GuesP1 (35,2)/    -0.9549160D0/
      Data GuesP2 (35,2)/     5.9447030D0/
      Data GuesP3 (35,2)/     2.3281420D0/
      Data (GuesP1(35,J),J=3,10)/8*0.0d0/
      Data (GuesP2(35,J),J=3,10)/8*0.0d0/
      Data (GuesP3(35,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 37
      Data USSPM3 (37)/      -4.4841970D0/
      Data UPPPM3 (37)/      -3.0722520D0/
      Data BETASP (37)/      -0.7442240D0/
      Data BETAPP (37)/      -0.5677060D0/
      Data ZSPM3  (37)/       0.6801340D0/
      Data ZPPM3  (37)/       0.9658140D0/
      Data ALPPM3 (37)/       5.7668440D0/
      Data EISOLP (37)/      -4.4841970D0/
      Data GSSPM3 (37)/       4.9848710D0/
      Data GSPPM3 (37)/       4.9983590D0/
      Data GPPPM3 (37)/       6.9857170D0/
      Data GP2PM3 (37)/       5.5079160D0/
      Data HSPPM3 (37)/       0.6211430D0/
      Data DDPM3  (37)/       3.2610134D0/
      Data QQPM3  (37)/       2.6599806D0/
      Data AMPM3  (37)/       0.1832000D0/
      Data ADPM3  (37)/       0.1544789D0/
      Data AQPM3  (37)/       0.3073402D0/
      Data GuesP1 (37,1)/     4.7724480D0/
      Data GuesP2 (37,1)/     5.8039230D0/
      Data GuesP3 (37,1)/     1.4442230D0/
      Data GuesP1 (37,2)/     0.5062010D0/
      Data GuesP2 (37,2)/     2.8303820D0/
      Data GuesP3 (37,2)/     2.6433260D0/
      Data (GuesP1(37,J),J=3,10)/8*0.0d0/
      Data (GuesP2(37,J),J=3,10)/8*0.0d0/
      Data (GuesP3(37,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 38
      Data USSPM3 (38)/     -11.4870450D0/
      Data UPPPM3 (38)/      -7.7152160D0/
      Data BETASP (38)/      -0.7533420D0/
      Data BETAPP (38)/      -0.5760720D0/
      Data ZSPM3  (38)/       1.5592340D0/
      Data ZPPM3  (38)/       1.2129820D0/
      Data ALPPM3 (38)/       4.5905860D0/
      Data EISOLP (38)/     -16.9990050D0/
      Data GSSPM3 (38)/       5.9750850D0/
      Data GSPPM3 (38)/       4.8869600D0/
      Data GPPPM3 (38)/       7.3407630D0/
      Data GP2PM3 (38)/      15.8140850D0/
      Data HSPPM3 (38)/       4.0512180D0/
      Data DDPM3  (38)/       2.1011077D0/
      Data QQPM3  (38)/       2.1179593D0/
      Data AMPM3  (38)/       0.2195915D0/
      Data ADPM3  (38)/       0.5137005D0/
      Data AQPM3  (38)/       0.1752812D0/
      Data GuesP1 (38,1)/     0.4398870D0/
      Data GuesP2 (38,1)/     5.7510950D0/
      Data GuesP3 (38,1)/     1.6731940D0/
      Data GuesP1 (38,2)/     0.3957700D0/
      Data GuesP2 (38,2)/     4.4993190D0/
      Data GuesP3 (38,2)/     2.3281550D0/
      Data (GuesP1(38,J),J=3,10)/8*0.0d0/
      Data (GuesP2(38,J),J=3,10)/8*0.0d0/
      Data (GuesP3(38,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 48        CADMIUM
      Data REFPM3(48)( 1:40)/'Cd: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(48)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (48)/     -15.8285840D0/
      Data UPPPM3 (48)/       8.7497950D0/
      Data BETASP (48)/      -8.5819440D0/
      Data BETAPP (48)/      -0.6010340D0/
      Data ZSPM3  (48)/       1.6793510D0/
      Data ZPPM3  (48)/       2.0664120D0/
      Data ALPPM3 (48)/       1.5253820D0/
      Data EISOLP (48)/     -22.4502080D0/
      Data GSSPM3 (48)/       9.2069600D0/
      Data GSPPM3 (48)/       8.2315390D0/
      Data GPPPM3 (48)/       4.9481040D0/
      Data GP2PM3 (48)/       4.6696560D0/
      Data HSPPM3 (48)/       1.6562340D0/
      Data DDPM3  (48)/       1.5982681D0/
      Data QQPM3  (48)/       1.2432402D0/
      Data AMPM3  (48)/       0.3383668D0/
      Data ADPM3  (48)/       0.3570290D0/
      Data AQPM3  (48)/       0.2820582D0/
      Data (GuesP1(48,J),J=1,10)/10*0.0d0/
      Data (GuesP2(48,J),J=1,10)/10*0.0d0/
      Data (GuesP3(48,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 49        INDIUM
      Data REFPM3(49)( 1:40)/'In: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(49)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (49)/     -26.1762050D0/
      Data UPPPM3 (49)/     -20.0058220D0/
      Data BETASP (49)/      -2.9933190D0/
      Data BETAPP (49)/      -1.8289080D0/
      Data ZSPM3  (49)/       2.0161160D0/
      Data ZPPM3  (49)/       1.4453500D0/
      Data ALPPM3 (49)/       1.4183850D0/
      Data EISOLP (49)/     -51.9750470D0/
      Data GSSPM3 (49)/       6.5549000D0/
      Data GSPPM3 (49)/       8.2298730D0/
      Data GPPPM3 (49)/       6.2992690D0/
      Data GP2PM3 (49)/       4.9842110D0/
      Data HSPPM3 (49)/       2.6314610D0/
      Data DDPM3  (49)/       1.5766241D0/
      Data QQPM3  (49)/       1.7774563D0/
      Data AMPM3  (49)/       0.2409004D0/
      Data ADPM3  (49)/       0.4532655D0/
      Data AQPM3  (49)/       0.3689812D0/
      Data GuesP1( 49,1)/    -0.3431380D0/
      Data GuesP2( 49,1)/     1.9940340D0/
      Data GuesP3( 49,1)/     1.6255160D0/
      Data GuesP1( 49,2)/    -0.1095320D0/
      Data GuesP2( 49,2)/     5.6832170D0/
      Data GuesP3( 49,2)/     2.8670090D0/
      Data (GuesP1(49,J),J=3,10)/8*0.0d0/
      Data (GuesP2(49,J),J=3,10)/8*0.0d0/
      Data (GuesP3(49,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 50        TIN
      Data REFPM3(50)( 1:40)/'Sn: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(50)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (50)/     -34.5501920D0/
      Data UPPPM3 (50)/     -25.8944190D0/
      Data BETASP (50)/      -2.7858020D0/
      Data BETAPP (50)/      -2.0059990D0/
      Data ZSPM3  (50)/       2.3733280D0/
      Data ZPPM3  (50)/       1.6382330D0/
      Data ALPPM3 (50)/       1.6996500D0/
      Data EISOLP (50)/     -78.8877790D0/
      Data GSSPM3 (50)/      10.1900330D0/
      Data GSPPM3 (50)/       7.2353270D0/
      Data GPPPM3 (50)/       5.6738100D0/
      Data GP2PM3 (50)/       5.1822140D0/
      Data HSPPM3 (50)/       1.0331570D0/
      Data DDPM3  (50)/       1.3120038D0/
      Data QQPM3  (50)/       1.5681814D0/
      Data AMPM3  (50)/       0.3744959D0/
      Data ADPM3  (50)/       0.3218163D0/
      Data AQPM3  (50)/       0.2832529D0/
      Data GuesP1 (50,1)/    -0.1503530D0/
      Data GuesP2 (50,1)/     6.0056940D0/
      Data GuesP3 (50,1)/     1.7046420D0/
      Data GuesP1 (50,2)/    -0.0444170D0/
      Data GuesP2 (50,2)/     2.2573810D0/
      Data GuesP3 (50,2)/     2.4698690D0/
      Data (GuesP1(50,J),J=3,10)/8*0.0d0/
      Data (GuesP2(50,J),J=3,10)/8*0.0d0/
      Data (GuesP3(50,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 51        ANTIMONY
      Data REFPM3(51)( 1:40)/'Sb: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(51)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (51)/     -56.4321960D0/
      Data UPPPM3 (51)/     -29.4349540D0/
      Data BETASP (51)/     -14.7942170D0/
      Data BETAPP (51)/      -2.8179480D0/
      Data ZSPM3  (51)/       2.3430390D0/
      Data ZPPM3  (51)/       1.8999920D0/
      Data ALPPM3 (51)/       2.0343010D0/
      Data EISOLP (51)/    -148.9382890D0/
      Data GSSPM3 (51)/       9.2382770D0/
      Data GSPPM3 (51)/       5.2776800D0/
      Data GPPPM3 (51)/       6.3500000D0/
      Data GP2PM3 (51)/       6.2500000D0/
      Data HSPPM3 (51)/       2.4244640D0/
      Data DDPM3  (51)/       1.4091903D0/
      Data QQPM3  (51)/       1.3521354D0/
      Data AMPM3  (51)/       0.3395177D0/
      Data ADPM3  (51)/       0.4589010D0/
      Data AQPM3  (51)/       0.2423472D0/
      Data GuesP1 (51,1)/     3.0020280D0/
      Data GuesP2 (51,1)/     6.0053420D0/
      Data GuesP3 (51,1)/     0.8530600D0/
      Data GuesP1 (51,2)/    -0.0188920D0/
      Data GuesP2 (51,2)/     6.0114780D0/
      Data GuesP3 (51,2)/     2.7933110D0/
      Data (GuesP1(51,J),J=3,10)/8*0.0d0/
      Data (GuesP2(51,J),J=3,10)/8*0.0d0/
      Data (GuesP3(51,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 52        TELLURIUM
      Data REFPM3(52)( 1:40)/'Te: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(52)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (52)/     -44.9380360D0/
      Data UPPPM3 (52)/     -46.3140990D0/
      Data BETASP (52)/      -2.6651460D0/
      Data BETAPP (52)/      -3.8954300D0/
      Data ZSPM3  (52)/       4.1654920D0/
      Data ZPPM3  (52)/       1.6475550D0/
      Data ALPPM3 (52)/       2.4850190D0/
      Data EISOLP (52)/    -168.0945925D0/
      Data GSSPM3 (52)/      10.2550730D0/
      Data GSPPM3 (52)/       8.1691450D0/
      Data GPPPM3 (52)/       7.7775920D0/
      Data GP2PM3 (52)/       7.7551210D0/
      Data HSPPM3 (52)/       3.7724620D0/
      Data DDPM3  (52)/       0.3484177D0/
      Data QQPM3  (52)/       1.5593085D0/
      Data AMPM3  (52)/       0.3768862D0/
      Data ADPM3  (52)/       1.1960743D0/
      Data AQPM3  (52)/       0.2184786D0/
      Data GuesP1 (52,1)/     0.0333910D0/
      Data GuesP2 (52,1)/     5.9563790D0/
      Data GuesP3 (52,1)/     2.2775750D0/
      Data GuesP1 (52,2)/    -1.9218670D0/
      Data GuesP2 (52,2)/     4.9732190D0/
      Data GuesP3 (52,2)/     0.5242430D0/
      Data (GuesP1(52,J),J=3,10)/8*0.0d0/
      Data (GuesP2(52,J),J=3,10)/8*0.0d0/
      Data (GuesP3(52,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 53      IODINE
      Data REFPM3(53)( 1:40)/' I: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(53)(41:80)/'EM. 10, 209 (1989).'/
      Data USSPM3 (53)/     -96.4540370D0/
      Data UPPPM3 (53)/     -61.0915820D0/
      Data BETASP (53)/     -14.4942340D0/
      Data BETAPP (53)/      -5.8947030D0/
      Data ZSPM3  (53)/       7.0010130D0/
      Data ZPPM3  (53)/       2.4543540D0/
      Data ZDPM3  (53)/       1.0000000D0/
      Data ALPPM3 (53)/       1.9901850D0/
      Data EISOLP (53)/    -288.3160860D0/
      Data GSSPM3 (53)/      13.6319430D0/
      Data GSPPM3 (53)/      14.9904060D0/
      Data GPPPM3 (53)/       7.2883300D0/
      Data GP2PM3 (53)/       5.9664070D0/
      Data HSPPM3 (53)/       2.6300350D0/
      Data DDPM3  (53)/       0.1581469D0/
      Data QQPM3  (53)/       1.0467302D0/
      Data AMPM3  (53)/       0.5009902D0/
      Data ADPM3  (53)/       1.6699104D0/
      Data AQPM3  (53)/       0.5153082D0/
      Data GuesP1 (53,1)/    -0.1314810D0/
      Data GuesP2 (53,1)/     5.2064170D0/
      Data GuesP3 (53,1)/     1.7488240D0/
      Data GuesP1 (53,2)/    -0.0368970D0/
      Data GuesP2 (53,2)/     6.0101170D0/
      Data GuesP3 (53,2)/     2.7103730D0/
      Data (GuesP1(53,J),J=3,10)/8*0.0d0/
      Data (GuesP2(53,J),J=3,10)/8*0.0d0/
      Data (GuesP3(53,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 55
      Data USSPM3 (55)/      -4.3207390D0/
      Data UPPPM3 (55)/      -3.3475800D0/
      Data BETASP (55)/      -0.9775130D0/
      Data BETAPP (55)/      -0.5484990D0/
      Data ZSPM3  (55)/       0.9188020D0/
      Data ZPPM3  (55)/       1.5522530D0/
      Data ALPPM3 (55)/       2.5439870D0/
      Data EISOLP (55)/      -4.3207390D0/
      Data GSSPM3 (55)/       4.9615320D0/
      Data GSPPM3 (55)/       4.9943210D0/
      Data GPPPM3 (55)/       5.0115170D0/
      Data GP2PM3 (55)/      16.0009950D0/
      Data HSPPM3 (55)/       0.5830050D0/
      Data DDPM3  (55)/       1.9526201D0/
      Data QQPM3  (55)/       1.9433821D0/
      Data AMPM3  (55)/       0.1823422D0/
      Data ADPM3  (55)/       0.2009234D0/
      Data AQPM3  (55)/       0.1864057D0/
      Data GuesP1 (55,1)/     3.3960620D0/
      Data GuesP2 (55,1)/     6.0007540D0/
      Data GuesP3 (55,1)/     1.6762910D0/
      Data GuesP1 (55,2)/     0.6677680D0/
      Data GuesP2 (55,2)/     3.0084010D0/
      Data GuesP3 (55,2)/     2.7188720D0/
      Data (GuesP1(55,J),J=3,10)/8*0.0d0/
      Data (GuesP2(55,J),J=3,10)/8*0.0d0/
      Data (GuesP3(55,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 56
      Data USSPM3 (56)/     -13.5286010D0/
      Data UPPPM3 (56)/     -10.1346170D0/
      Data BETASP (56)/      -3.9502870D0/
      Data BETAPP (56)/      -0.5749730D0/
      Data ZSPM3  (56)/       1.7594520D0/
      Data ZPPM3  (56)/       1.4916880D0/
      Data ALPPM3 (56)/       1.6584300D0/
      Data EISOLP (56)/     -16.9427990D0/
      Data GSSPM3 (56)/      10.1144030D0/
      Data GSPPM3 (56)/       4.8886780D0/
      Data GPPPM3 (56)/       8.8650930D0/
      Data GP2PM3 (56)/       8.7503680D0/
      Data HSPPM3 (56)/       1.4151370D0/
      Data DDPM3  (56)/       2.2086835D0/
      Data QQPM3  (56)/       2.0222866D0/
      Data AMPM3  (56)/       0.3717164D0/
      Data ADPM3  (56)/       0.2800750D0/
      Data AQPM3  (56)/       0.1811672D0/
      Data GuesP1 (56,1)/     0.3501710D0/
      Data GuesP2 (56,1)/     5.8729050D0/
      Data GuesP3 (56,1)/     1.9750730D0/
      Data GuesP1 (56,2)/     0.3339000D0/
      Data GuesP2 (56,2)/     4.0328190D0/
      Data GuesP3 (56,2)/     2.6276060D0/
      Data (GuesP1(56,J),J=3,10)/8*0.0d0/
      Data (GuesP2(56,J),J=3,10)/8*0.0d0/
      Data (GuesP3(56,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 80        MERCURY
      Data REFPM3(80)( 1:40)/'Hg: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(80)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (80)/     -17.7622290D0/
      Data UPPPM3 (80)/     -18.3307510D0/
      Data BETASP (80)/      -3.1013650D0/
      Data BETAPP (80)/      -3.4640310D0/
      Data ZSPM3  (80)/       1.4768850D0/
      Data ZPPM3  (80)/       2.4799510D0/
      Data ALPPM3 (80)/       1.5293770D0/
      Data EISOLP (80)/     -28.8997380D0/
      Data GSSPM3 (80)/       6.6247200D0/
      Data GSPPM3 (80)/      10.6392970D0/
      Data GPPPM3 (80)/      14.7092830D0/
      Data GP2PM3 (80)/      16.0007400D0/
      Data HSPPM3 (80)/       2.0363110D0/
      Data DDPM3  (80)/       1.2317811D0/
      Data QQPM3  (80)/       1.2164033D0/
      Data AMPM3  (80)/       0.2434664D0/
      Data ADPM3  (80)/       0.4515472D0/
      Data AQPM3  (80)/       0.2618394D0/
      Data GuesP1 (80,1)/     1.0827200D0/
      Data GuesP2 (80,1)/     6.4965980D0/
      Data GuesP3 (80,1)/     1.1951460D0/
      Data GuesP1 (80,2)/    -0.0965530D0/
      Data GuesP2 (80,2)/     3.9262810D0/
      Data GuesP3 (80,2)/     2.6271600D0/
      Data (GuesP1(80,J),J=3,10)/8*0.0d0/
      Data (GuesP2(80,J),J=3,10)/8*0.0d0/
      Data (GuesP3(80,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 81        THALLIUM
      Data REFPM3(81)( 1:40)/'Tl: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(81)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (81)/     -30.0531700D0/
      Data UPPPM3 (81)/     -26.9206370D0/
      Data BETASP (81)/      -1.0844950D0/
      Data BETAPP (81)/      -7.9467990D0/
      Data ZSPM3  (81)/       6.8679210D0/
      Data ZPPM3  (81)/       1.9694450D0/
      Data ALPPM3 (81)/       1.3409510D0/
      Data EISOLP (81)/     -56.6492050D0/
      Data GSSPM3 (81)/      10.4604120D0/
      Data GSPPM3 (81)/      11.2238830D0/
      Data GPPPM3 (81)/       4.9927850D0/
      Data GP2PM3 (81)/       8.9627270D0/
      Data HSPPM3 (81)/       2.5304060D0/
      Data DDPM3  (81)/       0.0781362D0/
      Data QQPM3  (81)/       1.5317110D0/
      Data AMPM3  (81)/       0.3844326D0/
      Data ADPM3  (81)/       2.5741815D0/
      Data AQPM3  (81)/       0.2213264D0/
      Data GuesP1 (81,1)/    -1.3613990D0/
      Data GuesP2 (81,1)/     3.5572260D0/
      Data GuesP3 (81,1)/     1.0928020D0/
      Data GuesP1 (81,2)/    -0.0454010D0/
      Data GuesP2 (81,2)/     2.3069950D0/
      Data GuesP3 (81,2)/     2.9650290D0/
      Data (GuesP1(81,J),J=3,10)/8*0.0d0/
      Data (GuesP2(81,J),J=3,10)/8*0.0d0/
      Data (GuesP3(81,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 82        LEAD
      Data REFPM3(82)( 1:40)/'Pb: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(82)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (82)/     -30.3227560D0/
      Data UPPPM3 (82)/     -24.4258340D0/
      Data BETASP (82)/      -6.1260240D0/
      Data BETAPP (82)/      -1.3954300D0/
      Data ZSPM3  (82)/       3.1412890D0/
      Data ZPPM3  (82)/       1.8924180D0/
      Data ALPPM3 (82)/       1.6200450D0/
      Data EISOLP (82)/     -73.4660775D0/
      Data GSSPM3 (82)/       7.0119920D0/
      Data GSPPM3 (82)/       6.7937820D0/
      Data GPPPM3 (82)/       5.1837800D0/
      Data GP2PM3 (82)/       5.0456510D0/
      Data HSPPM3 (82)/       1.5663020D0/
      Data DDPM3  (82)/       0.9866290D0/
      Data QQPM3  (82)/       1.5940562D0/
      Data AMPM3  (82)/       0.2576991D0/
      Data ADPM3  (82)/       0.4527678D0/
      Data AQPM3  (82)/       0.2150175D0/
      Data GuesP1 (82,1)/    -0.1225760D0/
      Data GuesP2 (82,1)/     6.0030620D0/
      Data GuesP3 (82,1)/     1.9015970D0/
      Data GuesP1 (82,2)/    -0.0566480D0/
      Data GuesP2 (82,2)/     4.7437050D0/
      Data GuesP3 (82,2)/     2.8618790D0/
      Data (GuesP1(82,J),J=3,10)/8*0.0d0/
      Data (GuesP2(82,J),J=3,10)/8*0.0d0/
      Data (GuesP3(82,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 83        BISMUTH
      Data REFPM3(83)( 1:40)/'Bi: (PM3): J. J. P. STEWART, J. COMP. CH'/
      Data REFPM3(83)(41:80)/'EM.(ACCEPTED)'/
      Data USSPM3 (83)/     -33.4959380D0/
      Data UPPPM3 (83)/     -35.5210260D0/
      Data BETASP (83)/      -5.6072830D0/
      Data BETAPP (83)/      -5.8001520D0/
      Data ZSPM3  (83)/       4.9164510D0/
      Data ZPPM3  (83)/       1.9349350D0/
      Data ALPPM3 (83)/       1.8574310D0/
      Data EISOLP (83)/    -109.2774910D0/
      Data GSSPM3 (83)/       4.9894800D0/
      Data GSPPM3 (83)/       6.1033080D0/
      Data GPPPM3 (83)/       8.6960070D0/
      Data GP2PM3 (83)/       8.3354470D0/
      Data HSPPM3 (83)/       0.5991220D0/
      Data DDPM3  (83)/       0.2798609D0/
      Data QQPM3  (83)/       1.5590294D0/
      Data AMPM3  (83)/       0.1833693D0/
      Data ADPM3  (83)/       0.6776013D0/
      Data AQPM3  (83)/       0.2586520D0/
      Data GuesP1 (83,1)/     2.5816930D0/
      Data GuesP2 (83,1)/     5.0940220D0/
      Data GuesP3 (83,1)/     0.4997870D0/
      Data GuesP1 (83,2)/     0.0603200D0/
      Data GuesP2 (83,2)/     6.0015380D0/
      Data GuesP3 (83,2)/     2.4279700D0/
      Data (GuesP1(83,J),J=3,10)/8*0.0d0/
      Data (GuesP2(83,J),J=3,10)/8*0.0d0/
      Data (GuesP3(83,J),J=3,10)/8*0.0d0/
C                    Data FOR ELEMENT 90
      Data USSPM3 (90)/     -40.5682920D0/
      Data UPPPM3 (90)/     -28.0891870D0/
      Data BETASP (90)/      -4.2562180D0/
      Data BETAPP (90)/      -4.2562180D0/
      Data ZSPM3  (90)/       1.4353060D0/
      Data ZPPM3  (90)/       1.4353060D0/
      Data ZDPM3  (90)/       1.0000000D0/
      Data ALPPM3 (90)/       2.1961078D0/
      Data EISOLP (90)/       0.0000000D0/
      Data GSSPM3 (90)/       9.8200000D0/
      Data GSPPM3 (90)/       8.3600000D0/
      Data GPPPM3 (90)/       7.3100000D0/
      Data GP2PM3 (90)/       6.5400000D0/
      Data HSPPM3 (90)/       1.3200000D0/
      Data DDPM3  (90)/       0.2011245D0/
      Data QQPM3  (90)/       0.2203208D0/
      Data AMPM3  (90)/       0.3608967D0/
      Data ADPM3  (90)/       1.1124287D0/
      Data AQPM3  (90)/       1.3153602D0/
      Data (GuesP1(90,J),J=1,10)/10*0.0d0/
      Data (GuesP2(90,J),J=1,10)/10*0.0d0/
      Data (GuesP3(90,J),J=1,10)/10*0.0d0/
C                    Data FOR ELEMENT 91
      Data USSPM3 (91)/     -75.2391520D0/
      Data UPPPM3 (91)/     -57.8320130D0/
      Data BETASP (91)/     -11.1422310D0/
      Data BETAPP (91)/     -11.1422310D0/
      Data ZSPM3  (91)/       2.6135910D0/
      Data ZPPM3  (91)/       2.0343930D0/
      Data ZDPM3  (91)/       1.0000000D0/
      Data ALPPM3 (91)/       2.4916445D0/
      Data EISOLP (91)/       0.0000000D0/
      Data HSPPM3 (91)/       0.0000001D0/
      Data DDPM3  (91)/       0.1232470D0/
      Data QQPM3  (91)/       0.1554408D0/
      Data AMPM3  (91)/       0.0000000D0/
      Data ADPM3  (91)/       0.0062312D0/
      Data AQPM3  (91)/       1.2399130D0/
      Data (GuesP1(91,J),J=1,10)/10*0.0d0/
      Data (GuesP2(91,J),J=1,10)/10*0.0d0/
      Data (GuesP3(91,J),J=1,10)/10*0.0d0/
C
      Data REFPM3(102)( 1:40)/'Cb: (PM3):  Capped Bond  (Hydrogen-like'/
      Data REFPM3(102)(41:80)/', takes on a  zero charge.)'/
      Data USSPM3(102)/     -11.9062760D0/
      Data BETASP(102)/-9999999.0000000D0/
      Data ZSPM3 (102)/       4.0000000D0/
      Data ZPPM3 (102)/       0.3000000D0/
      Data ZDPM3 (102)/       0.3000000D0/
      Data ALPPM3(102)/       2.5441341D0/
      Data EISOLP(102)/       4.0000000D0/
      Data GSSPM3(102)/      12.8480000D0/
      Data HSPPM3(102)/       0.1000000D0/
      Data DDPM3 (102)/       0.0684105D0/
      Data QQPM3 (102)/       1.0540926D0/
      Data AMPM3 (102)/       0.4721793D0/
      Data ADPM3 (102)/       0.9262742D0/
      Data AQPM3 (102)/       0.2909059D0/
      Data (GuesP1(102,J),J=1,10)/10*0.0d0/
      Data (GuesP2(102,J),J=1,10)/10*0.0d0/
      Data (GuesP3(102,J),J=1,10)/10*0.0d0/
      End
