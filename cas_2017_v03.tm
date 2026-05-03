KPL/MK

   This meta-kernel lists the CASSINI SPICE kernels providing coverage
   for 2017. All of the kernels listed below are archived in the CAS
   SPICE data set (DATA_SET_ID = "CO-S/J/E/V-SPICE-6-V1.0"). This set
   of files and the order in which they are listed were picked to
   provide the best available data and the most complete coverage for
   the specified year based on the information about the kernels
   available at the time this meta-kernel was made. For detailed
   information about the kernels listed below refer to the internal
   comments included in the kernels and the documentation accompanying
   the CASSINI SPICE data set.

   WARNING: the list of kernels in this meta-kernel includes the latest
   generic PCK file available at the time when this meta-kernel was
   created. This PCK file does not provide all parameters included in
   the CASSINI-specific PCKs and, for some parameters present in this
   file and CASSINI-specific PCKs, provides different values compared
   to the values in the CASSINI-specific PCKs. It is the user's
   responsibility to identify whether a CASSINI-specific PCK(s) is(are)
   needed for his/her application and to use this(these) PCK(s) in
   addition to or instead of the generic PCK file in this meta-kernel.

   It is recommended that users make a local copy of this file and
   modify the value of the PATH_VALUES keyword to point to the actual
   location of the CAS SPICE data set's ``data'' directory on their
   system. Replacing ``/'' with ``\'' and converting line terminators
   to the format native to the user's system may also be required if
   this meta-kernel is to be used on a non-UNIX workstation.

   This file was created on September 27, 2018 by Boris Semenov,
   NAIF/JPL. The original name of this file was cas_2017_v03.tm.

   \begindata

      PATH_VALUES     = ( '_Cassini_kernels' )

      PATH_SYMBOLS    = ( 'KERNELS' )

      KERNELS_TO_LOAD = (

                          '$KERNELS/lsk/naif0012.tls'

                          '$KERNELS/pck/pck00010.tpc'

                          '$KERNELS/fk/cas_rocks_v18.tf'

                          '$KERNELS/fk/cas_mimi_v202.tf'

                          '$KERNELS/fk/cas_dyn_v03.tf'

                          '$KERNELS/fk/cas_v41.tf'

                          '$KERNELS/ik/cas_caps_v03.ti'
                          '$KERNELS/ik/cas_cda_v01.ti'
                          '$KERNELS/ik/cas_cirs_v09.ti'
                          '$KERNELS/ik/cas_inms_v02.ti'
                          '$KERNELS/ik/cas_iss_v10.ti'
                          '$KERNELS/ik/cas_mag_v01.ti'
                          '$KERNELS/ik/cas_mimi_v11.ti'
                          '$KERNELS/ik/cas_radar_v11.ti'
                          '$KERNELS/ik/cas_rpws_v01.ti'
                          '$KERNELS/ik/cas_rss_v03.ti'
                          '$KERNELS/ik/cas_sru_v02.ti'
                          '$KERNELS/ik/cas_uvis_v06.ti'
                          '$KERNELS/ik/cas_vims_v06.ti'

                          '$KERNELS/sclk/cas00172.tsc'

                          '$KERNELS/spk/180927AP_RE_90165_18018.bsp'
                          '$KERNELS/spk/140809BP_IRRE_00256_25017.bsp'

                          '$KERNELS/spk/170509R_SCPSE_16363_17061.bsp'
                          '$KERNELS/spk/170720R_SCPSE_17061_17104.bsp'
                          '$KERNELS/spk/170811R_SCPSE_17104_17117.bsp'
                          '$KERNELS/spk/171017R_SCPSE_17117_17146.bsp'
                          '$KERNELS/spk/171024R_SCPSE_17146_17177.bsp'
                          '$KERNELS/spk/171031R_SCPSE_17177_17224.bsp'
                          '$KERNELS/spk/171214R_SCPSE_17224_17258.bsp'

                          '$KERNELS/ck/17001_18001pa_gapfill_v14.bc'

                          '$KERNELS/ck/16328_17034py_as_flown.bc'
                          '$KERNELS/ck/17034_17104py_as_flown.bc'
                          '$KERNELS/ck/17104_17145py_as_flown.bc'
                          '$KERNELS/ck/17145_17191py_as_flown.bc'
                          '$KERNELS/ck/17191_17258py_as_flown.bc'

                          '$KERNELS/ck/16328_16333ra.bc'
                          '$KERNELS/ck/16333_16338ra.bc'
                          '$KERNELS/ck/16338_16343ra.bc'
                          '$KERNELS/ck/16343_16348ra.bc'
                          '$KERNELS/ck/16348_16353ra.bc'
                          '$KERNELS/ck/16353_16358ra.bc'
                          '$KERNELS/ck/16358_16363ra.bc'
                          '$KERNELS/ck/16363_17002ra.bc'
                          '$KERNELS/ck/17002_17007ra.bc'
                          '$KERNELS/ck/17007_17012ra.bc'
                          '$KERNELS/ck/17012_17017ra.bc'
                          '$KERNELS/ck/17017_17022ra.bc'
                          '$KERNELS/ck/17022_17027ra.bc'
                          '$KERNELS/ck/17027_17032ra.bc'
                          '$KERNELS/ck/17032_17037ra.bc'
                          '$KERNELS/ck/17037_17042ra.bc'
                          '$KERNELS/ck/17042_17047ra.bc'
                          '$KERNELS/ck/17047_17052ra.bc'
                          '$KERNELS/ck/17052_17057ra.bc'
                          '$KERNELS/ck/17057_17062ra.bc'
                          '$KERNELS/ck/17062_17067ra.bc'
                          '$KERNELS/ck/17067_17072ra.bc'
                          '$KERNELS/ck/17072_17077ra.bc'
                          '$KERNELS/ck/17077_17082ra.bc'
                          '$KERNELS/ck/17082_17087ra.bc'
                          '$KERNELS/ck/17087_17092ra.bc'
                          '$KERNELS/ck/17092_17097ra.bc'
                          '$KERNELS/ck/17097_17102ra.bc'
                          '$KERNELS/ck/17102_17107ra.bc'
                          '$KERNELS/ck/17107_17112ra.bc'
                          '$KERNELS/ck/17112_17117ra.bc'
                          '$KERNELS/ck/17117_17122ra.bc'
                          '$KERNELS/ck/17122_17127ra.bc'
                          '$KERNELS/ck/17127_17132ra.bc'
                          '$KERNELS/ck/17132_17137ra.bc'
                          '$KERNELS/ck/17137_17142ra.bc'
                          '$KERNELS/ck/17142_17147ra.bc'
                          '$KERNELS/ck/17147_17152ra.bc'
                          '$KERNELS/ck/17152_17157ra.bc'
                          '$KERNELS/ck/17157_17162ra.bc'
                          '$KERNELS/ck/17162_17167ra.bc'
                          '$KERNELS/ck/17167_17172ra.bc'
                          '$KERNELS/ck/17172_17177ra.bc'
                          '$KERNELS/ck/17177_17182ra.bc'
                          '$KERNELS/ck/17182_17187ra.bc'
                          '$KERNELS/ck/17187_17192ra.bc'
                          '$KERNELS/ck/17192_17197ra.bc'
                          '$KERNELS/ck/17197_17202ra.bc'
                          '$KERNELS/ck/17202_17207ra.bc'
                          '$KERNELS/ck/17207_17212ra.bc'
                          '$KERNELS/ck/17212_17217ra.bc'
                          '$KERNELS/ck/17217_17222ra.bc'
                          '$KERNELS/ck/17222_17227ra.bc'
                          '$KERNELS/ck/17227_17232ra.bc'
                          '$KERNELS/ck/17232_17237ra.bc'
                          '$KERNELS/ck/17237_17242ra.bc'
                          '$KERNELS/ck/17242_17247ra.bc'
                          '$KERNELS/ck/17247_17252ra.bc'
                          '$KERNELS/ck/17252_17257ra.bc'
                          '$KERNELS/ck/17257_17262ra.bc'

                          '$KERNELS/ck/cas_cda_20170627.bc'

                          '$KERNELS/ck/cas_lemms_05109_20001_v2.bc'

                        )

   \begintext

End of MK.
