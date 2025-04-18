| You can create observations using any of the filters below. Many of the filters below are based off the list from python-FSPS (https://python-fsps.readthedocs.io/en/latest/filters/), so please give them credit where credit is due.


| In addition to these FSPS filters, we have included the Chinese Space Station Telescope (CSST) and the DESI filters (BASS_g, BASS_r, MzLS_z). 


| If you need to use a filter that is not included in this list, it should be very easy to add it. In this case, we encourage you to open one of the filter files below to better understand their format. To add a new filter, you need to obtain its transmission curve, save it to a file with the appropriate format, and then add it to the pymgal/pymgal/filters directory. You should then be able to call it in the same way as any of the pre-installed filters.


.. list-table::
   :widths: 5 10 25
   :header-rows: 1

   * - Index
     - Filter name
     - Description
   * - 1
     - V
     - Johnson V (from Bessell 1990 via M. Blanton) - this defines the Vega system
   * - 2
     - U
     - Johnson U (from Bessell 1990 via M. Blanton)
   * - 3
     - B
     - Johnson B (from Bessell 1990 via M. Blanton)
   * - 4
     - Buser_B
     - Buser B2 (from BC03)
   * - 5
     - Cousins_R
     - Cousins R (from Bessell 1990 via M. Blanton)
   * - 6
     - Cousins_I
     - Cousins I (from Bessell 1990 via M. Blanton)
   * - 7
     - CFHT_B
     - CFHT B-band (from Blanton's kcorrect)
   * - 8
     - CFHT_R
     - CFHT R-band (from Blanton's kcorrect)
   * - 9
     - CFHT_I
     - CFHT I-band (from Blanton's kcorrect)
   * - 10
     - 2MASS_J
     - 2MASS J filter (total response w/atm)
   * - 11
     - 2MASS_H
     - 2MASS H filter (total response w/atm)
   * - 12
     - 2MASS_Ks
     - 2MASS Ks filter (total response w/atm)
   * - 13
     - SDSS_u
     - SDSS Camera u Response Function, airmass = 1.3 (June 2001)
   * - 14
     - SDSS_g
     - SDSS Camera g Response Function, airmass = 1.3 (June 2001)
   * - 15
     - SDSS_r
     - SDSS Camera r Response Function, airmass = 1.3 (June 2001)
   * - 16
     - SDSS_i
     - SDSS Camera i Response Function, airmass = 1.3 (June 2001)
   * - 17
     - SDSS_z
     - SDSS Camera z Response Function, airmass = 1.3 (June 2001)
   * - 18
     - WFPC2_F255W
     - HST WFPC2 F255W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 19
     - WFPC2_F300W
     - HST WFPC2 F300W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 20
     - WFPC2_F336W
     - HST WFPC2 F336W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 21
     - WFPC2_F439W
     - HST WFPC2 F439W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 22
     - WFPC2_F450W
     - HST WFPC2 F450W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 23
     - WFPC2_F555W
     - HST WFPC2 F555W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 24
     - WFPC2_F606W
     - HST WFPC2 F606W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 25
     - WFPC2_F814W
     - HST WFPC2 F814W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 26
     - WFPC2_F850LP
     - HST WFPC2 F850LP (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 27
     - WFC_ACS_F435W
     - HST ACS F435W  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 28
     - WFC_ACS_F475W
     - HST ACS F475W  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 29
     - WFC_ACS_F555W
     - HST ACS F555W (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 30
     - WFC_ACS_F606W
     - HST ACS F606W  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 31
     - WFC_ACS_F625W
     - HST ACS F625W  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 32
     - WFC_ACS_F775W
     - HST ACS F775W  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 33
     - WFC_ACS_F814W
     - HST ACS F814W  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 34
     - WFC_ACS_F850LP
     - HST ACS F850LP  (http://acs.pha.jhu.edu/instrument/photometry/)
   * - 35
     - WFC3_UVIS_F218W
     - HST WFC3 UVIS F218W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 36
     - WFC3_UVIS_F225W
     - HST WFC3 UVIS F225W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 37
     - WFC3_UVIS_F275W
     - HST WFC3 UVIS F275W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 38
     - WFC3_UVIS_F336W
     - HST WFC3 UVIS F336W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 39
     - WFC3_UVIS_F390W
     - HST WFC3 UVIS F390W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 40
     - WFC3_UVIS_F438W
     - HST WFC3 UVIS F438W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 41
     - WFC3_UVIS_F475W
     - HST WFC3 UVIS F475W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 42
     - WFC3_UVIS_F555W
     - HST WFC3 UVIS F555W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 43
     - WFC3_UVIS_F606W
     - HST WFC3 UVIS F606W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 44
     - WFC3_UVIS_F775W
     - HST WFC3 UVIS F775W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 45
     - WFC3_UVIS_F814W
     - HST WFC3 UVIS F814W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 46
     - WFC3_UVIS_F850LP
     - HST WFC3 UVIS F850LP (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/) Chip #1
   * - 47
     - WFC3_IR_F098M
     - HST WFC3 IR F098M (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)
   * - 48
     - WFC3_IR_F105W
     - HST WFC3 IR F105W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)
   * - 49
     - WFC3_IR_F110W
     - HST WFC3 IR F110W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)
   * - 50
     - WFC3_IR_F125W
     - HST WFC3 IR F125W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)
   * - 51
     - WFC3_IR_F140W
     - HST WFC3 IR F140W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)
   * - 52
     - WFC3_IR_F160W
     - HST WFC3 IR F160W (http://www.stsci.edu/~WFC3/UVIS/SystemThroughput/)
   * - 53
     - IRAC_1
     - Spitzer IRAC Channel 1 (3.6um)
   * - 54
     - IRAC_2
     - Spitzer IRAC Channel 2 (4.5um)
   * - 55
     - IRAC_3
     - Spitzer IRAC Channel 3 (5.8um)
   * - 56
     - IRAC_4
     - Spitzer IRAC Channel 4 (8.0um)
   * - 57
     - ISAAC_Ks
     - ISAAC Ks
   * - 58
     - FORS_V
     - FORS V
   * - 59
     - FORS_R
     - FORS R
   * - 60
     - NICMOS_F110W
     - HST NICMOS F110W
   * - 61
     - NICMOS_F160W
     - HST NICMOS F160W
   * - 62
     - GALEX_FUV
     - GALEX FUV
   * - 63
     - GALEX_NUV
     - GALEX NUV
   * - 64
     - DES_g
     - DES g  (from Huan Lin, for DES camera)
   * - 65
     - DES_r
     - DES r  (from Huan Lin, for DES camera)
   * - 66
     - DES_i
     - DES i  (from Huan Lin, for DES camera)
   * - 67
     - DES_z
     - DES z  (from Huan Lin, for DES camera)
   * - 68
     - DES_Y
     - DES Y  (from Huan Lin, for DES camera)
   * - 69
     - WFCAM_Z
     - WFCAM (UKIRT) Z  (from Hewett et al. 2006, via A. Smith)
   * - 70
     - WFCAM_Y
     - WFCAM (UKIRT) Y  (from Hewett et al. 2006, via A. Smith)
   * - 71
     - WFCAM_J
     - WFCAM (UKIRT) J  (from Hewett et al. 2006, via A. Smith)
   * - 72
     - WFCAM_H
     - WFCAM (UKIRT) H  (from Hewett et al. 2006, via A. Smith)
   * - 73
     - WFCAM_K
     - WFCAM (UKIRT) K  (from Hewett et al. 2006, via A. Smith)
   * - 74
     - Steidel_Un
     - Steidel Un (via A. Shapley; see Steidel et al. 2003)
   * - 75
     - Steidel_G
     - Steidel G  (via A. Shapley; see Steidel et al. 2003)
   * - 76
     - Steidel_Rs
     - Steidel Rs (via A. Shapley; see Steidel et al. 2003)
   * - 77
     - Steidel_I
     - Steidel I  (via A. Shapley; see Steidel et al. 2003)
   * - 78
     - MegaCam_u
     - CFHT MegaCam u (via 3DHST filter list; includes atm transmission)
   * - 79
     - MegaCam_g
     - CFHT MegaCam g (via 3DHST filter list; includes atm transmission)
   * - 80
     - MegaCam_r
     - CFHT MegaCam r (via 3DHST filter list; includes atm transmission)
   * - 81
     - MegaCam_i
     - CFHT MegaCam i (via 3DHST filter list; includes atm transmission)
   * - 82
     - MegaCam_z
     - CFHT MegaCam z (via 3DHST filter list; includes atm transmission)
   * - 83
     - WISE_W1
     - WISE W1, 3.4um (http://www.astro.ucla.edu/~wright/WISE/passbands.html)
   * - 84
     - WISE_W2
     - WISE W2, 4.6um (http://www.astro.ucla.edu/~wright/WISE/passbands.html)
   * - 85
     - WISE_W3
     - WISE W3, 12um (http://www.astro.ucla.edu/~wright/WISE/passbands.html)
   * - 86
     - WISE_W4
     - WISE W4, 22um (http://www.astro.ucla.edu/~wright/WISE/passbands.html)
   * - 87
     - UVOT_W2
     - UVOT W2 (from Erik Hoversten, 2011)
   * - 88
     - UVOT_M2
     - UVOT M2 (from Erik Hoversten, 2011)
   * - 89
     - UVOT_W1
     - UVOT W1 (from Erik Hoversten, 2011)
   * - 90
     - MIPS_24
     - Spitzer MIPS 24um
   * - 91
     - MIPS_70
     - Spitzer MIPS 70um
   * - 92
     - MIPS_160
     - Spitzer MIPS 160um
   * - 93
     - SCUBA_450WB
     - SCUBA 450WB (www.jach.hawaii.edu/JCMT/continuum/background/background.html)
   * - 94
     - SCUBA_850WB
     - SCUBA 850WB (www.jach.hawaii.edu/JCMT/continuum/background/background.html)
   * - 95
     - PACS_70
     - Herschel PACS   70um
   * - 96
     - PACS_100
     - Herschel PACS  100um
   * - 97
     - PACS_160
     - Herschel PACS  160um
   * - 98
     - SPIRE_250
     - Herschel SPIRE 250um
   * - 99
     - SPIRE_350
     - Herschel SPIRE 350um
   * - 100
     - SPIRE_500
     - Herschel SPIRE 500um
   * - 101
     - IRAS_12
     - IRAS 12um
   * - 102
     - IRAS_25
     - IRAS 25um
   * - 103
     - IRAS_60
     - IRAS 60um
   * - 104
     - IRAS_100
     - IRAS 100um
   * - 105
     - Bessell_L
     - Bessell L band  (Bessell & Brett 1988)
   * - 106
     - Bessell_LP
     - Bessell L' band (Bessell & Brett 1988)
   * - 107
     - Bessell_M
     - Bessell M band  (Bessell & Brett 1988)
   * - 108
     - Stromgren_u
     - Stromgren u (Bessell 2011)
   * - 109
     - Stromgren_v
     - Stromgren v (Bessell 2011)
   * - 110
     - Stromgren_b
     - Stromgren b (Bessell 2011)
   * - 111
     - Stromgren_y
     - Stromgren y (Bessell 2011)
   * - 112
     - I1500
     - Idealized 1500A bandpass with 15% bandwidth, FWHM = 225A from M. Dickinson
   * - 113
     - I2300
     - Idealized 2300A bandpass with 15% bandwidth, FWHM = 345A from M. Dickinson
   * - 114
     - I2800
     - Idealized 2800A bandpass with 15% bandwidth, FWHM = 420A from M. Dickinson
   * - 115
     - JWST_F070W
     - JWST F070W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 116
     - JWST_F090W
     - JWST F090W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 117
     - JWST_F115W
     - JWST F115W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 118
     - JWST_F150W
     - JWST F150W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 119
     - JWST_F200W
     - JWST F200W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 120
     - JWST_F277W
     - JWST F277W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 121
     - JWST_F356W
     - JWST F356W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 122
     - JWST_F444W
     - JWST F444W (https://jwst-docs.stsci.edu/jwst-near-infrared-camera/nircam-instrumentation/nircam-filters)
   * - 123
     - NEWFIRM_J1
     - NEWFIRM J1 (via 3DHST filter list)
   * - 124
     - NEWFIRM_J2
     - NEWFIRM J2 (via 3DHST filter list)
   * - 125
     - NEWFIRM_J3
     - NEWFIRM J3 (via 3DHST filter list)
   * - 126
     - NEWFIRM_H1
     - NEWFIRM H1 (via 3DHST filter list)
   * - 127
     - NEWFIRM_H2
     - NEWFIRM H2 (via 3DHST filter list)
   * - 128
     - NEWFIRM_K
     - NEWFIRM K  (via 3DHST filter list)
   * - 129
     - VISTA_Y
     - VISTA VIRCAM Y (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 130
     - VISTA_J
     - VISTA VIRCAM J (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 131
     - VISTA_H
     - VISTA VIRCAM H (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 132
     - VISTA_K
     - VISTA VIRCAM K (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 133
     - SUPRIMECAM_B
     - Subaru Suprime-Cam B (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 134
     - SUPRIMECAM_g
     - Subaru Suprime-Cam g+ (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 135
     - SUPRIMECAM_V
     - Subaru Suprime-Cam V (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 136
     - SUPRIMECAM_r
     - Subaru Suprime-Cam r+ (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 137
     - SUPRIMECAM_i
     - Subaru Suprime-Cam i+ (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 138
     - SUPRIMECAM_z
     - Subaru Suprime-Cam z+ (http://www.astro.caltech.edu/~capak/filters/index.html)
   * - 139
     - PS1_g
     - Pan-STARRS1 g (http://iopscience.iop.org/0004-637X/750/2/99/suppdata/apj425122t3_mrt.txt)
   * - 140
     - PS1_r
     - Pan-STARRS1 r (http://iopscience.iop.org/0004-637X/750/2/99/suppdata/apj425122t3_mrt.txt)
   * - 141
     - PS1_i
     - Pan-STARRS1 i (http://iopscience.iop.org/0004-637X/750/2/99/suppdata/apj425122t3_mrt.txt)
   * - 142
     - PS1_z
     - Pan-STARRS1 z (http://iopscience.iop.org/0004-637X/750/2/99/suppdata/apj425122t3_mrt.txt)
   * - 143
     - PS1_y
     - Pan-STARRS1 y (http://iopscience.iop.org/0004-637X/750/2/99/suppdata/apj425122t3_mrt.txt)
   * - 144
     - LSST_u
     - LSST u (version 1.5 https://github.com/lsst/throughputs/tree/master/baseline via SVO)
   * - 145
     - LSST_g
     - LSST g (version 1.5 https://github.com/lsst/throughputs/tree/master/baseline via SVO)
   * - 146
     - LSST_r
     - LSST r (version 1.5 https://github.com/lsst/throughputs/tree/master/baseline via SVO)
   * - 147
     - LSST_i
     - LSST i (version 1.5 https://github.com/lsst/throughputs/tree/master/baseline via SVO)
   * - 148
     - LSST_z
     - LSST z (version 1.5 https://github.com/lsst/throughputs/tree/master/baseline via SVO)
   * - 149
     - LSST_y
     - LSST y (version 1.5 https://github.com/lsst/throughputs/tree/master/baseline via SVO)
   * - 150
     - Euclid_VIS
     - Euclid VIS (Master Euclid mission database via SVO)
   * - 151
     - Euclid_Y
     - Euclid Y (Master Euclid mission database via SVO)
   * - 152
     - Euclid_J
     - Euclid J (Master Euclid mission database via SVO)
   * - 153
     - Euclid_H
     - Euclid H (Master Euclid mission database via SVO)
   * - 154
     - Roman_F062
     - Roman F062 (https://roman.gsfc.nasa.gov/science/Roman_Reference_Information.html via SVO)
   * - 155
     - Roman_F087
     - Roman F087 (https://roman.gsfc.nasa.gov/science/Roman_Reference_Information.html via SVO)
   * - 156
     - Roman_F106
     - Roman F106 (https://roman.gsfc.nasa.gov/science/Roman_Reference_Information.html via SVO)
   * - 157
     - Roman_F129
     - Roman F129 (https://roman.gsfc.nasa.gov/science/Roman_Reference_Information.html via SVO)
   * - 158
     - Roman_F158
     - Roman F158 (https://roman.gsfc.nasa.gov/science/Roman_Reference_Information.html via SVO)
   * - 159
     - Roman_F184
     - Roman F184 (https://roman.gsfc.nasa.gov/science/Roman_Reference_Information.html via SVO)
   * - 160
     - CSST_nuv
     - Chinese Space Station Telescope nuv 
   * - 161
     - CSST_u
     - Chinese Space Station Telescope g 
   * - 162
     - CSST_g
     - Chinese Space Station Telescope i 
   * - 163
     - CSST_r
     - Chinese Space Station Telescope r 
   * - 164
     - CSST_i
     - Chinese Space Station Telescope u 
   * - 165
     - CSST_z
     - Chinese Space Station Telescope y 
   * - 166
     - CSST_y
     - Chinese Space Station Telescope z 
   * - 167
     - BASS_g
     - Beijing Arizona Sky Survey g-band used in DESI (https://speclite.readthedocs.io/en/latest/filters.html)
   * - 168
     - BASS_r
     - Beijing Arizona Sky Survey r-band used in DESI (https://speclite.readthedocs.io/en/latest/filters.html)
   * - 169
     - MzLS_z
     - Mayall z-band Legacy Survey z-band used in DESI (https://speclite.readthedocs.io/en/latest/filters.html)
   
  
  
