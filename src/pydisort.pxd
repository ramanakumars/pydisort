cdef extern from "cdisort.h":
    # DISORT flags
    ctypedef struct disort_flag:
        int usrtau, usrang, ibcnd, lamber, planck, spher, onlyfl, brdf_type, quiet, intensity_correction, old_intensity_correction
        int prnt[5]

    # DISORT radiant information (output from running the model)
    ctypedef struct disort_radiant:
        double rfldir, rfldn, flup, dfdt, uavg, uavgdn, uavgup, uavgso
    
    # DISORT boundary conditions
    ctypedef struct disort_bc:
        double fbeam,  umu0, phi0, fisot, ttemp, btemp, temis, albedo

    # BRDF specifications for different models
    ctypedef struct rpv_brdf_spec:
        double rho0, k, theta, sigma, t1, t2, scale

    ctypedef struct ambrals_brdf_spec:
        double iso, vol, geo

    ctypedef struct cam_brdf_spec:
        double u10, pcl, xsal

    ctypedef struct disort_brdf:
      # brdf types */
      rpv_brdf_spec *rpv      # specification for rpv BRDF      */
      ambrals_brdf_spec *ambrals  # specification for ambrals BRDF  */
      cam_brdf_spec *cam       # specification for Cox&Munk BRDF */

    # Input state to the DISORT calculation
    ctypedef struct disort_state:
        char header[128]
        disort_flag flag
        disort_bc bc
        disort_brdf brdf
        int nlyr, nmom, nstr, nmom_nstr, ntau, numu, nphi, nphase         
        double wvnmlo, wvnmhi, accur, radius
        double *dtauc
        double *ssalb
        double *pmom
        double *temper
        double *utau
        double *umu
        double *phi
        double *zd
        double *mu_phase
        double *phase
     
    # Output from the DISORT model
    ctypedef struct disort_output:
        disort_radiant *rad       # See typedef disort_radiant
        double *albmed
        double *trnmed
        double *uu
        double *u0u

    cdef int SPECIAL_BC

    cdef int ISOTROPIC
    cdef int RAYLEIGH
    cdef int HENYEY_GREENSTEIN
    cdef int HAZE_GARCIA_SIEWERT
    cdef int CLOUD_GARCIA_SIEWERT

    extern void c_disort(disort_state *, disort_output *)
    extern void c_disort_out_alloc(disort_state *, disort_output *)
    extern void c_disort_state_alloc(disort_state *)
    extern void c_disort_out_free(disort_state *, disort_output *)
    extern void c_getmom(int, double, int, double *);
