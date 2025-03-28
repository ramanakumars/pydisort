from pydisort cimport disort_state, disort_output, c_disort_state_alloc, c_disort_out_alloc, c_disort_out_free, SPECIAL_BC
cimport cython
from libc.stdlib cimport malloc, free
import numpy as np

cdef copy_array_to_c(array_src, double* array_dest):
    """
    Copy an array from a Python list to a C array pointer

    :param array_src: source array
    :param array_dest: pointer to the destination array
    """
    for i in xrange(len(array_src)):
        array_dest[i] = array_src[i]

cdef copy_array_from_c(double *array_src, int size):
    """
    Copy an array from a C array pointer to a Python list

    :param array_src: pointer to the array
    :param size: size of the array

    :returns: Python list containing the array elements
    """
    array_dest = []
    for i in xrange(size):
        array_dest.append(array_src[i])

    return array_dest

cdef class PyDisort:
    cdef disort_state state
    cdef disort_output output

    def __init__(self):
        pass

    def set_flags(self, usrtau: cython.int, usrang: cython.int, ibcnd: cython.int, lamber: cython.int,
                  onlyfl: cython.int, quiet: cython.int, spher: cython.int,
                  old_intensity_correction: cython.int, planck: cython.int) -> cython.int:
        """
        Set the flags for DISORT as follows:

        :param usrtau: TRUE=> radiant quantities returned at user-specified optical depths  
        :param usrang: TRUE=> radiant quantities returned at user-specified polar angles    
        :param ibcnd: SPECIAL_BC => return only albedo and transmis., see Ref S2           
        :param lamber: TRUE=> isotropically reflecting bottom boundary, FALSE=>bi-dir.      
        :param planck: TRUE=>incl. thermal emission                                         
        :param spher: TRUE=>pseudo-spherical geometry, otherwise plane-parallel            
        :param onlyfl: FALSE=>return intensities in addition to other radiant quantities    
        :param quiet: quiet output                                                         
        :param intensity_correction: apply intensity correction                           
        :param old_intensity_correction: use original intensity correction routine            
        """
        self.state.flag.usrtau = usrtau
        self.state.flag.usrang = usrang
        self.state.flag.ibcnd = ibcnd
        self.state.flag.lamber = lamber
        self.state.flag.onlyfl = onlyfl
        self.state.flag.quiet = quiet
        self.state.flag.spher = spher
        self.state.flag.old_intensity_correction = old_intensity_correction
        self.state.flag.planck = planck

        return 0

    def set_boundary(self, umu0: cython.double, phi0: cython.double, fbeam: cython.double, fisot: cython.double, ttemp: cython.double,
                     btemp: cython.double, temis: cython.double, albedo: cython.double):
        """
        Set the boundary conditions for the atmosphere

        :param umu0: Polar angle cosine of incident beam (positive)                      
        :param phi0: Azimuth angle of incident beam (0 to 360 deg)                       
        :param fbeam: Intensity of incident parallel beam at top boundary                 
        :param fisot: Intensity of top-boundary isotropic illumination                    
        :param ttemp: Temperature [K] of top boundary                                     
        :param btemp: Temperature [K] of bottom boundary                                  
        :param temis: Emissivity of top boundary. Needed if planck = TRUE                 
        :param albedo: Albedo of bottom boundary, needed if lamber = TRUE                  
        """
        self.state.bc.fbeam = fbeam
        self.state.bc.fisot = fisot
        self.state.bc.albedo = albedo
        self.state.bc.umu0 = umu0
        self.state.bc.phi0 = phi0

        return 0

    def set_input(self, tau: iterable, ssa: iterable, umu: iterable, phi: iterable,
                  pmom: iterable, nstr: cython.int, nphase: cython.int, nmom: cython.int):
        """
        Set the inputs for the model, including the sizes of the grids

        :param tau: input optical depth profile (shape: nlayers + 1)
        :param ssa: input single-scattering albedo (shape: nlayers + 1)
        :param umu: input cosine of sampling incidence angles (shape: numu)
        :param uphi: input cosine of sampling phase angles (shape: nphi)
        :param pmom: phase functions for each of the moments (shape: (max(nmom, nstr) + 1) * nlayers)
        :param nstr: number of streams
        :param nphase: number of phase angles
        :param nmom: number of moments

        """
        self.state.accur = 0.

        # set the sizes of the inputs
        self.state.nstr = nstr
        self.state.nphase = nphase
        self.state.nlyr = len(tau) - 1
        self.state.nmom = nmom
        self.state.numu = len(umu)
        self.state.nphi = len(phi)

        # allocate the memory for the state arrays
        c_disort_state_alloc(&self.state)

        # copy the data over
        copy_array_to_c(tau, self.state.dtauc)
        copy_array_to_c(ssa, self.state.ssalb)
        copy_array_to_c(umu, self.state.umu)
        copy_array_to_c(phi, self.state.phi)
        copy_array_to_c(pmom, self.state.pmom)

        return 0

    def run(self):
        """
        Run the DISORT model
        :returns: dictionary containing the radiant, fluxes, albedo and transmittivity (latter two in the case of SPECIAL_BC)
        """
        c_disort_out_alloc(&self.state, &self.output)
        c_disort(&self.state, &self.output)

        radiant = {
            "rfldir": self.output.rad.rfldir, "rfldn": self.output.rad.rfldn, "flup": self.output.rad.flup,
            "dfdt": self.output.rad.dfdt, "uavg": self.output.rad.uavg, "uavgdn": self.output.rad.uavgdn,
            "uavgup": self.output.rad.uavgup, "uavgso": self.output.rad.uavgso
        }

        uu = np.asarray(copy_array_from_c(self.output.uu, self.state.numu * self.state.ntau * self.state.nphi)).reshape((self.state.ntau, self.state.numu))
        u0u = np.asarray(copy_array_from_c(self.output.u0u, self.state.numu * self.state.ntau * self.state.nphi)).reshape((self.state.ntau, self.state.numu))
        
        if self.state.flag.ibcnd == SPECIAL_BC:
            albmed = np.asarray(copy_array_from_c(self.output.albmed, self.state.numu))
            trnmed = np.asarray(copy_array_from_c(self.output.trnmed, self.state.numu))
        else:
            albmed = None
            trnmed = None

        c_disort_out_free(&self.state, &self.output)

        return {"radiant": radiant,
                "uu": uu, 
                "u0u": u0u,
                "albmed": albmed,
                "trnmed": trnmed}
