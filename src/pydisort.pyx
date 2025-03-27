from pydisort cimport disort_state, disort_output
cimport cython
from libc.stdlib cimport malloc, free

cdef copy_array_to_c(array_src, double* array_dest):
    for i in xrange(len(array_src)):
        array_dest[i] = array_src[i]

cdef copy_array_from_c(double *array_src, int size):
    array_dest = []
    for i in xrange(size):
        array_dest.append(array_src[i])

    return array_dest

cdef class PyDisort:
    cdef disort_state state
    cdef disort_output output

    def __init__(self):
        pass


    def set_flags(self, usrtau, usrang, ibcnd, lamber, onlyfl, quiet, spher, old_intensity_correction, planck):
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

    def set_boundary(self, fbeam, fisot, albedo, umu0, phi0):
        self.state.bc.fbeam = fbeam
        self.state.bc.fisot = fisot
        self.state.bc.albedo = albedo
        self.state.bc.umu0 = umu0
        self.state.bc.phi0 = phi0

        return 0

    def set_input(self, tau, ssa, umu, phi, nstr, nphase, pmom, nmom):
        self.state.accur = 0.

        self.state.nstr = nstr
        self.state.nphase = nphase
        self.state.nlyr = len(tau) - 1
        self.state.nmom = nmom
        self.state.numu = len(umu)
        self.state.nphi = len(phi)

        c_disort_state_alloc(&self.state)
        copy_array_to_c(tau, self.state.dtauc)
        copy_array_to_c(ssa, self.state.ssalb)
        copy_array_to_c(umu, self.state.umu)
        copy_array_to_c(phi, self.state.phi)
        copy_array_to_c(pmom, self.state.pmom)

        return 0

    def run(self):
        c_disort_out_alloc(&self.state, &self.output)
        c_disort(&self.state, &self.output)

        uu = np.asarray(copy_array_from_c(self.output.uu, self.state.numu * self.state.ntau * self.state.nphi)).reshape((self.state.nlyr, self.state.ntau, self.state.nphi))
        u0u = np.asarrary(copy_array_from_c(self.output.u0u, self.state.numu * self.state.ntau * self.state.nphi).reshape((self.state.nlyr, self.state.ntau, self.state.nphi))

        return self.output.rad.flup, uu, u0u
