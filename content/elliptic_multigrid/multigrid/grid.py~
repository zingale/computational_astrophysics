import numpy as np

class Grid:
    def __init__(self, nx, ng=1, xmin=0, xmax=1,
                 bc_left_type="dirichlet", bc_left_val=0.0,
                 bc_right_type="dirichlet", bc_right_val=0.0):

        self.xmin = xmin
        self.xmax = xmax
        self.ng = ng
        self.nx = nx

        self.bc_left_type = bc_left_type
        self.bc_left_val = bc_left_val

        self.bc_right_type = bc_right_type
        self.bc_right_val = bc_right_val

        # python is zero-based.  Make easy intergers to know where the
        # real data lives
        self.ilo = ng
        self.ihi = ng+nx-1

        # physical coords -- cell-centered
        self.dx = (xmax - xmin)/(nx)
        self.x = xmin + (np.arange(nx+2*ng)-ng+0.5)*self.dx

        # storage for the solution
        self.v = self.scratch_array()
        self.f = self.scratch_array()
        self.r = self.scratch_array()

    def scratch_array(self):
        """return a scratch array dimensioned for our grid """
        return np.zeros((self.nx+2*self.ng), dtype=np.float64)

    def norm(self, e):
        """compute the L2 norm of e that lives on our grid"""
        return np.sqrt(self.dx * np.sum(e[self.ilo:self.ihi+1]**2))

    def residual_norm(self):
        """compute the residual norm"""
        r = self.scratch_array()
        r[self.ilo:self.ihi+1] = self.f[self.ilo:self.ihi+1] - (self.v[self.ilo+1:self.ihi+2] -
                                                                2 * self.v[self.ilo:self.ihi+1] +
                                                                self.v[self.ilo-1:self.ihi]) / self.dx**2
        return self.norm(r)

    def source_norm(self):
        """compute the source norm"""
        return self.norm(self.f)

    def fill_bcs(self):
        """fill the boundary conditions on phi"""

        # we only deal with a single ghost cell here

        # left
        if self.bc_left_type.lower() == "dirichlet":
            self.v[self.ilo-1] = 2 * self.bc_left_val - self.v[self.ilo]
        elif self.bc_left_type.lower() == "neumann":
            self.v[self.ilo-1] = self.v[self.ilo] - self.dx * self.bc_left_val
        else:
            raise ValueError("invalid bc_left_type")

        # right
        if self.bc_right_type.lower() == "dirichlet":
            self.v[self.ihi+1] = 2 * self.bc_right_val - self.v[self.ihi]
        elif self.bc_right_type.lower() == "neumann":
            self.v[self.ihi+1] = self.v[self.ihi] - self.dx * self.bc_right_val
        else:
            raise ValueError("invalid bc_right_type")

    def restrict(self, comp="v"):

        # create a coarse array
        ng = self.ng
        nc = self.nx//2

        ilo_c = ng
        ihi_c = ng + nc - 1

        coarse_data = np.zeros((nc + 2*ng), dtype=np.float64)

        if comp == "v":
            fine_data = self.v
        elif comp == "f":
            fine_data = self.f
        elif comp == "r":
            fine_data = self.r
        else:
            raise ValueError("invalid component")

        coarse_data[ilo_c:ihi_c+1] = 0.5 * (fine_data[self.ilo:self.ihi+1:2] +
                                            fine_data[self.ilo+1:self.ihi+1:2])

        return coarse_data

    def prolong(self, comp="v"):

        """
        prolong the data in the current (coarse) grid to a finer
        (factor of 2 finer) grid.  Return an array with the resulting
        data (and same number of ghostcells).

        We will reconstruct the data in the zone from the
        zone-averaged variables using the centered-difference slopes

                  (x)
        f(x,y) = m    x/dx + <f>

        When averaged over the parent cell, this reproduces <f>.

        Each zone's reconstrution will be averaged over 2 children.

        |           |     |     |     |
        |    <f>    | --> |     |     |
        |           |     |  1  |  2  |
        +-----------+     +-----+-----+

        We will fill each of the finer resolution zones by filling all
        the 1's together, using a stride 2 into the fine array.  Then
        the 2's, this allows us to operate in a vector
        fashion.  All operations will use the same slopes for their
        respective parents.

        """

        if comp == "v":
            coarse_data = self.v
        elif comp == "f":
            coarse_data = self.f
        elif comp == "r":
            coarse_data = self.r
        else:
            raise ValueError("invalid component")


        # allocate an array for the coarsely gridded data
        ng = self.ng
        nf = self.nx * 2

        fine_data = np.zeros((nf + 2*ng), dtype=np.float64)

        ilo_f = ng
        ihi_f = ng + nf - 1

        # slopes for the coarse data
        m_x = self.scratch_array()
        m_x[self.ilo:self.ihi+1] = 0.5 * (coarse_data[self.ilo+1:self.ihi+2] -
                                          coarse_data[self.ilo-1:self.ihi])

        # fill the '1' children
        fine_data[ilo_f:ihi_f+1:2] = \
            coarse_data[self.ilo:self.ihi+1] - 0.25 * m_x[self.ilo:self.ihi+1]

        # fill the '2' children
        fine_data[ilo_f+1:ihi_f+1:2] = \
            coarse_data[self.ilo:self.ihi+1] + 0.25 * m_x[self.ilo:self.ihi+1]

        return fine_data


class Multigrid:
    """
    The main multigrid class for cell-centered data.

    We require that nx be a power of 2 for simplicity
    """

    def __init__(self, nx, xmin=0.0, xmax=1.0,
                 bc_left_type="dirichlet", bc_right_type="dirichlet",
                 nsmooth=10, nsmooth_bottom=50,
                 verbose=0,
                 true_function=None):

        self.nx = nx
        self.ng = 1

        self.xmin = xmin
        self.xmax = xmax

        self.nsmooth = nsmooth
        self.nsmooth_bottom = nsmooth_bottom

        self.max_cycles = 100

        self.verbose = verbose

        self.bc_left_type = bc_left_type
        self.bc_right_type = bc_right_type

        # a function that gives the analytic solution (if available)
        # for diagnostics only
        self.true_function = true_function

        # a small number used in computing the error, so we don't divide by 0
        self.small = 1.e-16

        # assume that self.nx = 2^(nlevels-1)
        # this defines nlevels such that we end exactly on a 2 zone grid
        self.nlevels = int(np.log(self.nx)/np.log(2.0))

        # a multigrid object will be a list of grids
        self.grids = []

        # create the grids.  Here, self.grids[0] will be the coarsest
        # grid and self.grids[nlevel-1] will be the finest grid we
        # store the solution, v, the rhs, f.

        nx_t = 2
        for _ in range(self.nlevels):

            # add a grid for this level
            self.grids.append(Grid(nx_t, xmin=self.xmin, xmax=self.xmax,
                                   bc_left_type=self.bc_left_type,
                                   bc_right_type=self.bc_right_type))

            nx_t *= 2

        # provide coordinate and indexing information for the solution mesh
        self.soln_grid = self.grids[self.nlevels-1]

        self.ilo = self.soln_grid.ilo
        self.ihi = self.soln_grid.ihi

        self.x = self.soln_grid.x
        self.dx = self.soln_grid.dx

        # store the source norm
        self.source_norm = 0.0

        # after solving, keep track of the number of cycles taken, the
        # relative error from the previous cycle, and the residual error
        # (normalized to the source norm)

        self.num_cycles = 0
        self.residual_error = 1.e33
        self.relative_error = 1.e33

    def get_solution(self):
        return self.grids[self.nlevels-1].v.copy()

    def get_solution_object(self):
        return self.grids[self.nlevels-1]

    def init_solution(self):
        """
        initialize the solution to the elliptic problem as zero
        """
        self.soln_grid.v[:] = 0.0

    def init_rhs(self, data):
        self.soln_grid.f[:] = data.copy()

        # store the source norm
        self.source_norm = self.soln_grid.norm(self.soln_grid.f)

        if self.verbose:
            print("Source norm = ", self.source_norm)

        # note: if we wanted to do inhomogeneous Dirichlet BCs, we
        # would modify the source term, f, here to include a boundary
        # charge

    def _compute_residual(self, level):
        """ compute the residual and store it in the r variable"""

        myg = self.grids[level]

        # compute the residual
        # r = f - L phi
        myg.r[myg.ilo:myg.ihi+1] = myg.f[myg.ilo:myg.ihi+1] - \
            (myg.v[myg.ilo-1:myg.ihi] + myg.v[myg.ilo+1:myg.ihi+2] -
             2.0*myg.v[myg.ilo:myg.ihi+1]) / (myg.dx * myg.dx)

    def smooth(self, level, nsmooth):
        """ use Gauss-Seidel iterations to smooth """

        myg = self.grids[level]

        myg.fill_bcs()

        # do red-black G-S
        for _ in range(nsmooth):

            myg.v[myg.ilo:myg.ihi+1:2] = 0.5 * (
                -myg.dx * myg.dx * myg.f[myg.ilo:myg.ihi+1:2] +
                myg.v[myg.ilo+1:myg.ihi+2:2] + myg.v[myg.ilo-1:myg.ihi:2])

            myg.fill_bcs()

            myg.v[myg.ilo+1:myg.ihi+1:2] = 0.5 * (
                -myg.dx * myg.dx * myg.f[myg.ilo+1:myg.ihi+1:2] +
                myg.v[myg.ilo+2:myg.ihi+2:2] + myg.v[myg.ilo:myg.ihi:2])

            myg.fill_bcs()

    def solve(self, rtol=1.e-11):
        """do V-cycles util the L2 norm of the relative solution difference is
        < rtol

        """

        if self.verbose:
            print("source norm = ", self.source_norm)

        old_soln = self.soln_grid.v.copy()

        residual_error = 1.e33
        cycle = 1

        # diagnostics that are returned -- residual error norm and true
        # error norm (if possible) for each cycle
        rlist = []
        elist = []

        while residual_error > rtol and cycle <= self.max_cycles:

            # zero out the solution on all but the finest grid
            for level in range(self.nlevels-1):
                self.grids[level].v[:] = 0.0

            # descending part
            if self.verbose:
                print(f"<<< beginning V-cycle (cycle {cycle}) >>>\n")

            self.v_cycle(self.nlevels-1)

            # compute the error with respect to the previous solution
            # this is for diagnostic purposes only

            diff = (self.soln_grid.v - old_soln) / (self.soln_grid.v + self.small)
            relative_error = self.soln_grid.norm(diff)

            old_soln = self.soln_grid.v.copy()

            # compute the residual error, relative to the source norm
            self._compute_residual(self.nlevels-1)

            residual_error = self.soln_grid.norm(self.soln_grid.r)
            if self.source_norm != 0.0:
                residual_error /= self.source_norm

            if residual_error < rtol:
                self.num_cycles = cycle
                self.relative_error = relative_error
                self.residual_error = residual_error
                self.soln_grid.fill_bcs()

            if self.verbose:
                print(f"cycle {cycle}: relative err = {relative_error:11.6g}, residual err = {residual_error:11.6g}\n")

            rlist.append(residual_error)

            if self.true_function is not None:
                elist.append(self.soln_grid.norm(old_soln - self.true_function(self.soln_grid.x)))

            cycle += 1

        return elist, rlist

    def v_cycle(self, level):

        if level > 0:
            fp = self.grids[level]
            cp = self.grids[level-1]

            if self.verbose:
                self._compute_residual(level)
                old_res_norm = fp.norm(fp.r)

            # smooth on the current level
            self.smooth(level, self.nsmooth)

            # compute the residual
            self._compute_residual(level)

            if self.verbose:
                print(f"  level = {level}, nx = {fp.nx:4}, residual change: {old_res_norm:11.6g} -> {fp.norm(fp.r):11.6g}")

            # restrict the residual down to the RHS of the coarser level
            cp.f[:] = fp.restrict("r")

            # solve the coarse problem
            self.v_cycle(level-1)

            # prolong the error up from the coarse grid
            fp.v += cp.prolong("v")

            if self.verbose:
                self._compute_residual(level)
                old_res_norm = fp.norm(fp.r)

            # smooth
            self.smooth(level, self.nsmooth)

            if self.verbose:
                self._compute_residual(level)
                print(f"  level = {level}, nx = {fp.nx:4}, residual change: {old_res_norm:11.6g} -> {fp.norm(fp.r):11.6g}")

        else:
            # solve the discrete coarse problem just via smoothing
            if self.verbose:
                print("  bottom solve")

            bp = self.grids[0]

            self.smooth(0, self.nsmooth_bottom)

            bp.fill_bcs()
