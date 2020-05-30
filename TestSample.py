# --------------------------------------------------------------------------------- #
# CASE DESCRIPTION
#
#               _____________________             ___
#              |                     |             |
#              |                     |             |
#              |                     |             |
#           >>||                     |             |
#        Li >>||                     |             | Ly
#           >>||                     |             |
#              |                     ||>>          |
#              |                     ||>> Lo       |
#              |_____________________||>>         _|_
#
#              |---------------------|
#                         Lx
#
# ASSUMPTIONS
# 1. Steady state
# 2. Slip wall
# 3. Incompressible


# --------------------------------------------------------------------------------- #
import numpy as np

import cores.CoreMethods as CoreMethod
import cores.PlotMethods as PlotMethod


# --------------------------------------------------------------------------------- #
# INPUT
# Domain dimension
Lx = 0.1                        # Length in x-direction
Ly = 0.1                        # Length in y-direction
Li = 0.02                       # Inlet dimension
Lo = 0.02                       # Outlet dimension
Loc = 0.04                      # Inlet position
Nx = 21                         # number of nodes in x-direction
Ny = 21                         # number of nodes in y-direction

# Fluid properties
mu = 0.01                       # Dynamic viscosity, Pa.s
rho = 1                         # Density, kg/m3
pref = 101325                   # Reference pressure, Pa

# Boundary conditions
uInlet = 1.0                    # Inlet velocity, m/s

# Solver properties
alpha_u = 0.7                   # Under-relaxation factor of u-momentum
alpha_v = 0.7                   # Under-relaxation factor of v-momentum
alpha_p = 0.3                   # Under-relaxation factor of pressure
error_continuity = 1e-3         # Target error of continuity equation
error_momentum_u = 1e-3         # Target error of u-momentum equation
error_momentum_v = 1e-3         # Target error of v-momentum equation


# --------------------------------------------------------------------------------- #
# PROCESS
# Geometry and grid generation
xc, yc = CoreMethod.Cartesian(Lx, Ly, Nx, Ny)
xs, ys = CoreMethod.Scalar(xc, yc, Nx, Ny)
xu, yu = CoreMethod.VectorU(xc, yc, Nx, Ny)
xv, yv = CoreMethod.VectorV(xc, yc, Nx, Ny)


# Initialize flow variables
ps, pp, us, vs, uOutlet, dx, dy, nInlet, nOutlet, nLoc = CoreMethod.InitializeVariables(uInlet, Lx, Ly, Li, 
                                                                                        Lo, Loc, Nx, Ny)


# Applying boundary condition
us, vs, ps = CoreMethod.BoundaryCondition(us, vs, ps, uInlet, uOutlet, nInlet, nOutlet, nLoc)


# Start iteration
resCont = 0.01
resMomu = 0.01
resMomv = 0.01
ITER = 0
ITERMAX = 150

iterCount = np.zeros(ITERMAX)
rmseCont = np.ones(ITERMAX)
rmseMomu = np.ones(ITERMAX)
rmseMomv = np.ones(ITERMAX)

print("")

#while (resCont>error_continuity) or (resMomu>error_momentum_u) or (resMomv>error_momentum_v):
while (ITER<ITERMAX):
    # SOLVE WITH SIMPLE ALGORITHM
    us, vs, ps, resCont, resMomu, resMomv = CoreMethod.SIMPLEAlgorithm(us, vs, ps, pp, rho, mu, 
                                                                        alpha_u, alpha_v, alpha_p, 
                                                                        uInlet, uOutlet, nInlet, nOutlet,
                                                                        nLoc, Lx, Ly, xu, yu, xv, yv, 
                                                                        xc, yc, xs, ys, Nx, Ny)

    # SAVE RESIUDAL
    
    rmseCont[ITER] = resCont
    rmseMomu[ITER] = resMomu
    rmseMomv[ITER] = resMomv
    iterCount[ITER] = ITER + 1

    # UPDATE POINTER
    ITER = ITER + 1
    CoreMethod.updateProgress((ITER)/ITERMAX, resCont, resMomu, resMomv)

    #print(f"Error-{ITER} = {np.round(resCont,3)} \t {np.round(resMomu,3)} \t {np.round(resMomv,3)}")


# --------------------------------------------------------------------------------- #
# RESULTS VISUALIZATION
# Residuals
PlotMethod.Residue(iterCount, rmseCont, rmseMomu, rmseMomv)

# Grid
#PlotMethod.Grid(xc, yc, xs, ys, xu, yu, xv, yv)

# CFD Results
u_nodes, v_nodes, p_nodes, v_res = CoreMethod.NodesVal(us, vs, ps, Nx, Ny)
PlotMethod.Velocity(xs[1:-1,1:-1], ys[1:-1,1:-1], v_res, u_nodes, v_nodes, "resultant", 100)
PlotMethod.Velocity(xs[1:-1,1:-1], ys[1:-1,1:-1], v_res, u_nodes, v_nodes, "u", 100)
PlotMethod.Velocity(xs[1:-1,1:-1], ys[1:-1,1:-1], v_res, u_nodes, v_nodes, "v", 100)
PlotMethod.Pressure(xs[1:-1,1:-1], ys[1:-1,1:-1], p_nodes, 100)
PlotMethod.Streamline(xs[1:-1,1:-1], ys[1:-1,1:-1], v_res, u_nodes, v_nodes, 100)
PlotMethod.Compare(xs[1:-1,1:-1], ys[1:-1,1:-1], v_res, u_nodes, v_nodes, p_nodes, Nx, Ny)
# --------------------------------------------------------------------------------- #