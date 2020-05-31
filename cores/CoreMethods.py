import math
import numpy as np
import sys
import time


# --------------------------------------------------------------------------------- #
# GRID GENERATION
def Cartesian(Lx, Ly, Nx, Ny):
    # Generate Cartesian/Physical Grids
    xc = np.zeros((Ny, Nx))
    yc = np.zeros((Ny, Nx))

    yc[0,:] = Ly
    for i in range(1, Nx):
        xc[:,i] = xc[:,i-1] + Lx/(Nx-1)
    
    for j in range(1, Ny):
        yc[j,:] = yc[j-1,:] - Ly/(Ny-1)

    return (xc, yc)

def Scalar(Xc, Yc, Nx, Ny):
    # Generate Staggered Grids (Scalar)
    xs = np.zeros((Ny+1, Nx+1))
    ys = np.zeros((Ny+1, Nx+1))

    # Inner nodes
    for j in range(1, Ny):
        for i in range(1, Nx):
            xs[j,i] = (Xc[j-1,i-1] + Xc[j-1,i] + Xc[j,i-1] + Xc[j,i])/4
            ys[j,i] = (Yc[j-1,i-1] + Yc[j-1,i] + Yc[j,i-1] + Yc[j,i])/4

    # Additions
    # West
    xs[:,0] = xs[:,1] - abs(xs[:,2] - xs[:,1])
    ys[:,0] = ys[:,1] - abs(ys[:,2] - ys[:,1])
    
    # East
    xs[:,-1] = xs[:,-2] + abs(xs[:,-2] - xs[:,-3])
    ys[:,-1] = ys[:,-2] + abs(ys[:,-2] - ys[:,-3])

    # North
    xs[0,:]  = xs[1,:] + abs(xs[2,:] - xs[1,:])
    ys[0,:]  = ys[1,:] + abs(ys[2,:] - ys[1,:])

    # South
    xs[-1,:] = xs[-2,:] - abs(xs[-2,:] - xs[-3,:])
    ys[-1,:] = ys[-2,:] - abs(ys[-2,:] - ys[-3,:])

    #CORNER POINTS
    xs[0,0]= xs[1,0]
    ys[0,0]= ys[0,1]

    xs[0,-1] = xs[1,-1]
    ys[0,-1] = ys[0,-2]

    xs[-1,0] = xs[-2,0]
    ys[-1,0] = ys[-1,1]

    xs[-1,-1] = xs[-2,-1]
    ys[-1,-1] = ys[-1,-2]

    return (xs, ys)

def VectorU(Xc, Yc, Nx, Ny):
    xu = np.zeros((Ny+1, Nx))
    yu = np.zeros((Ny+1, Nx))

    # Inner nodes
    for j in range(1, Ny):
        for i in range(Nx):
            xu[j,i] = 0.5*(Xc[j-1,i] + Xc[j,i])
            yu[j,i] = 0.5*(Yc[j-1,i] + Yc[j,i])

    # Outer nodes
    xu[0,:] = xu[1,:]
    yu[0,:] = yu[1,:] + 2*abs(Yc[0,:] - yu[1,:])

    xu[-1,:] = xu[-2,:]
    yu[-1,:] = yu[-2,:] - 2*abs(Yc[-1,:] - yu[-2,:])

    return (xu, yu)

def VectorV(Xc, Yc, Nx, Ny):
    xv = np.zeros((Ny, Nx+1))
    yv = np.zeros((Ny, Nx+1))

    # Inner nodes
    for j in range(Ny):
        for i in range(1, Nx):
            xv[j,i] = 0.5*(Xc[j,i-1] + Xc[j,i])
            yv[j,i] = 0.5*(Yc[j,i-1] + Yc[j,i])

    # Outer nodes
    xv[:,0] = xv[:,1] - 2*abs(Xc[:,0]-xv[:,1])
    yv[:,0] = yv[:,1] 

    xv[:,-1] = xv[:,-2] +  2*abs(Xc[:,-1] - xv[:,-2])
    yv[:,-1] = yv[:,-2] 

    return (xv, yv)


# PRE-PROCESSING
def InitializeVariables(uInlet, Lx, Ly, Linlet, Loutlet, Loc, Nx, Ny):
    u = np.zeros((Ny+1, Nx)).astype(float)
    v = np.zeros((Ny, Nx+1)).astype(float)
    p = np.zeros((Ny+1, Nx+1)).astype(float)
    pp = np.zeros((Ny+1, Nx+1)).astype(float)

    u[:,:] = 0.0#1e-20
    v[:,:] = 0.0#1e-20

    dx = Lx/(Nx-1)
    dy = Ly/(Ny-1)
    nInlet = math.ceil(Linlet/dy)
    nOutlet = math.ceil(Loutlet/dy)
    nLoc = math.ceil(Loc/dy)

    uOutlet = uInlet*Linlet/Loutlet         # Outlet velocity based on continuity

    return (p, pp, u, v, uOutlet, dx, dy, nInlet, nOutlet, nLoc)

def BoundaryCondition(u, v, p, uInlet, uOutlet, nInlet, nOutlet, nLoc):
    # Boundari Condition
    # Inlet     : Velocity
    # Outlet    : Velocity
    # Wall      : No-slip
    
    # u-velocity
    # Left domain: Stationary Wall + Inlet
    u[1:(nLoc+1), 0] = 0                        # for wall
    u[(nLoc+1):(nLoc+1+nInlet), 0] = uInlet     # for inlet

    # Right domain: Stationary Wall + Outlet
    u[-(1+nOutlet):-1, -1] = uOutlet              # for outlet
    u[1:-(1+nOutlet),-1] = 0                      # for wall

    # Top & bottom
    u[0,:] = -u[1,:]
    u[-1,:] = -u[-2,:]

    # v-velocity
    # Top & Bottom: Stationary Wall
    v[0,:] = 0
    v[-1,:] = 0

    # Right & Left
    v[:,0] = -v[:,1]
    v[:,-1] = -v[:,-2]

    # Pressure    
    # At wall (dp/dy = 0)
    # Top and Bottom
    #p[0,:] = -p[1,:]
    #p[-1,:] = -p[-2,:]

    # Left and Right
    #p[:,0] = -p[:,1]
    #p[:,-1] = -p[:,-2]
    #p[:,-2] = p[:,-3]
    #p[:,-1] = -p[:,-2]

    #pMin = np.min(p)

    # Outlet
    #p[-(1+nOutlet):-1, -2] = 0.0
    #p[-(1+nOutlet):-1, -1] = 0.0

    #p[-(1+nOutlet):-1, -2] = p[-(1+nOutlet):-1, -2] - pMin
    #p[-(1+nOutlet):-1, -1] = p[-(1+nOutlet):-1, -1] - pMin

    return (u, v, p)


# SOLVER
def CorrectContinuity(p, pp, alpha_p, Nx, Ny):
    for i in range(1,Nx):
        for j in range(1,Ny):
            p[j,i] = p[j,i] + alpha_p*pp[j,i]

    return (p)

def CorrectMomentum(u, v, pp, Apu, Apv, alpha_u, alpha_v, alpha_p, xu, yu, xv, yv, Nx, Ny):
    S = 1e-8
    
    # Correct x
    for i in range(1, Nx-1):
        for j in range(1, Ny):
            # Calculate Cell Faces Area
            A = 0.5*(abs(yv[j-1,i] - yv[j,i]) + abs(yv[j-1,i+1] - yv[j,i+1]))

            # Calculate Coefficient d
            d = A/(Apu[j,i] + S)

            # Calculate Corrected velocity
            u[j,i] = u[j,i] + d*(pp[j,i] - pp[j,i+1])

    # Correct y
    for i in range(1, Nx):
        for j in range(1, Ny-1):
            # Calculate Cell Faces Area
            A = 0.5*(abs(xu[j,i] - xu[j,i-1]) + abs(xu[j+1,i] - xu[j+1,i-1]))

            # Calculate Coefficient d
            d = A/(Apv[j,i] + S)

            # Calculate Corrected velocity
            v[j,i] = v[j,i] + d*(pp[j+1,i] - pp[j,i])

    
    return (u, v)

def ErrorContinuity(Nx, Ny, R):
    error_mat = np.zeros((Ny-1, Nx-1))
    for j in range(Ny-1):
        for i in range(Nx-1):
            error_mat[j,i] = np.sqrt((R[j+1, i+1])**2)
    
    sum_error_mat = sum(error_mat)
    error = sum(sum_error_mat)

    return (error)

def ErrorMomentum(Nx, Ny, Ru, Rv):
    # u-momentum
    error_mat_u = np.zeros((Ny+1, Nx))
    for j in range(Ny+1):
        for i in range(Nx):
            error_mat_u[j,i] = np.sqrt((Ru[j,i])**2)

    # v-momentum
    error_mat_v = np.zeros((Ny, Nx+1))
    for j in range(Ny):
        for i in range(Nx+1):
            error_mat_v[j,i] = np.sqrt((Rv[j,i])**2)

    sum_error_mat_u = sum(error_mat_u)
    sum_error_mat_v = sum(error_mat_v)

    error_u = sum(sum_error_mat_u)
    error_v = sum(sum_error_mat_v)

    return (error_u, error_v)

def NodesVal(u, v, p, Nx, Ny):
    uNodes = np.zeros((Ny-1, Nx-1))
    vNodes = np.zeros((Ny-1, Nx-1))
    pNodes = np.zeros((Ny-1, Nx-1))
    vRes = np.zeros((Ny-1, Nx-1))

    for i in range(Nx-1):
        uNodes[:,i] = 0.5*(u[1:-1,i] + u[1:-1,i+1])

    for j in range(Ny-1):
        vNodes[j,:] = 0.5*(v[j,1:-1] + v[j+1,1:-1])

    for j in range(Ny-1):
        for i in range(Nx-1):
            vRes[j,i] = np.sqrt(uNodes[j,i]**2 + vNodes[j,i]**2)

    #p[:,-2] = p[:,-3]
    #p[-2,:] = p[-3,:]
    pNodes = p[1:-1, 1:-1]

    pNodes = pNodes - np.min(pNodes)
    #pNodes = pNodes*1.5816

    return (uNodes, vNodes, pNodes, vRes)

def Solveu(u, v, p, rho, mu, alpha_u, xc, yc, xu, yu, xv, yv, Nx, Ny):
    Ap = np.zeros((Ny+1, Nx))
    A = np.zeros((Nx-2, Nx-2))
    B = np.zeros((Nx-1, 1))
    rmu = np.zeros((Ny+1, Nx))

    #u momentum
    for j in range (1,Ny):
        for i in range (1,Nx-1):
            #face area
            AreaW = abs(yv[j-1,i] - yv[j,i])
            AreaE = abs(yv[j-1,i+1] - yv[j,i+1])
            AreaN = abs(xv[j-1,i+1] - xv[j-1,i])
            AreaS = abs(xv[j,i+1] - xv[j,i])
            
            dxu = 0.5*(AreaN + AreaS)
            dyu = 0.5*(AreaW + AreaE)
            
            Vol = dxu*dyu
            
            dWP = abs(xu[j,i] - xu[j,i-1])
            dPE = abs(xu[j,i+1] - xu[j,i])
            dNP = abs(yu[j-1,i] - yu[j,i])
            dPS = abs(yu[j,i] - yu[j+1,i])
            
            #diffusi coef
            G = mu
            Dw = G*AreaW/dWP
            De = G*AreaE/dPE
            Ds = G*AreaS/dPS
            Dn = G*AreaN/dNP
            
            #convective coef
            uw  = 0.5*(u[j,i] 	+ u[j,i-1])
            ue  = 0.5*(u[j,i]   + u[j,i+1])
            vn  = 0.5*(v[j-1,i] + v[j-1,i+1])
            vs  = 0.5*(v[j,i]   + v[j,i+1])
            
            Fe = rho*ue*AreaE
            Fw = rho*uw*AreaW
            Fn = rho*vn*AreaN
            Fs = rho*vs*AreaS
            
            deltaF  = Fe - Fw + Fn - Fs
            
            #Hybrid scheme
            Aw = max(Fw, (Dw+Fw/2), 0)
            Ae = max(-Fe, (De-Fe/2), 0)
            An = max(-Fn, (Dn-Fn/2), 0)
            As = max(Fs, (Ds+Fs/2), 0)
            
            Ap[j,i] = Ae + Aw + As + An + deltaF
            
            #CALCULATE SOURCE TERMS
            #Source due to pressure gradient
            SuPreGrad = (p[j,i]-p[j,i+1])*Vol/dxu

            #Source due to relaxation factor
            relaxation = ((1-alpha_u)*Ap[j,i]/alpha_u)*u[j,i]
            
            #UNDER RELAXATION CALCULATION
            Ap[j,i] = Ap[j,i]/alpha_u
            
            #ASSEMBLY SOURCE(B) & A MATRIX
            A[i-1,i-1]  = Ap[j,i]
            
            if i==1:
                A[i-1,i] = -Ae
                B[i-1,0] = SuPreGrad + relaxation + An*u[j-1,i] + As*u[j+1,i] + Aw*u[j,i-1]
                
            elif i==Nx-2:
                A[i-1,i-2] = -Aw
                B[i-1,0] = SuPreGrad + relaxation + An*u[j-1,i] + As*u[j+1,i] + Ae*u[j,i+1]
            else:
                A[i-1,i-2] = -Aw
                A[i-1,i]   = -Ae
                B[i-1,0] = SuPreGrad + relaxation + An*u[j-1,i] + As*u[j+1,i]

            rmu[j,i] = Ap[j,i]*u[j,i] - (Aw*u[j,i-1] + Ae*u[j,i+1] + An*u[j-1,i] + As*u[j+1,i]) - (SuPreGrad + relaxation)
    
        #SOLVE DISCRETIZED U MOMENTUM EQUATION
        u[j,1:-1]  = (TDMA(A,B)).transpose()

    return (u, Ap, rmu)

def Solvev(u, v, p, rho, mu, alpha_v, xc, yc, xu, yu, xv, yv, Nx, Ny):
    Ap      = np.zeros((Ny,Nx+1))
    A		= np.zeros((Ny-2,Nx-2))
    B		= np.zeros((Ny-2,1))
    rmv     = np.zeros((Ny,Nx+1))   
    
    #v-momentum
    for i in range (1,Nx):
        for j in range (1,Ny-1):
            #CALCULATE CELL FACES AREA
            AreaW = abs(yu[j,i-1] - yu[j+1,i-1])
            AreaE = abs(yu[j,i]   - yu[j+1,i]) 
            AreaN = abs(xu[j,i]   - xu[j,i-1])
            AreaS = abs(xu[j+1,i] - xu[j+1,i-1])
            
            dWP = abs(xv[j,i]   - xv[j,i-1])
            dPE = abs(xv[j,i+1] - xv[j,i])
            dNP = abs(yv[j-1,i] - yv[j,i])
            dPS = abs(yv[j,i]   - yv[j+1,i])
            
            dxv = 0.5*(AreaN+AreaS)	   
            dyv = 0.5*(AreaW+AreaE)
            
            vol	= dxv*dyv
            
            #CALCULATE DIFFUSIVE COEFFICIENTS
            G       = mu
            Dw    	= G *AreaW/dWP
            De    	= G *AreaE/dPE
            Dn    	= G *AreaN/dNP
            Ds    	= G *AreaS/dPS
            
            #CALCULATE CONVECTIVE COEFFICIENTS
            uw = 0.5*(u[j,i-1] + u[j+1,i-1])
            ue = 0.5*(u[j,i]   + u[j+1,i])
            vn = 0.5*(v[j-1,i] + v[j,i])
            vs = 0.5*(v[j,i]   + v[j+1,i])
            
            Fw  = rho*uw*AreaW
            Fe  = rho*ue*AreaE
            Fs  = rho*vs*AreaS
            Fn  = rho*vn*AreaN
            
            deltaF  = Fe - Fw + Fn - Fs

            # Hybrid scheme
            Aw = max(Fw,(Dw+Fw/2),0)
            Ae = max(-Fe,(De-Fe/2),0)
            An = max(-Fn,(Dn-Fn/2),0)
            As = max(Fs,(Ds+Fs/2),0)
            
            Ap[j,i] = Aw + Ae + As + An + deltaF

            
            #CALCULATE SOURCE TERMS
            #Source due to pressure gradient
            SvPreGrad = (p[j+1,i] - p[j,i])*vol/dyv
                    
            #Source due to relaxation factor
            relaxation  = ((1-alpha_v)*Ap[j,i]/alpha_v)*v[j,i]
            
            #UNDER RELAXATION CALCULATION
            Ap[j,i] = Ap[j,i]/alpha_v
            
            #ASSEMBLY SOURCE(B) & A MATRIX
            A[j-1,j-1] = Ap[j,i]
            
            if j==1:
                A[j-1,j] = -As
                B[j-1,0] = SvPreGrad + relaxation + Aw*v[j,i-1] + Ae*v[j,i+1] + An*v[j-1,i]
            
            elif j==Ny-2:
                A[j-1,j-2] = -An
                B[j-1,0]   = SvPreGrad + relaxation + Aw*v[j,i-1] + Ae*v[j,i+1] + As*v[j+1,i]
            
            else : 
                A[j-1,j-2] = -An
                A[j-1,j]   = -As
                B[j-1,0]   = SvPreGrad + relaxation + Aw*v[j,i-1] + Ae*v[j,i+1]
            
            #CALCULATE ERROR
            rmv[j,i] = Ap[j,i]*v[j,i] - (Aw*v[j,i-1] + Ae*v[j,i+1] + An*v[j-1,i] + As*v[j+1,i]) - (SvPreGrad + relaxation)
            
        #SOLVE DISCRETIZED V MOMENTUM EQUATION
        v[1:-1,i]  = (TDMA(A,B)).transpose()
    
    return (v, Ap, rmv)

def Solvep(u, v, Apu, Apv, rho, alpha_u, alpha_v, alpha_p, xc, yc, Nx, Ny):
    # Initialization
    pp = np.zeros((Ny+1, Nx+1))
    b = np.zeros((Ny+1, Nx+1))
    s = np.zeros((Ny+1, Nx+1))
    Aw = np.zeros((Ny+1, Nx+1))
    Ae = np.zeros((Ny+1, Nx+1))
    As = np.zeros((Ny+1, Nx+1))
    An = np.zeros((Ny+1, Nx+1))
    Ap = np.zeros((Ny+1, Nx+1))
    S = 1e-8

    # Start sweeping the whole nodes
    for i in range(1, Nx):
        for j in range(1, Ny):
            # GEOMETRY PREPARATION
            # Face area
            AreaW = abs(yc[j,i-1] - yc[j-1,i-1])
            AreaE = abs(yc[j,i] - yc[j-1,i])
            AreaN = abs(xc[j-1,i] - xc[j-1,i-1])
            AreaS = abs(xc[j,i-1] - xc[j,i])

            # Average cell length
            deltaXS = 0.5*(AreaN + AreaS)
            deltaYS = 0.5*(AreaW + AreaE)

            vol = deltaXS*deltaYS                   # Cell volumes


            # CALCULATE d (AP/A)
            # Coefficient manipulation for boundary
            dw = AreaW*alpha_u/(Apu[j,i-1] + S)
            de = AreaE*alpha_u/(Apu[j,i] + S)
            dn = AreaN*alpha_v/(Apv[j-1,i] + S)
            ds = AreaS*alpha_v/(Apv[j,i] + S)

            # Assembly coefficient
            Aw[j,i] = rho*dw*AreaW
            Ae[j,i] = rho*de*AreaE
            An[j,i] = rho*dn*AreaN
            As[j,i] = rho*ds*AreaS

            # Boundary correction
            Aw[:,1]  = 0
            Ae[:,-2] = 0 
            An[1,:]  = 0
            As[-2,:] = 0

            Ap[j,i] = Aw[j,i] + Ae[j,i] + An[j,i] + As[j,i]


            # CALCULATE CONTINUITY
            Fw = rho*u[j,i-1]*AreaW
            Fe = rho*u[j,i]*AreaE
            Fn = rho*v[j-1,i]*AreaN
            Fs = rho*v[j,i]*AreaS

            

            b[j,i] = Fw - Fe + Fs - Fn


    # Calculate x
    for i in range(1, Nx):
        for j in range(1, Ny):
            # CALCULATE SOURCE
            s[j,i] = b[j,i] + Aw[j,i]*pp[j,i-1] + Ae[j,i]*pp[j,i+1]

            if (j==1):
                s[j,i] = s[j,i] + An[j,i]*pp[j-1,i]
            elif (j==(Ny-1)):
                s[j,i] = s[j,i] + As[j,i]*pp[j+1,i]

            # ASSEMBLY MATRICES (A)
            A = Amatrix(An[1:-1,i], Ap[1:-1,i], As[1:-1,i])
            pp[1:-1,i] = TDMA(A,s[1:-1,i]).transpose()

    # Calculate x
    s[:,:] = 0
    for j in range(1, Ny):
        for i in range(1, Nx):
            # CALCULATE SOURCE
            s[j,i] = b[j,i] + An[j,i]*pp[j-1,i] + As[j,i]*pp[j+1,i]

            if (i==1):
                s[j,i] = s[j,i] + Aw[j,i]*pp[j,i-1]
            elif (i==(Nx-1)):
                s[j,i] = s[j,i] + Ae[j,i]*pp[j,i+1]

            # ASSEMBLY MATRICES (A)
            A = Amatrix(Aw[j,1:-1], Ap[j,1:-1], Ae[j,1:-1])
            pp[j,1:-1] = TDMA(A,s[j,1:-1]).transpose()

    
    return (pp, b)

def SIMPLEAlgorithm(us, vs, ps, pp, rho, mu, alpha_u, alpha_v, alpha_p, uInlet, uOutlet, nInlet, nOutlet,
                    nLoc, Lx, Ly, xu, yu, xv, yv, xc, yc, xs, ys, Nx, Ny):
    # 1. Solve discretised momentum equations -> us, vs
    un, Apu, rmu = Solveu(us, vs, ps, rho, mu, alpha_u, xc, yc, xu, yu, xv, yv, Nx, Ny)
    vn, Apv, rmv = Solvev(us, vs, ps, rho, mu, alpha_v, xc, yc, xu, yu, xv, yv, Nx, Ny)

    # 2. Solve pressure correction equation -> pp
    pp, rmp = Solvep(un, vn, Apu, Apv, rho, alpha_u, alpha_v, alpha_p, xc, yc, Nx, Ny)

    # 3. Correct pressure and velocities
    ps = CorrectContinuity(ps, pp, alpha_p, Nx, Ny)
    us, vs = CorrectMomentum(un, vn, pp, Apu, Apv, alpha_u, alpha_v, alpha_p,
                                        xu, yu, xv, yv, Nx, Ny)

    # 4. Update boundary condition
    us, vs, ps = BoundaryCondition(us, vs, ps, uInlet, uOutlet, nInlet, nOutlet, nLoc)

    # 5. Error calculation
    resCont = ErrorContinuity(Nx, Ny, rmp)
    resMomu, resMomv = ErrorMomentum(Nx, Ny, rmu, rmv)

    return (us, vs, ps, resCont, resMomu, resMomv)

def updateProgress(progress, rmseCont, rmseMomu, rmseMomv):
    barLength = 30 # Modify this to change the length of the progress bar
    status = ""

    if isinstance(progress, int):
        progress = float(progress)
    
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    
    block = int(round(barLength*progress))
    text = "\rSolving SIMPLE Algorithm      : |{0}| {1}%, (RMSE: p = {2}, u = {3}, v = {4}), {5}".format( "█"*block + " "*(barLength-block), round(100*progress,2), round(rmseCont,9), round(rmseMomu,9), round(rmseMomv,9), status)
    
    sys.stdout.write(text)
    sys.stdout.flush()


# LINEAR SOLUTION
def TDMA(A, b):
    [nA,nb]= np.shape(A)
    a      = np.array([nA,nb])
    n      = a[0]
    alpha  = np.zeros((n,1))
    beta   = np.zeros((n,1))
    D      = np.zeros((n,1))
    C      = np.zeros((n,1))
    AA     = np.zeros((n,1))
    Cprime = np.zeros((n,1))
    x      = np.zeros((n,1))

    for j in range (n):
        D [j,0]    = A[j,j]
        C [j,0]    = b[j]

        if j==0:
            alpha [j,0] = -A[j,j+1]
            beta  [j,0] = 0
            AA    [j,0] = alpha[j,0] / (D[j,0])
            Cprime[j,0] = (0+ C[j,0])/ (D[j,0] - beta[j,0]*0)

        elif j==(n-1):
            alpha[j,0] = 0
            beta [j,0] = -A[j,j-1]
            AA   [j,0] = 0
            Cprime[j,0] = (beta[j,0]*Cprime[j-1,0]+ C[j,0])/ (D[j,0] - beta[j,0]*AA[j-1,0])
        else :
            alpha[j,0] = -A[j,j+1]
            beta [j,0] = -A[j,j-1]
            AA   [j,0] = alpha[j,0] / (D[j,0] - beta[j,0]*AA[j-1,0])
            Cprime[j,0] = (beta[j,0]*Cprime[j-1,0]+ C[j,0])/ (D[j,0] - beta[j,0]*AA[j-1,0])

    for j in range (n):
        if j==0:
            x[n-1-j,0] = Cprime[n-1-j,0]

        else:
            x[n-1-j,0] = AA[n-1-j,0]*x[n-j,0] + Cprime[n-1-j,0] 

    return (x)

def Amatrix(beta, D, alpha):
    [n]    = np.shape(beta)
    a      = np.array([n])
    n      = a[0]
    A      = np.zeros((n,n))

    for j in range (n):
        for i in range (n):
            if i==j:
                A[j,i] = D[j]
            elif i==j+1:
                A[j,i] = -alpha[j]
            elif i==j-1:
                A[j,i] = -beta[j]

    return (A)

# --------------------------------------------------------------------------------- #
