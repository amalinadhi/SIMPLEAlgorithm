import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------- #
# GRID
def Grid(xc, yc, xs, ys, xu, yu, xv, yv):
    # Initialize figure
    plt.figure(figsize=(8, 8), dpi=100)
    plt.axis("equal")

    # Physical grids
    plt.scatter(xc, yc, color="k", marker="o", s=10)
    plt.plot(xc, yc, color="k")
    plt.plot(xc.transpose(), yc.transpose(), color="k")

    # Staggered grids
    plt.scatter(xv, yv, color="b", marker="^", s=10)                  # V-grid
    plt.scatter(xu, yu, color="r", marker=">", s=10)                  # U-grid
    plt.scatter(xs, ys, color="g", marker="o", s=10)                  # p-grid

    # Plot properties
    plt.title("STAGGERED GRID ARRANGEMENT")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.show()


# SOLVING PROCESS
def Residue(iters, rmseCont, rmseu, rmsev):
    minIter = np.min(iters)
    maxIter = np.max(iters)
    
    # Generating plot
    fig = plt.figure(figsize=(8, 6))
    ax = plt.subplot(111)

    plt.plot(iters, rmseCont, color='r', linestyle='-', linewidth='2', label=r'$p$')
    plt.plot(iters, rmseu, color='g', linestyle='--', linewidth='2', label=r'$u$')
    plt.plot(iters, rmsev, color='b', linestyle='-.', linewidth='2', label=r'$v$')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel("RMSE")
    plt.xlabel("Iterations")
    plt.xlim(minIter, maxIter)
    plt.yscale("log")
    plt.show()


# VISUALIZATION
def Velocity(x, y, vel, u, v, types="resultant", nLegend=21):
    # Initialize legend properties
    if (types=="resultant"):
        maxV = 1.362
        minV = 0.000
        dx = 0.136
    elif (types=="u"):
        maxV = 1.306
        minV = -0.004
        dx = 0.131
        minV = minV - dx
    elif (types=="v"):
        maxV = 0.279
        minV = -0.625
        dx = 0.09
        maxV = maxV + dx

    # Initialize figure
    fig, ax = plt.subplots(figsize=(10,8), dpi=100)
    leg = np.arange(minV, maxV, dx)

    if (types=="resultant"):
        cf = plt.contourf(x, y, vel, leg, cmap=cm.jet)         # Contour plot
        plt.title(r"RESULTANT VELOCITY, $V_{RESULTANT} (m/s)$")
    elif (types=="u"):
        cf = plt.contourf(x, y, u, leg, cmap=cm.jet)         # Contour plot
        plt.title(r"u-VELOCITY, $u, (m/s)$")
    elif (types=="v"):
        cf = plt.contourf(x, y, v, leg, cmap=cm.jet)         # Contour plot
        plt.title(r"v-VELOCITY, $v, (m/s)$")
    cbaxes = fig.add_axes()
    cb = plt.colorbar(cf, cbaxes)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    plt.xlabel(r"$X (m)$")
    plt.ylabel(r"$Y (m)$")
    
    plt.show()

def Pressure(x, y, p, nLegend=21):
    p = p - np.min(p)
    p = p*1.5816
    
    #p[:,-2] = p[:,-3]

    # Initialize legend properties
    maxV = 5.694
    minV = 0.000
    dx = 0.1

    # Initialize figure
    fig, ax = plt.subplots(figsize=(10,8), dpi=100)
    leg = np.arange(minV, maxV, dx)                         # Legend range
    
    cf = plt.contourf(x, y, p, leg, cmap=cm.jet)         # Contour plot
    #cf = plt.contourf(x, y, p, cmap=cm.jet) 
    plt.title(r"Relative Static Pressure, $p_{static}, (Pa)$")
    cbaxes = fig.add_axes()
    cb = plt.colorbar(cf, cbaxes)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    plt.xlabel(r"$X (m)$")
    plt.ylabel(r"$Y (m)$")
    
    plt.show()

def Streamline(x, y, vel, u, v, nLegend=21):
    maxV = 1.362
    minV = 0.000
    dx = 0.136
    leg = np.arange(minV, maxV+dx, dx)
    

    fig, ax = plt.subplots(figsize=(10,8), dpi=100)
    cf = plt.contourf(x, y, vel, leg, cmap=cm.jet)
    cbaxes = fig.add_axes()
    cb = plt.colorbar(cf, cbaxes)
    plt.xlabel(r"$X (m)$")
    plt.ylabel(r"$Y (m)$")
    
    u = u.astype(np.float32)
    v = v.astype(np.float32)
    vel = vel.astype(np.float32)
    
    Ny, Nx = u.shape
    for j in range(Ny):
        for i in range(Nx):
            u[j,i] = u[j,i]/vel[j,i]
            v[j,i] = v[j,i]/vel[j,i]

    ax.quiver(x,y,u,v, color="k", width=0.002)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    ax.set_aspect('equal')
    
    plt.tight_layout()
    plt.show()

def Compare(xs, ys, v_res, u_nodes, v_nodes, p_nodes, Nx, Ny):    
    #indeks = 9
    indeks = int((Nx-3)/2)


    # Horizontal
    ph = (p_nodes[indeks]+p_nodes[indeks+1])/2
    uh = (u_nodes[indeks]+u_nodes[indeks+1])/2
    #uh = u_nodes[indeks]
    vh = (v_nodes[indeks]+v_nodes[indeks+1])/2
    xh = xs[0]

    # Vertical
    pv = []
    uv = []
    vv = []
    yv = []
    for i in range(Nx-1):
        yv.append(ys[i][indeks])
        pv.append((p_nodes[i][indeks] + p_nodes[i][indeks+1])/2)
        uv.append((u_nodes[i][indeks] + u_nodes[i][indeks+1])/2)
        vv.append((v_nodes[i][indeks] + v_nodes[i][indeks+1])/2)

    yv = np.array(yv)
    pv = np.array(pv)
    uv = np.array(uv)
    vv = np.array(vv)

    # Load Data
    data = pd.read_csv("datas/ANSYS.csv", delimiter=",")
    X_AN = data["X"]
    uh_AN = data["uh"]
    vh_AN = data["vh"]
    ph_AN = data["ph"]
    Y_AN = data["Y"]
    uv_AN = data["uv"]
    vv_AN = data["vv"]
    pv_AN = data["pv"]


    data = pd.read_csv("datas/CFDSOF.csv", delimiter=",")
    X_CFD = data["X"]
    uh_CFD = data["uh"]
    vh_CFD = data["vh"]
    ph_CFD = data["ph"]
    Y_CFD = data["Y"]
    uv_CFD = data["uv"]
    vv_CFD = data["vv"]
    pv_CFD = data["pv"]


    fig = plt.figure(figsize=(9, 6))
    plt.plot(xh, uh, color='r', linestyle='-', linewidth='2', label=r'SIMPLE Algorithm')
    plt.plot(X_AN, uh_AN, color='g', linestyle='--', linewidth='1', label=r'ANSYS')
    plt.plot(X_CFD, uh_CFD, color='b', linestyle='-.', linewidth='1', label=r'CFDSOF')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.xlim(0, 0.1)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel(r"$u (m/s)$")
    plt.xlabel(r"$x (m)$")
    plt.show()

    fig = plt.figure(figsize=(9, 6))
    plt.plot(xh, vh, color='r', linestyle='-', linewidth='2', label=r'SIMPLE Algorithm')
    plt.plot(X_AN, vh_AN, color='g', linestyle='--', linewidth='1', label=r'ANSYS')
    plt.plot(X_CFD, vh_CFD, color='b', linestyle='-.', linewidth='1', label=r'CFDSOF')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.xlim(0, 0.1)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel(r"$v (m/s)$")
    plt.xlabel(r"$x (m)$")
    plt.show()

    fig = plt.figure(figsize=(9, 6))
    plt.plot(yv, uv, color='r', linestyle='-', linewidth='2', label=r'SIMPLE Algorithm')
    plt.plot(Y_AN, uv_AN, color='g', linestyle='--', linewidth='1', label=r'ANSYS')
    plt.plot(Y_CFD, uv_CFD, color='b', linestyle='-.', linewidth='1', label=r'CFDSOF')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.xlim(0, 0.1)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel(r"$u (m/s)$")
    plt.xlabel(r"$y (m)$")
    plt.show()

    fig = plt.figure(figsize=(9, 6))
    plt.plot(yv, vv, color='r', linestyle='-', linewidth='2', label=r'SIMPLE Algorithm')
    plt.plot(Y_AN, vv_AN, color='g', linestyle='--', linewidth='1', label=r'ANSYS')
    plt.plot(Y_CFD, vv_CFD, color='b', linestyle='-.', linewidth='1', label=r'CFDSOF')
    plt.legend(prop={"size":10}, loc='upper right')
    plt.xlim(0, 0.1)
    plt.minorticks_on()
    plt.grid(b=True, which='major', color='0.65', linestyle='-')
    plt.grid(b=True, which='minor', color='0.65', linestyle='--')
    plt.ylabel(r"$v (m/s)$")
    plt.xlabel(r"$y (m)$")
    plt.show()

# ---------------------------------------------------------------------------- #