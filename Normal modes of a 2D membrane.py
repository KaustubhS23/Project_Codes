import matplotlib.pyplot as plt
from numpy import sin, cos, linspace, meshgrid, pi
from scipy.special import jn, jn_zeros

def circnm_spt(n,m,a=1):
    '''
    Parameters
    ----------
    m : Positive integer
        Radial index of the normal mode.
    n : Non negative integer
        Angular index of the normal mode.
    a : Positive float
        The radius of the circular rim.
    
    Returns
    -------
    Two square grids and the spatial component of the solution

    '''
    
    # Create the mesh in polar coordinates.
    r = linspace(0, a, 50)
    theta = linspace(0, 2*pi, 50)
    R, THETA = meshgrid(r, theta)

    # Express the mesh in the cartesian system.
    X, Y = R*cos(THETA), R*sin(THETA)
   
    lst_zro=jn_zeros(n,m)
    #List of the first m positive zeros of the Bessel function of the first kind of order n.
    j_mn=lst_zro[-1]
    
    spt=cos(n*THETA)*jn(n,j_mn*R/a)
    #The spatial part of the normal mode. It is time independent
    #It has been computed and stored outisde the plotting loop so that we can save on computational resources.
    
    return(X,Y,spt)

def circnm_tdep(n,m,a,v,A):
    '''
    Parameters
    ----------
    n : Non negative integer
        The index for the angular part of the solution.
    m : Positive integer
        The index for the radial part of the solution.
    a : Positive float
        Radius of the circular rim.
    v : Positive float
        'Speed of oscillation' in the differential equation.
    A : Positive float
        Directly proportional to the maximum amplitude of the normal mode

    Returns
    -------
    The function which governs the time dependent part of the solution.

    '''
    lst_zro=jn_zeros(n,m)
    #List of the first m positive zeros of the Bessel function of the first kind of order n.
    j_mn=lst_zro[-1]
    def T(t):
        return(A*sin(j_mn*v*t/a))
    return(T)
    

def draw(X,Y,spt,t_dep):
    '''
    Parameters
    ----------
    X : Array
        A square numpy array which contains the X values of the grid.
    Y : Array
        A square numpy array which contains the Y values of the grid.
    spt : Array
        A square numpy array which containes the values of the spatial component of the solutions (Time independent part)
    t_dep : Function
        The function which governs the time dependent part of the solution.

    Returns
    -------
    None.
    Animates the solution of the 2D wave equation which has spt as the spatial part and t_dep as the time dependent part.
    '''
    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    def draw_plot(x, y, u):
        ax.clear()
        ax.set_zlim(-2.01, 2.01)
        ax.plot_surface(x,y,u,cmap='inferno',antialiased=True)
        ax.set_xlabel(r'$X$')
        ax.set_ylabel(r'$Y$')
        ax.set_zlabel(r'$Amp$')
        plt.pause(1E-5)
    nt=int(10E3)
    t=0
    dt=.005
    for i in range(nt):
        u_step=t_dep(t)*spt
        t+=dt
        draw_plot(X, Y, u_step)
    
    plt.show()
    
def circ_nm(n,m,a=1,v=5,A=2):
    '''
    Compiles all the previous functions.

    Parameters
    ----------
    n : Non negative integer
        The index for the angular part of the solution.
    m : Positive integer
        The index for the radial part of the solution.
    a : Positive float, optional
        Radius of the circular rim. The default is 1.
    v : Positive float, optional
        'Speed of oscillation' in the differential equation. The default is 5.
    A : Positive float, optional
        Directly proportional to the maximum amplitude of the normal mode. The default is 2.

    Returns
    -------
    None.
    Animates the normal mode with these values.
    '''
    X,Y,spt=circnm_spt(n,m,a)
    tdep=circnm_tdep(n,m,a,v,A)
    draw(X,Y,spt,tdep)
    
def circ_superpos(nm1,nm2,a=1,v=5,A=2):
    '''
    Animates the superposition of two normal modes of a vibrating circular membrane

    Parameters
    ----------
    nm1 : List
        Containes the values of n and m for the first normal mode.
    nm2 : TYPE
        Containes the values of n and m for the second normal mode.
    a : Positive float, optional
        Radius of the circular rim. The default is 1.
    v : Positive float, optional
        'Speed of oscillation' in the differential equation. The default is 5.
    A : TYPE, optional
        Directly proportional to the maximum amplitude of the normal mode. The default is 2.

    Returns
    -------
    None.

    '''
    n1=nm1[0]
    m1=nm1[1]
    n2=nm2[0]
    m2=nm2[1]
    X,Y,spt1=circnm_spt(n1,m1,a)
    X,Y,spt2=circnm_spt(n2,m2,a)
    tdep1=circnm_tdep(n1,m1,a,v,A)
    tdep2=circnm_tdep(n2,m2,a,v,A)
    spt=spt1+spt2
    def tdep(t):
        return(tdep1(t)+tdep2(t))
    draw(X,Y,spt,tdep)

def rectnm_spt(n,m,a,b):
    '''
    Parameters
    ----------
    n : Positive integer
        Index of the normal mode
    m : Positve integer
        Index of the normal mode.
    a : Positive float
        The length of the rectangular rim.
    b : Positive float
        Width of the rectangular rim.
    
    Returns
    -------
    The spatial component of the solution and two square grids specifying the base.

    '''
    
    # Create the mesh in polar coordinates.
    x = linspace(0, a, 50)
    y = linspace(0, b, 50)
    X, Y = meshgrid(x, y)
    
    spt=sin(n*pi*X/(a))*sin(m*pi*Y/(b))
    #The spatial part of the normal mode. It is time independent
    #It has been computed and stored outisde the plotting loop so that we can save on computational resources.

    return(X,Y,spt)

def rectnm_tdep(n,m,a,b,v,A):
    '''
    Parameters
    ----------
    n : Positive integer
        Index of the normal mode
    m : Positve integer
        Index of the normal mode.
    a : Positive float
        The length of the rectangular rim.
    b : Positive float
        Width of the rectangular rim.
    v : Positve float
        'The speed of oscillation' in the differential equation.
    A : Positve float
        Maximum amplitude of the normal mode.

    Returns
    -------
    The time dependent part of the solution

    '''
    def T(t):
        return( A*sin(pi*v*t*((m/b)**2+((n/a)**2)**0.5)))
    return(T)

def rect_nm(n,m,a=1,b=2,v=5,A=2):
    '''
    Plots the specified normal mode with a rectangular base.

    Parameters
    ----------
    n : Positive integer
        Index of the normal mode
    m : Positve integer
        Index of the normal mode.
    a : Positive float
        The length of the rectangular rim. The default is 1.
    b : Positive float
        Width of the rectangular rim. The default is 2.
    v : Positve float, optional
        'The speed of oscillation' in the differential equation. The default is 5.
    A : Positve float, optional
        Maximum amplitude of the normal mode. The default is 2.

    Returns
    -------
    None.

    '''
    X,Y,spt=rectnm_spt(n,m,a,b)
    tdep=rectnm_tdep(n,m,a,b,v,A)
    draw(X,Y,spt,tdep)
    
def rect_superpos(nm1,nm2,a=1,b=2,v=5,A=.5):
    '''
    Plots the superposition of two normal modes.
    
    Parameters
    ----------
    nm1 : List
        Contains the values of n and m for the first normal mode.
    nm2 : List
        Containes the values of n and m for the second normal mode.
    a : Positve float, optional
        The length of the rectangle. The default is 1.
    b : Positive float, optional
        The width of the rectangle. The default is 2.
    v : Positve float, optional
        The 'speed of oscillation' in the differential equation. The default is 5.
    A : Positive float, optional
        The maximum amplitude of the solution. The default is 2.

    Returns
    -------
    None.

    '''
    n1=nm1[0]
    m1=nm1[1]
    n2=nm2[0]
    m2=nm2[1]
    
    X,Y,spt1=rectnm_spt(n1,m1,a,b)
    X,Y,spt2=rectnm_spt(n2,m2,a,b)
    
    tdep1=rectnm_tdep(n1,m1,a,b,v,A)
    tdep2=rectnm_tdep(n2,m2,a,b,v,A)
    spt=spt1+spt2
    
    def tdep(t):
        return(tdep1(t)+tdep2(t))
    draw(X,Y,spt,tdep)
    

rect_nm(1,4)
#circ_nm(2,3)

#rect_superpos([1,3],[5,6])
#circ_superpos([0,1],[6,6])
