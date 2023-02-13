import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits import mplot3d
from matplotlib import cm

def psi1(x):
    return(np.exp(-(x/2)**2)/4)
def psi2(x):
    return(0)
def psi3(x):
    return(np.sin(x))

def v1(x):
    return(x)
def v2(z):
    return(2)
def v3(z):
    return(0)

def solve_explicit(end_t=20,type_init=1,bc_type='Dirichlet',d_bc=[0,0],psi_type=1,v_type=1,xi=-10,xf=10):
    '''
    Solves the Convection - Diffusion equation in one dimension numerically by the Finite Difference Method
    
    end_t : Time till which the solution is computed (starting from 0)
    init_c : Initial Conditions (The value of the temprature at every discreete x value), given in the form of a 1*n array (n is the number of discrete x values which is (xi-xf)/h_x))
    xi : First x value at which the value of temprature is calculated
    xf : Last x value at which the value of temprature is caluculated
    psi : Function of position which returns the value of the source term (Time independent)
    v : The convection field. Function which returns the value of the velocity of the medium at a certain position. (Time independent)
    bc_type : The kind of boundary condition that is needed. If Dirichlet boundary condition is chosen, then d_bc will be a list containing the values at xi and xf respectively.
              If Neumann boundary condition is chosen, then the case will be when the heat flux is zero at both ends.

    Returns the value of temprature at every point x domain for every discrete time value in the form of an array,
    in which the rows have the value of tempratures at every x value, for a fixed time, and successive rows have tempratures for successive time values ( y[time index,position index] )
    '''
    h_x=1E-1
    h_t=1E-2
    alpha=1/6

    x_ran=np.arange(xi,xf+.1,h_x)
    
    k1=1
    k2=.5
    if type_init==1:
        init_c=k1*np.exp(-(x_ran*k2)**2)
    elif type_init==2:
        init_c=np.sin(x_ran*np.pi/(2))
    elif type_init==3:
        init_c=np.sin(x_ran*np.pi)

    if psi_type==1:
        psi=psi1
    elif psi_type==2:
        psi=psi2
    elif psi_type==3:
        psi=psi3

    if v_type==1:
        v=v1
    elif v_type==2:
        v=v2
    elif v_type==3:
        v=v3
    
    #len(init_c) would serve as the number of discrete values of x
    ndt=int(end_t/h_t)
    #ndt would be the number of discrete values of t, for a small enough h_t, end_t/h_t will be an integer
    
    y=np.zeros([ndt,len(init_c)])
    #a row of y would contain the value of temprature at increasing value of x for constant time. different rows correspond to different values of time.
    
    y[0]=init_c
    #Setting the time at t=0 to be the initial conditions
    
    if bc_type=="Dirichlet":
        y[:,0]=d_bc[0]
        y[:,-1]=d_bc[1]
        for t in range(1,ndt):
            #t takes all the indeces of all the dicreetized time value (except t=0, since we already know what the value of temprature is at all x at that time)
        
            for x in range(1,len(init_c)-1):
                #x takes the values of all the indices of the discreetized x values (except the ones at the boundaries)
            
                y[t,x]=alpha*( (y[t-1,x+1] - 2*y[t-1,x] + y[t-1,x-1]) /(h_x**2) + psi(x*h_x+xi) -v(x*h_x+xi)*( y[t-1,x+1] - y[t-1,x-1] )/(2*h_x) -( v((x+1)*h_x+xi) - v((x-1)*h_x+xi) )*y[t-1,x]/(2*h_x) )*h_t + y[t-1,x]
           
    if bc_type=="Neumann":
        for t in range(1,ndt):        
            for x in range(1,len(init_c)-1):
            
                y[t,x]=alpha*( (y[t-1,x+1] - 2*y[t-1,x] + y[t-1,x-1]) /(h_x**2) + psi(x*h_x+xi) -v(x*h_x+xi)*( y[t-1,x+1] - y[t-1,x-1] )/(2*h_x) -( v((x+1)*h_x+xi) - v((x-1)*h_x+xi) )*y[t-1,x]/(2*h_x) )*h_t + y[t-1,x]

                y[t,0]=(y[t-1,0]+y[t-1,1])/2
                y[t,-1]=(y[t-1,-1]+y[t-1,-2])/2
                #Setting the values of the function at the boundary to be the average of themselves and their immidiate neighbour. 

    fig, ax = plt.subplots()
    plt.grid()
    l, = plt.plot(x_ran, init_c, lw=2,color='r')
    plt.xlabel("Position")
    plt.ylabel("Temprature")

    ax = plt.axis([xi,xf,-1,1])

    axamp = plt.axes([0.25, .005, 0.50, 0.02])
    
    samp = Slider(axamp, 'Time', 0, end_t, valinit=0)

    def update(val):
    
        ct = samp.val
        # ct is the current value of the slider
    
        itr=int((ct-.001)/h_t)
        l.set_ydata(y[itr,:])
        # update curve
    
        fig.canvas.draw_idle()
        # redraw canvas while idle
    
    samp.on_changed(update)
    # call update function on slider value change

    plt.show()
    

def solve_CN(end_t=20,type_init=2,d_bc=[0,0],psi_type=1,v_type=1,xi=-10,xf=10):
    '''
    Solves the Convection - Diffusion equation in one dimension numerically by the Finite Difference Method
    
    end_t : Time till which the solution is computed (starting from 0)
    init_c : Initial Conditions (The value of the temprature at every discreete x value), given in the form of a 1*n array (n is the number of discrete x values which is (xi-xf)/h_x))
    xi : First x value at which the value of temprature is calculated
    xf : Last x value at which the value of temprature is caluculated
    psi : Function of position which returns the value of the source term (Time independent)
    v : The convection field. Function which returns the value of the velocity of the medium at a certain position. (Time independent)
    bc_type : The kind of boundary condition that is needed. If Dirichlet boundary condition is chosen, then d_bc will be a list containing the values at xi and xf respectively.
              If Neumann boundary condition is chosen, then the case will be when the heat flux is zero at both ends.

    Returns the value of temprature at every point x domain for every discrete time value in the form of an array,
    in which the rows have the value of tempratures at every x value, for a fixed time, and successive rows have tempratures for successive time values ( y[time index,position index] )
    '''
    def tridiag(A):               
    #uses a slightly modified dollittle algorithm
    
        n=len(A)
        c=np.zeros([n-1,1])            
        d=np.zeros([n,1])
        e=np.zeros([n-1,1])
    
        d[0,0]=A[0,0]
        for i in range(1,n):        
            c[i-1]=A[i,i-1]
            d[i]=A[i,i]
            e[i-1]=A[i-1,i]
            
        #initializing the coeffecient columns
        
        for k in range(1,n):
            d[k]+= -c[k-1]*e[k-1]/d[k-1]
            #d has the diagnoal elements of U,e remains unchanged
        
            c[k-1]=c[k-1]/d[k-1]
            #since the c column would become zero, we instead store the lamda values in it
        return(c,d,e)

    def sub(b,c,d,e):
        n=len(b)
        y=np.zeros([n,1])
        y[0]=b[0] 
        for i in range(1,n):
            y[i]=b[i]-y[i-1]*c[i-1]
        x=np.zeros([n,1])
        x[n-1]=y[n-1]/d[n-1]
        #forward substitution phase
    
        for i in range(n-2,-1,-1):
            x[i]=(y[i]-e[i]*x[i+1])/d[i]
        #backward substitution phase

        return(x)
    
    h_x=1E-1
    h_t=1E-1
    alpha=1/6

    x_ran=np.arange(xi,xf+.1,h_x)
    
    k1=1
    k2=.5
    if type_init==1:
        init_c=k1*np.exp(-(x_ran*k2)**2)
    elif type_init==2:
        init_c=np.sin(x_ran*np.pi/(2))
    elif type_init==3:
        init_c=np.sin(x_ran*np.pi)

    if psi_type==1:
        psi=psi1
    elif psi_type==2:
        psi=psi2
    elif psi_type==3:
        psi=psi3

    if v_type==1:
        v=v1
    elif v_type==2:
        v=v2
    elif v_type==3:
        v=v3
    
    ndt=int(end_t/h_t)
    #ndt would be the number of discrete values of t, for a small enough h_t, end_t/h_t will be an integer

    ndx=len(init_c)
    #ndx would be the number of discrete values of x, for a small enough h_t, end_t/h_x will be an integer

    y=np.zeros([ndt,len(init_c)])
    #a row of y would contain the value of temprature at increasing value of x for constant time. different rows correspond to different values of time.
    
    y[0]=init_c
    #Setting the time at t=0 to be the initial conditions
    
    y[:,0]=d_bc[0]
    y[:,-1]=d_bc[1]
    
    A=np.zeros([ndx,ndx])
    zeta=h_t/(4*h_x)
    xsi=alpha*h_t/(2*h_x**2)
    
    #filling in A
    A[0,0]=1
    A[-1,-1]=1
    
    x=xi
    for i in range(1,len(A)-1):
        x+=h_x
        A[i,i-1]=-xsi-zeta*v(x-h_x)
        A[i,i]=1+2*xsi
        A[i,i+1]=zeta*v(x+h_x)-xsi
        

    c,d,e=tridiag(A)
    #The beauty of LU decompostion comes into play here. Since A is not changing (Steady State case), the decompostion would remain the same for every time step. Only b would change.
    #Therefore we can use the same L and U matrices over and over again, and we only have to so the forward and the back substitutions in the main loop.

    b=np.zeros([ndx,1])
    b[0,0]=y[0,0]
    b[-1,0]=y[0,-1]

    for t in range(1,ndt):
        x=xi
        for i in range(1,ndx-1):
            x+=h_x
            b[i,0]=psi(x)*h_t/2 + ( (xsi+zeta*v(x-h_x))*y[t-1,i-1] + (1-2*xsi)*y[t-1,i] + (xsi-zeta*v(x+h_x))*y[t-1,i+1] )
        x=sub(b,c,d,e)

        y[t,:]=x.T        

    fig, ax = plt.subplots()
    plt.grid()
    
    l, = plt.plot(x_ran, init_c, lw=2,color='r')

    ax = plt.axis([xi,xf,-1,1])

    axamp = plt.axes([0.25, .005, 0.50, 0.02])
    
    samp = Slider(axamp, 'Time', 0, end_t, valinit=0)

    def update(val):
        # ct is the current value of the slider
        ct = samp.val
        # update curve
        itr=int((ct-.001)/h_t)
        l.set_ydata(y[itr,:])
        # redraw canvas while idle
        fig.canvas.draw_idle()

    # call update function on slider value change
    samp.on_changed(update)

    plt.show()

def solve_3d(end_t=10,type_init=2,type_bc=1123,xi=-2,xf=2):
    
    alpha=1/6
    
    hx=1E-1
    nx= int((xf-xi)/hx)
    #number of grid points in x (and the y) axes alone.

    ht=1E-2
    nt=int(end_t/ht)

    mesh_range = np.arange(xi, xf, hx)
    temp=np.ones_like(mesh_range)
    
    x, y = np.meshgrid(mesh_range, mesh_range,sparse=False)

    # Initial condition and the boundary conditions
    if type_init==1:
        init_u = np.exp(-5 *(x**2 + y**2))
    elif type_init==2:
        init_u=np.sin(4*x**2+4*y**2)/(x**2+y**2)
    elif type_init==3:
        init_u=np.sin(4*x**2+4*y**2)
        
    if type_bc==1:
        bc1=temp #left wall
        bc2=0*temp #bottom wall
        bc3=temp #right wall
        bc4=0*temp #top wall

        init_u[:,0]=bc1
        init_u[-1,:]=bc2
        init_u[:,-1]=bc3
        init_u[0,:]=bc4

    if type_bc==2:
        bc1=0*abs(mesh_range) #left wall
        bc2=0*abs(mesh_range) #bottom wall
        bc3=0*abs(mesh_range) #right wall
        bc4=0*abs(mesh_range) #top wall
        
        init_u[:,0]=bc1
        init_u[-1,:]=bc2
        init_u[:,-1]=bc3
        init_u[0,:]=bc4

    elif type_bc==3:
        bc1=1*abs(mesh_range) #left wall
        bc2=0*abs(mesh_range) #bottom wall
        bc3=1*abs(mesh_range) #right wall
        bc4=0*abs(mesh_range) #top wall

        init_u[:,0]=bc1
        init_u[-1,:]=bc2
        init_u[:,-1]=bc3
        init_u[0,:]=bc4

    def pde_step(u):
        '''
        given the tempratutre distribution u at time t, returns the temprature distribution at t+ht 
        '''
        #for the time being, the given boundary conditions are not used, the values of the initial conditions at the boundaries is used as the bc instaed
        ut=u.copy()
        for i in range(2,len(u)-2):
            for j in range(2,len(u)-2):
                ut[i,j]=alpha*(u[i+2,j]+u[i-2,j]+u[i,j-2]+u[i,j+2]-4*u[i,j])*ht/(2*hx)**2+u[i,j]
        return(ut)
        
    def draw_plot(x, y, u):
        ax.clear()
        ax.set_zlim(-2.01, 2.01)
        ax.plot_surface(x,y,u,cmap=cm.coolwarm,antialiased=True)

        plt.pause(1e-5)
    

    plt.ion()
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_zlim(-2.01, 2.01)
    ax.set_xlabel('X',fontsize=20)
    ax.set_ylabel('Y',fontsize=20)
 

    draw_plot(x, y, init_u)

    u_step = init_u
    for t in range(nt):
        u_step = pde_step(u_step)

        draw_plot(x, y, u_step)

#*********************************************************** Main code complete, Menu Driven Part starts ***************************************************************************

while True:
    dim=str(input("1d or 2d : "))
    
    if dim=='1d':
        method=str(input("CN or Explicit (E) : "))
        
        if method=='CN':
            psi_type=int(input('psi type : '))
            v_type=int(input("v type : "))
            bc=eval(input("Enter the boundary conditions : "))
            ic_type=int(input("Type of initial conditions : "))
            end_time=float(input("Enter the time till which the solution will be calculated : "))
            
            solve_CN(end_time,ic_type,bc,psi_type,v_type)
            
        elif method=='E':
            psi_type=int(input('psi type : '))
            v_type=int(input("v type : "))            
            bc_type=str(input("Dirichlet (D) or Neumann (N) : "))
            
            if bc_type=='D':
                bc_type="Dirichlet"
                bc=eval(input("Enter the boundary conditions : "))
                
            elif bc_type=='N':
                bc_type='Neumann'
                bc=[0,0]
                
            ic_type=int(input("Type of initial conditions : "))
            end_time=float(input("Enter the time till which the solution will be calculated : "))
            
            solve_explicit(end_time,ic_type,bc_type,bc,psi_type,v_type)
            
    elif dim=='2d':
        end_time=float(input("Enter the time till which the solution will be calculated : "))
        ic_type=int(input("Type of initial conditions : "))
        bc_type=int(input("Type of boundary conditions : "))

        solve_3d(end_time,ic_type,bc_type)
    done=str(input("Are we done? (y/n) : "))
    if done=='y':
        break
    else:
        continue

