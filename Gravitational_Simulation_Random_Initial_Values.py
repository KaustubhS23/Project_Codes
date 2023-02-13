import random as ran
import vpython as vp

vp.scene.background = vp.color.black

G = 1E9
dt = 0.0001
t = 0 # so that we can output time values aswell
radii = 100
# classes and functions

class Vector:

    #defines the initialisation function
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z
    #defines a function to give the magnitude of the vector
    def mag(self):
        return (((self.x)**2+(self.y)**2+(self.z)**2)**0.5)
    def printvector(self):
        return str(self.x)+"i +",str(self.y)+"j +",str(self.z)+"k"

    def unit(self):
        return(self.multiply(1/self.mag()))
    #summation of vectors:

    def sum(self,v):
        return Vector(self.x+v.x,self.y+v.y,self.z+v.z)

    #dot product:
    def dot(self,v):
        return (v.x*self.x+v.y*self.y+v.z*self.z)
    def dif(self,v):
        return(self.sum(v.multiply(-1)))
    def multiply(self,n):
        return Vector(self.x*n, self.y*n, self.z*n)


class Particle:
    def __init__(self,position,v,m,r,f):
        self.position = position
        self.v = v
        self.m = m
        self.r = r
        self.f = f
    # gives force on v due to self
    def calcforce(self,p):
        try:
            r = p.position.dif(self.position)
            F = G*self.m*p.m/(r.mag()**3)
            FF = r.multiply(F)
        except ZeroDivisionError:
            return Vector(0, 0, 0)
        return FF
    def velocity(self):
        a = (self.f.multiply(1/self.m))
        self.v = self.v.sum(a.multiply(dt))
    def displace(self):
        self.position = self.position.sum(self.v.multiply(dt))
    def collide(self,p):
        unit_id = self.position.dif(p.position).unit()
        M = self.m+p.m
        u1par = unit_id.multiply(self.v.dot(unit_id))
        u2par = unit_id.multiply(p.v.dot(unit_id))
        self.v = self.v.dif(u1par)
        p.v = p.v.dif(u2par)
        u1par.printvector()
        v1 = Vector(u1par.x,u1par.y,u1par.z)
        v2 = Vector(u2par.x,u2par.y,u2par.z)
        u2par.x = (2*self.m*v1.x+(p.m-self.m)*v2.x)/M
        u1par.x = (2*p.m*v2.x-(p.m-self.m)*v1.x)/M
        u2par.y = (2*self.m*v1.y+(p.m-self.m)*v2.y)/M
        u1par.y = (2*p.m*v2.y-(p.m-self.m)*v1.y)/M
        u2par.z = (2*self.m*v1.z+(p.m-self.m)*v2.z)/M
        u1par.z = (2*p.m*v2.z-(p.m-self.m)*v1.z)/M
        self.v = self.v.sum(u1par)
        p.v = p.v.sum(u2par)
        self.velo = vp.vector(self.v.x, self.v.y, self.v.z)
        self.f = self.calcforce(p).multiply(-1)
    
    def docollide(self,p):
        if self.position.dif(p.position).mag()<= 2*radii:
            return True
        return False

def total_force(i, n):
    sumly = Vector(0, 0, 0)
    for j in range(n):
        if j == i:
            continue
        sumly = sumly.sum(particles_list[i].calcforce(particles_list[j]))
    return sumly

def combine(i, j):
    if vpython_spheres[i].mass + vpython_spheres[j].mass == 0:
        return None
    vpython_spheres.append(vp.sphere(pos=(vpython_spheres[i].pos + vpython_spheres[j].pos) / 2, mass=vpython_spheres[i].mass + vpython_spheres[j].mass, radius=radii, velo=(vpython_spheres[i].velo * vpython_spheres[i].mass + vpython_spheres[j].velo * vpython_spheres[j].mass) / (vpython_spheres[i].mass + vpython_spheres[j].mass), color=vp.color.blue, make_trail=False))
    particles_list.append(Particle(particles_list[i].position.sum(particles_list[j].position).multiply(0.5), (particles_list[i].v.multiply(particles_list[i].m).sum(particles_list[j].v.multiply(particles_list[j].m))).multiply(1/(particles_list[i].m + particles_list[j].m)), particles_list[i].m + particles_list[j].m, 50, particles_list[i].f.sum(particles_list[j].f)))
    particles_list[i].m = 0
    particles_list[j].m = 0
    particles_list[i].position = Vector(ran.uniform(100000,200000), ran.uniform(-100000, -200000), ran.uniform(-100000, -200000))
    particles_list[j].position = Vector(ran.uniform(100000,200000), ran.uniform(-100000, -200000), ran.uniform(-100000, -200000))
    particles_list[i].f = Vector(0, 0, 0)
    particles_list[j].f = Vector(0, 0, 0)
    vpython_spheres[i].visible = False
    vpython_spheres[j].visible = False
    vpython_spheres[i].mass = 0
    vpython_spheres[j].mass = 0
    vpython_spheres[i].pos = vp.vector(ran.uniform(100000, 200000), ran.uniform(-100000, -200000), ran.uniform(-100000, -200000))
    vpython_spheres[j].pos = vp.vector(ran.uniform(100000, 200000), ran.uniform(-100000, -200000), ran.uniform(-100000, -200000))

def initialize_particles(n, x_pos, y_pos, z_pos, input_mass, velocity_x, velocity_y, velocity_z):
    vpython_spheres[n] = vp.sphere(pos=vp.vector(x_pos, y_pos, z_pos), velo=vp.vector(velocity_x, velocity_y, velocity_z), mass=input_mass, radius=radii, color=vp.color.white, make_trail=False)
    particles_list[n] = Particle(Vector(x_pos, y_pos, z_pos), Vector(velocity_x, velocity_y, velocity_z), input_mass, 50, Vector(0, 0, 0))



# getting the inputs and setting up the objects
opt = input("1.elastic\n2.inelastic\nSelect 1 or 2 >>> ")
noh = int(input("n ___::::>>> "))
timee = float(input("the time uptill which you want the objects to move::>> "))  
vpython_spheres = [0 * i for i in range(noh)]
p_force = [0 * i for i in range(noh)]
particles_list = [0*i for i in range(noh)]
#name = input('give the name of the file: ')
#file = open(name, 'r')
for i in range(noh):
    x=ran.uniform(-1000,1000)
    y=ran.uniform(-1000,1000)
    z=ran.uniform(-1000,1000)
    xd=ran.uniform(-100,100)
    yd=ran.uniform(-100,100)
    zd=ran.uniform(-100,100)
    m=10
    initialize_particles(i,x,y,z,10,xd,yd,zd)
'''
for li in file:
    li.rstrip()
    liis = li.split()
    initialize_particles(int(liis[0]), float(eval(liis[1])), float(eval(liis[2])), float(eval(liis[3])), float(eval(liis[4])), float(eval(liis[5])), float(eval(liis[6])), float(eval(liis[7])))
'''
#outputfile = open(name +'_output.txt', '+w')

'''
x_axis = vp.arrow(pos=vp.vector(0,0,0), axis=vp.vector(1000, 0, 0), shaftwidth=1, color = vp.color.red)
y_axis = vp.arrow(pos=vp.vector(0,0,0), axis=vp.vector(0, 1000, 0), shaftwidth=1, color= vp.color.green)
z_axis = vp.arrow(pos=vp.vector(0,0,0), axis=vp.vector(0, 0, 1000), shaftwidth=1, color = vp.color.blue)
'''

# the loop for calculating and displaying
while True:
    vp.rate(5000)
    for i in range(noh):
        particles_list[i].f = total_force(i, noh)
        for j in range(int(noh)):
            if i == j:
                continue
            if particles_list[i].docollide(particles_list[j]):
                if opt == '1':
                    particles_list[i].collide(particles_list[j])
                elif opt == '2':
                    combine(i, j) 
                    noh += 1
    for i in range(noh):
        try:
            particles_list[i].velocity()
            particles_list[i].displace()
        except ZeroDivisionError:
            continue
    for i in range(noh):
        vpython_spheres[i].momentum = vp.vector(particles_list[i].v.x, particles_list[i].v.y, particles_list[i].v.z)
        vpython_spheres[i].pos = vp.vector(particles_list[i].position.x, particles_list[i].position.y, particles_list[i].position.z)
        t += dt
        #outputfile.write('{}, {}, {} \n'.format(i, particles_list[i].position.printvector(), particles_list[i].v.printvector()))
    if t >= timee:
        break
#file.close()
#outputfile.close()
print('finished')
