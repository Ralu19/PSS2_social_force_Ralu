t = 0 # start time
Tf = 40 # end time
# variables for modeling
tau = 0.5 
mass = 80
F = 2000
Fwall = 2000
lambda_ = 0.5
delta = 0.08
kappa = 120000
eta = 240000
dmax = 0.1 # max distance to take in account for contacts
drawper = 5000 # generate plot for 1 per 1000 iterations of dt

nn = 10 # number of people
box = [120,130,10,20] # coordinates of the box that will be populated [xmin, xmax, ymin, ymax]
dest_name = "door" # name given to by domain.add_destination function
radius_distribution = ["uniform",0.4,0.6] # distribution variable 
velocity_distribution = ["normal",1.2,0.1] # distribution varible
rng = 0 # some seed value for the distribution, if =0 then random value will be chosen on run
dt = 0.0005 # timestep
dmin_people=0 # minimal disired distance to other people
dmin_walls=0 # minimal disired distance to walls

#need to be intitailized to play nice
draw = False
cc = 0
itermax=10 # max number of uzawa projectsions, only intressting that is used as projection method
counter = 0

# sys and Os are used to interact with the terminal and filesystem
import sys, os
# a trick to add the root folder of the project to the list of folders that will be searched when importing libs
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..',)
sys.path.append( mymodule_dir )

from cromosim import *
from cromosim.micro import *
from matplotlib.patches import Circle
from matplotlib.lines import Line2D

plt.ion()

#create a domain object from picture walls_trainstation
dom = Domain(name = 'trainstation', background = 'walls_trainstation.png', pixel_size = 0.1)
## To define the color for the walls
wall_color = [0,0,0]

## To define the color for the issue of the room
door_color = [255,0,0]


dom.build_domain()
dom.plot(id=3, title="Domain", savefig=False)
## To create a Destination object towards the door
dest = Destination(name='door', colors=[door_color],
                   excluded_colors=[wall_color])
dom.add_destination(dest)



dom.plot_wall_dist(id=1, step=20,
    title="Distance to walls and its gradient",
    savefig=False, filename="room_wall_distance.png")

dom.plot_desired_velocity('door',id=2, step=20,
    title="Distance to the destination and desired velocity",
    savefig=False, filename="room_desired_velocity.png")

#intialize people

groups = [{"nb":nn, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, "box":box, "destination":dest_name}] #create dict bundeling above values
# has to be a list of a dict for some reason
#curseddatatype
people = people_initialization(dom, groups, dt, dmin_people, dmin_walls, rng, itermax, projection_method='cvxopt')
contacts = None
colors = people["xyrv"][:,2]
plot_people(0,dom,people,contacts,colors)
plt.show()


# main calculating loop

while(t<Tf):
    print("\n===> Time = "+str(t))
    print("===> Compute desired velocity for domain ",name)
    I, J, Vd = dom.people_desired_velocity(people["xyrv"], people["destinations"])
    people["Vd"] = Vd
    people["I"] = I
    people["J"] = J

    print("===> Compute social forces for domain ",name)  
    contacts = compute_contacts(dom, people["xyrv"], dmax)
    print("     Number of contacts: ",contacts.shape[0])
    Forces = compute_forces(F, Fwall, people["xyrv"], contacts, people["Uold"], Vd, lambda_, delta, kappa, eta)            
    nn = people["xyrv"].shape[0]
    people["U"] = dt*(Vd[:nn,:]-people["Uold"][:nn,:])/tau + people["Uold"][:nn,:] + dt*Forces[:nn,:]/mass

    people, sensors = move_people(t, dt,people,sensors = {})
    #people = people_update_destination(people["xyrv"],domains = {"dom"},dom.pixel_size)

    people["Uold"] = people["U"]

    if(draw):
        colors =  people["xyrv"][:,2]
                ## coloring people according to their destinations
                # colors = np.zeros(all_people[name]["xyrv"].shape[0])
                # for i,dest_name in enumerate(all_people[name]["destinations"]):
                #     ind = np.where(all_people[name]["destinations"]==dest_name)[0]
                #     colors[ind]=i
        plot_people(20, dom, people, contacts,
                            colors, time=t,
                            plot_people=True, plot_contacts=False,
                            plot_paths=True, plot_velocities=False,
                            plot_desired_velocities=False, plot_sensors=False, savefig=False,
                            filename = "results/"+"domain_trainstation" +str(counter).zfill(6)+".png")
        plt.pause(0.01)

    t += dt
    cc +=1
    counter += 1
    if (cc>=drawper):
        draw = True
        cc = 0
        
    else:
        draw = False

plt.ioff()
plt.show()
sys.exit()