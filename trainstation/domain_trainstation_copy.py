t = 0 # start time
Tf = 30 # end time
# variables for modeling
tau = 0.5 
mass = 80
F = 2000
Fwall = 2000
lambda_ = 0.5
delta = 0.08
kappa = 120000
eta = 240000
dmax = 0.1 # max distance to take ,vin account for contacts
drawper = 1000 # generate plot for 1 per 1000 iterations of dt

nn = 10 # number of people
# N_1 = 10 # stationary group in zone 1 
# N_2 = 30 # stationary group in zone 2
# N_3 = 40 # stationary group in zone 3
# N_4 = 20 # stationary group in zone 4
# N_stationary = N_1 + N_2 + N_3 + N_4 
N_stationary = 10 
ntot = nn + N_stationary
box = [120,130,10,20] # coordinates of the box that will be populated [xmin, xmax, ymin, ymax]
dest_name = "door" # name given to by domain.add_destination function
radius_distribution = ["uniform",0.4,0.6] # distribution variable 
velocity_distribution = ["normal",1.2,0.1] # distribution varible
rng = 0 # some seed value for the distribution, if =0 then random value will be chosen on run
dt = 0.001 # timestep
dmin_people=0.01 # minimal disired distance to other people
dmin_walls=0.01 # minimal disired distance to walls

## DZ stuff
Fdz = 575
dzy_up = 51.0
dzy_down = 5.6

## ZONES 
## we can adjust the boxes how we like
zone_1 = [1.0, 50.2 , 1.0, 52.3 ]
zone_2 = [50.21, 125.0, 1.0, 52.3]
zone_3 = [[125.01, 216.0 , 1.0, 52.3]]
zone_4 =[[216.01 , 296.0, 1.0, 52.3]]

## Spawning zones

escalators_down = [11,33,14.1,17.9]
escalators_up = [5.5,33,38.7,41.7]
stairs = [6,33,20,36.7]


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


## Awareness stuff ##
#The awareness values should be between 0 and 1
mean_awr = 0.5 # mean of awareness
stdev_awr = 0.2 # stdev of awareness
awr = np.random.normal(mean_awr, stdev_awr, nn+N_stationary) # awareness[i] is the awareness value of person i
for i in range(awr.shape[0]): #Just in case you get unlucky
    if awr[i] < 0:
        awr[i] = 0
    elif awr[i] > 1:
        awr[i] = 1
#awr = np.ones(nn) #for testing 

#create a domain object from picture walls_trainstation
dom = Domain(name = 'trainstation', background = 'trainstation/trainstation2.png', pixel_size = 0.2)
## To define the color for the walls
wall_color = [0,0,0]

## To define the color for the issue of the room
door_color = [255,0,0]

## To add an obstacle using a matplotlib shape colored with wall_color :
##     Circle( (center_x,center_y), radius )

circle = Circle((40,2.0), 1.0)
dom.add_shape(circle, outline_color=door_color, fill_color=door_color)




dom.build_domain()
#dom.plot(id=3, title="Domain", savefig=False)
## To create a Destination object towards the door
dest = Destination(name='door', colors=[door_color],
                   excluded_colors=[wall_color])
dom.add_destination(dest)
#dom.plot_wall_dist(id=1, step=20,title="Distance to walls and its gradient",savefig=False, filename="room_wall_distance.png")

#dom.plot_desired_velocity('door',id=2, step=20,title="Distance to the destination and desired velocity",savefig=False, filename="room_desired_velocity.png")



# create dict for sensors
sensors = [{'name':'dzdown', 'domain':'dom','line':[1.0, 4.5, 296.0, 4.5],'id':[],'xy':[],'dir':[],'times':[]}]
sensors.append({'name':'dzup', 'domain':'dom','line':[1.0, 51.0, 296.0, 51.0],'id':[],'xy':[],'dir':[],'times':[]})

#intialize people

groups = [{"nb":nn, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, "box":stairs, "destination":dest_name}] #create dict bundeling above values
# has to be a list of a dict for some reason
#curseddatatype
people = people_initialization(dom, groups, dt, dmin_people, dmin_walls, rng, itermax, projection_method='cvxopt')
contacts = None
colors = people["xyrv"][:,2]

plt.show()
# group1 = [{"nb":N_1, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, "box":zone_1, "destination":dest_name}]
# people1 = people_initialization(dom, group1,dt,dmin_people,dmin_walls,rng,itermax,projection_method='cvxopt')

group2 = [{"nb":N_stationary, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, "box":[40,70,10,15], "destination":dest_name}]
people2 = people_initialization(dom, group2,dt,dmin_people,dmin_walls,rng,itermax,projection_method='cvxopt')

# group3 = [{"nb":N_3, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, "box":zone_3, "destination":dest_name}]
# people3 = people_initialization(dom, group3,dt,dmin_people,dmin_walls,rng,itermax,projection_method='cvxopt')

# group4 = [{"nb":N_4, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, "box":zone_4, "destination":dest_name}]
# people4 = people_initialization(dom, group4,dt,dmin_people,dmin_walls,rng,itermax,projection_method='cvxopt')

#makes a new destination in the domain and changes the destination to the new one in the dict
# for i in range(N_1):
#         position = people1["xyrv"][i][0:2]
#         circle = Circle(position, radius=1)
#         dom.add_shape(circle,outline_color=[0,0,i+100],fill_color=[0,0,i+100])
#         dest = Destination(name='stationary '+str(i), colors=[[0,0,i+100]])
#         dom.add_destination(dest)
#         people2["destinations"][i] = "stationary "+str(i)

for i in range(N_stationary):
        position = people2["xyrv"][i][0:2]
        circle = Circle(position, radius=1)
        dom.add_shape(circle,outline_color=[0,0,i+100],fill_color=[0,0,i+100])
        dest = Destination(name='stationary '+str(i), colors=[[0,0,i+100]])
        dom.add_destination(dest)
        people2["destinations"][i] = "stationary "+str(i)

# for i in range(N_3):
#         position = people3["xyrv"][i][0:2]
#         circle = Circle(position, radius=1)
#         dom.add_shape(circle,outline_color=[0,0,i+100],fill_color=[0,0,i+100])
#         dest = Destination(name='stationary '+str(i), colors=[[0,0,i+100]])
#         dom.add_destination(dest)
#         people2["destinations"][i] = "stationary "+str(i)

# for i in range(N_4):
#         position = people4["xyrv"][i][0:2]
#         circle = Circle(position, radius=1)
#         dom.add_shape(circle,outline_color=[0,0,i+100],fill_color=[0,0,i+100])
#         dest = Destination(name='stationary '+str(i), colors=[[0,0,i+100]])
#         dom.add_destination(dest)
#         people2["destinations"][i] = "stationary "+str(i)        

#merge multiple dicts with same key values
all_people = {}
for k,v in people.items():
    try:
        all_people[k] = np.concatenate((people[k],people2[k]),axis = 0)
    except:
        all_people[k] = people[k]
people = all_people

safety = np.zeros((ntot,Tf))


# main calculating loop
#plot_people(0,dom,people,contacts,colors)
while(t<Tf):
    print("\n===> Time = "+str(t))
    print("===> Compute desired velocity for domain ",name)
    I, J, Vd = dom.people_desired_velocity(people["xyrv"], people["destinations"])
    people["Vd"] = Vd
    people["I"] = I
    people["J"] = J

    ## Safety stuff'
    t_int = int(t)
    for i in range(ntot):
        if people["xyrv"][i,1] >= dzy_up:
            safety[i,t_int] = ((people["xyrv"][i,1]-dzy_up)/(53.9-dzy_up))**2
        elif people["xyrv"][i,1] <= dzy_down:
            safety[i,t_int] = ((-people["xyrv"][i,1]+dzy_down)/(dzy_down))**2
        else:
            safety[i,t_int] = 0


    print("===> Compute social forces for domain ",name)  
    contacts = compute_contacts(dom, people["xyrv"], dmax)
    print("     Number of contacts: ",contacts.shape[0])
    Forces = compute_forces_dz(F, Fwall, people["xyrv"], contacts, people["Uold"], Vd, lambda_, delta, kappa, eta,Fdz,dzy_down,dzy_up,awr)            
    nn = people["xyrv"].shape[0]
    people["U"] = dt*(Vd[:nn,:]-people["Uold"][:nn,:])/tau + people["Uold"][:nn,:] + dt*Forces[:nn,:]/mass

    people, sensors = move_people(t, dt,people,sensors)
    #people = people_update_destination(people["xyrv"],domains = {"dom"},dom.pixel_size)

    people["Uold"] = people["U"]

    if(draw):
        colors =  people["xyrv"][:,2]
                ## coloring people according to their destinations
                # colors = np.zeros(all_people[name]["xyrv"].shape[0])
                # for i,dest_name ,vin enumerate(all_people[name]["destinations"]):
                #     ind = np.where(all_people[name]["destinations"]==dest_name)[0]
                #     colors[ind]=i
        plot_people(20, dom, people, contacts,
                            colors, time=t,
                            plot_people=True, plot_contacts=False,
                            plot_paths=True, plot_velocities=False,
                            plot_desired_velocities=False, plot_sensors=True, sensors= sensors,savefig=False,
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

plot_sensors(100, sensors, t, savefig=False)
plt.pause(0.01)



safety_total = np.sum(safety,axis=0) 

plt.figure(23489237)
plt.plot(safety_total)


plt.ioff()
plt.show()
sys.exit()






