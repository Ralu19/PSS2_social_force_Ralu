# Authors:
#     Sylvain Faure <sylvain.faure@universite-paris-saclay.fr>
#     Bertrand Maury <bertrand.maury@universite-paris-saclay.fr>
#
#      cromosim/examples/micro/social/micro_social.py
#      python micro_social.py --json input.json
#
# License: GPL


import sys, os
script_dir = os.path.dirname( __file__ )
mymodule_dir = os.path.join( script_dir, '..', '..', '..',)
sys.path.append( mymodule_dir )

from cromosim import *
from cromosim.micro import *
from optparse import OptionParser
import json

plt.ion()

"""
    python micro_social_modified.py --json input.json
"""
parser = OptionParser(usage="usage: %prog [options] filename",
    version="%prog 1.0")
parser.add_option('--json',dest="jsonfilename",default="input.json",
    type="string",
                  action="store",help="Input json filename")
opt, remainder = parser.parse_args()
print("===> JSON filename = ",opt.jsonfilename)
with open(opt.jsonfilename) as json_file:
    try:
        input = json.load(json_file)
    except json.JSONDecodeError as msg:
        print(msg)
        print("Failed to load json file ",opt.jsonfilename)
        print("Check its content \
            (https://fr.wikipedia.org/wiki/JavaScript_Object_Notation)")
        sys.exit()

"""
    Get parameters from json file :
    prefix: string
        Folder name to store the results
    with_graphes: bool
        true if all the graphes are shown and saved in png files,
        false otherwise
    seed: integer
        Random seed which can be used to reproduce a random selection
        if >0
    For each domain :
    |    name: string
    |        Domain name
    |    background: string
    |        Image file used as background
    |    px: float
    |        Pixel size in meters (also called space step)
    |    width: integer
    |        Domain width (equal to the width of the background image)
    |    height: integer
    |        Domain height (equal to the height of the background image)
    |    wall_colors: list
    |        rgb colors for walls
    |        [ [r,g,b],[r,g,b],... ]
    |    shape_lines: list
    |        Used to define the Matplotlib Polyline shapes,
    |        [
    |          {
    |             "xx": [x0,x1,x2,...],
    |             "yy": [y0,y1,y2,...],
    |             "linewidth": float,
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |          },...
    |        ]
    |    shape_circles: list
    |        Used to define the Matplotlib Circle shapes,
    |        [
    |           {
    |             "center_x": float,
    |             "center_y": float,
    |             "radius": float,
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    shape_ellipses: list
    |        Used to define the Matplotlib Ellipse shapes,
    |        [
    |           {
    |             "center_x": float,
    |             "center_y": float,
    |             "width": float,
    |             "height": float,
    |             "angle_in_degrees_anti-clockwise": float (degre),
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    shape_rectangles: list
    |        Used to define the Matplotlib Rectangle shapes,
    |        [
    |           {
    |             "bottom_left_x": float,
    |             "bottom_left_y": float,
    |             "width": float,
    |             "height": float,
    |             "angle_in_degrees_anti-clockwise": float (degre),
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    shape_polygons: list
    |        Used to define the Matplotlib Polygon shapes,
    |        [
    |           {
    |             "xy": float,
    |             "outline_color": [r,g,b],
    |             "fill_color": [r,g,b]
    |            },...
    |        ]
    |    destinations: list
    |        Used to define the Destination objects,
    |        [
    |           {
    |             "name": string,
    |             "colors": [[r,g,b],...],
    |             "excluded_colors": [[r,g,b],...],
    |             "desired_velocity_from_color": [] or
    |             [
    |                {
    |                   "color": [r,g,b],
    |                   "desired_velocity": [ex,ey]
    |                },...
    |             ],
    |             "velocity_scale": float,
    |             "next_destination": null or string,
    |             "next_domain": null or string,
    |             "next_transit_box": null or [[x0,y0],...,[x3,y3]]
    |            },...
    |        ]
    |--------------------
    For each group of persons, required for the initialization process:
    |    nb:
    |        Number of people in the group
    |    domain:
    |        Name of the domain where people are located
    |    radius_distribution:
    |        Radius distribution used to create people
    |        ["uniform",min,max] or ["normal",mean,sigma]
    |    velocity_distribution:
    |        Velocity distribution used to create people
    |        ["uniform",min,max] or ["normal",mean,sigma]
    |    box:
    |        Boxe to randomly position people at initialization
    |        [ [x0,y0],[x1,y1],...]
    |    destination:
    |        Initial destination for the group
    |--------------------
    For each sensor:
    |    domain:
    |        Name of the domain where the sensor is located
    |    line:
    |        Segment through which incoming and outgoing flows are measured
    |        [x0,y0,x1,y1]
    |--------------------
    Tf: float
        Final time
    dt: float
        Time step
    drawper: integer
        The results will be displayed every "drawper" iterations
    mass: float
        Mass of one person (typically 80 kg)
    tau: float
        (typically 0.5 s)
    F: float
        Coefficient for the repulsion force between individuals
        (typically 2000 N)
    kappa: float
        Stiffness constant to handle overlapping (typically
        120000 kg s^-2)
    delta: float
        To maintain a certain distance from neighbors (typically 0.08 m)
    Fwall: float
        Coefficient for the repulsion force between individual and
        walls (typically 2000 N, like for F)
    lambda: float
        Directional dependence (between 0 and 1 = fully isotropic case)
    eta: float
        Friction coefficient (typically 240000 kg m^-1 s^-1)
    projection_method: string
        Name of the projection method : cvxopt(default),
        mosek(a licence is needed) or uzawa
    dmax: float
        Maximum distance used to detect neighbors
    dmin_people: float
        Minimum distance allowed between individuals
    dmin_walls: float
        Minimum distance allowed between an individual and a wall
    plot_people: boolean
        If true, people are drawn
    plot_contacts: boolean
        If true, active contacts between people are drawn
    plot_desired_velocities: boolean
        If true, people desired velocities are drawn
    plot_velocities: boolean
        If true, people velocities are drawn
    plot_sensors: boolean
        If true, plot sensor lines on people graph and sensor data graph
    plot_paths: boolean
        If true, people paths are drawn
"""
#boundaries of domain:
# x (without the black borders) \in ~(0.6, 296.6) (297.6 for the entire pic)
# y ( --//--) \in ~(1.1, 53.9) (54 for the entire pic)

#stairs: [(11,33),(14.1,17.9)] (escalators down), [(5.5,33),(38.7,41.7)] (escalators up),
#[(6, 33),(20, 36.7)] (normal stairs)

prefix = input["prefix"]
if not os.path.exists(prefix):
    os.makedirs(prefix)
seed = input["seed"]
with_graphes = input["with_graphes"]
json_domains = input["domains"]
#print("===> JSON data used to build the domains : ",json_domains)
json_people_init = input["people_init"]
#print("===> JSON data used to create the groups : ",json_people_init)
json_sensors = input["sensors"]
#print("===> JSON data used to create sensors : ",json_sensors)
nn = 20
N_stationary = input["N_stationary"] #number of stationary people
Tf = input["Tf"]
dt = input["dt"]
drawper = input["drawper"]
mass = input["mass"]
tau = input["tau"]
F = input["F"]

#DZ specific inputs
Fdz = input["Fdz"] 
dzy_up = input["dzy_up"] 
dzy_down = input["dzy_down"]
#box_dz_up = [[1.0,296.0, 51.0, 52.3]]
box_dz_down = [[1.0,296.0, 1.0, 4.5]]

## ZONES 
## we can adjust the boxes how we like
zone_1 = [[1.0, 50.2 , 1.0, 52.3 ]]
zone_2 = [[50.21, 125.0, 1.0, 52.3]]
zone_3 = [[125.01, 216.0 , 1.0, 52.3]]
zone_4 =[[216.01 , 296.0, 1.0, 52.3]]

kappa = input["kappa"]
delta = input["delta"]
Fwall = input["Fwall"]
lambda_ = input["lambda"]
eta = input["eta"]
projection_method = input["projection_method"]
dmax = input["dmax"]
dmin_people = input["dmin_people"]
dmin_walls = input["dmin_walls"]
plot_p = input["plot_people"]
plot_c = input["plot_contacts"]
plot_v = input["plot_velocities"]
plot_vd = input["plot_desired_velocities"]
plot_pa = input["plot_paths"]
plot_s = input["plot_sensors"]
plot_pa = input["plot_paths"]
radius_distribution = ["uniform",0.4,0.6]
velocity_distribution = ["normal",1.2,0.1] # distribution varible
print("===> Final time, Tf = ",Tf)
print("===> Time step, dt = ",dt)
print("===> To draw the results each drawper iterations, \
    drawper = ",drawper)
print("===> Maximal distance to find neighbors, dmax = ",
    dmax,", example : 2*dt")
print("===> ONLY used during initialization ! Minimal distance between \
       persons, dmin_people = ",dmin_people)
print("===> ONLY used during initialization ! Minimal distance between a \
       person and a wall, dmin_walls = ",dmin_walls)

 
## Awareness stuff ##
#The awareness values should be between 0 and 1
mean_awr = 0.5 # mean of awareness
stdev_awr = 0.2 # stdev of awareness
awr = np.random.normal(mean_awr, stdev_awr, nn) # awareness[i] is the awareness value of person i
for i in range(awr.shape[0]): #Just in case you get unlucky
    if awr[i] < 0:
        awr[i] = 0
    elif awr[i] > 1:
        awr[i] = 1
#awr = np.ones(nn) #for testing 
  
## Safety stuff


"""
    Build the Domain objects
"""
domains = {}
for i,jdom in enumerate(json_domains):
    jname = jdom["name"]
    print("===> Build domain number ",i," : ",jname)
    jbg = jdom["background"]
    jpx = jdom["px"]
    jwidth = jdom["width"]
    jheight = jdom["height"]
    jwall_colors = jdom["wall_colors"]
    if (jbg==""):
        dom = Domain(name=jname, pixel_size=jpx, width=jwidth,
                     height=jheight, wall_colors=jwall_colors)
    else:
        dom = Domain(name=jname, background=jbg, pixel_size=jpx,
                     wall_colors=jwall_colors)
    ## To add lines : Line2D(xdata, ydata, linewidth)
    for sl in jdom["shape_lines"]:
        line = Line2D(sl["xx"],sl["yy"],linewidth=sl["linewidth"])
        dom.add_shape(line,outline_color=sl["outline_color"],
                      fill_color=sl["fill_color"])
    ## To add circles : Circle( (center_x,center_y), radius )
    for sc in jdom["shape_circles"]:
        circle = Circle( (sc["center_x"], sc["center_y"]), sc["radius"] )
        dom.add_shape(circle,outline_color=sc["outline_color"],
                      fill_color=sc["fill_color"])
    ##To add the circles for the stationary people 
    ##stationary people = obstacles; hence the color black 
    ##(which is the wall color)
    for c in range(N_stationary):
        ## boundaries of kiosk: [124,215, 17.3,38.3]
        ## boundaries of some other thingy at the end of the platform: [286,297, 16.5,38.1]  
        ## y's of waiting places in zone_1: [(1,9),(46.5,54)]
        dummy = random.choices([1,2,3,4], [1,5,20,10])
        if dummy[0] == 1:
            x_rnd = random.uniform(zone_1[0][0],zone_1[0][1]) 
            y_rnd = random_from_intervals(((1,9), (46.5,54)))
        elif dummy[0] == 2:
            x_rnd = random.uniform(zone_2[0][0],zone_2[0][1]) 
            y_rnd = random.uniform(zone_2[0][2],zone_2[0][3])
        elif dummy[0] == 3:
            x_rnd = random.uniform(zone_3[0][0],zone_3[0][1]) 
            y_rnd = random_from_intervals(((1,17.3),(38.3,54)))
            #move ppl closer to walls depending on the mean awareness
            if y_rnd < 7:
                chance = random.uniform(0,1)
                if chance > mean_awr:
                    y_rnd = y_rnd + 10 
            elif y_rnd > 48.3:
                chance = random.uniform(0,1)
                if chance > mean_awr:
                    y_rnd = y_rnd - 10      
        else:    
            x_rnd = random.uniform(zone_4[0][0],zone_4[0][1])
            if (x_rnd>286.0 and x_rnd<297):
                y_rnd = random_from_intervals(((1,17.3),(38.3,54)))
            else:
                y_rnd = random.uniform(1,54)           
        circles = Circle((x_rnd, y_rnd), 0.35) #random location of dots
        dom.add_shape(circles, outline_color=[0,0,0],fill_color=[0,0,0])

    ## To add ellipses : Ellipse( (center_x,center_y), width, height,
    ##                            angle_in_degrees_anti-clockwise )
    for se in jdom["shape_ellipses"]:
        ellipse = Ellipse( (se["center_x"], se["center_y"]),
                            se["width"], se["height"],
                            se["angle_in_degrees_anti-clockwise"])
        dom.add_shape(ellipse,outline_color=se["outline_color"],
                      fill_color=se["fill_color"])
    ## To add rectangles : Rectangle( (bottom_left_x,bottom_left_y),
    ##                       width, height, angle_in_degrees_anti-clockwise )
    for sr in jdom["shape_rectangles"]:
        rectangle = Rectangle( (sr["bottom_left_x"],sr["bottom_left_y"]),
                               sr["width"], sr["height"],
                               sr["angle_in_degrees_anti-clockwise"])
        dom.add_shape(rectangle,outline_color=sr["outline_color"],
                      fill_color=sr["fill_color"])
    ## To add polygons : Polygon( [[x0,y0],[x1,y1],...] )
    for spo in jdom["shape_polygons"]:
        polygon = Polygon(spo["xy"])
        dom.add_shape(polygon,outline_color=spo["outline_color"],
                      fill_color=spo["fill_color"])
    ## To build the domain : background + shapes
    dom.build_domain()
    
    ## To add all the available destinations
    for j,dd in enumerate(jdom["destinations"]):
        desired_velocity_from_color=[]
        for gg in dd["desired_velocity_from_color"]:
            desired_velocity_from_color.append(
                np.concatenate((gg["color"],gg["desired_velocity"])))
        dest = Destination(name=dd["name"],colors=dd["colors"],
        excluded_colors=dd["excluded_colors"],
        desired_velocity_from_color=desired_velocity_from_color,
        velocity_scale=dd["velocity_scale"],
        next_destination=dd["next_destination"],
        next_domain=dd["next_domain"],
        next_transit_box=dd["next_transit_box"])
        print("===> Destination : ",dest)
        dom.add_destination(dest)
        #if (with_graphes):
        #    dom.plot_desired_velocity(dd["name"],id=100*i+10+j,step=20)
    

    print("===> Domain : ",dom)
    #if (with_graphes):
    #    dom.plot(id=100*i)
    #    dom.plot_wall_dist(id=100*i+1,step=20)

    domains[dom.name] = dom

print("===> All domains = ",domains)


"""
    To create the sensors to measure the pedestrian flows
"""

all_sensors = {}
for domain_name in domains:
    all_sensors[domain_name] = []
for s in json_sensors:
    s["id"] = []
    s["times"] = []
    s["xy"] = []
    s["dir"] = []
    all_sensors[s["domain"]].append(s)
    print("===> All sensors = ",all_sensors)

"""
    Initialization
"""

## Current time
t = 0.0
counter = 0



# plt.show()
# group2 = [{"nb":5, "radius_distribution":radius_distribution, "velocity_distribution":velocity_distribution, \
#     "box":[100,120,10,20], "destination":"dest1"}]
# people2 = people_initialization(dom, group2,dt,dmin_people=dmin_people,dmin_walls=dmin_walls,seed=seed,\
#         itermax=10,projection_method=projection_method, verbose=True)

# #makes a new destination in the domain and changes the destination to the new one in the dict
# for i in range(5):
#     position = people2["xyrv"][i][0:2]
#     circle = Circle(position, radius=1)
#     dom.add_shape(circle,outline_color=[0,0,i+100],fill_color=[0,0,i+100])
#     desti = Destination(name='stationary '+str(i), colors=[[0,0,i+100]])
#     dom.add_destination(desti)
#     people2["destinations"][i] = "stationary "+str(i)

## Initialize people
all_people = {}
for i,peopledom in enumerate(json_people_init):
    dom = domains[peopledom["domain"]]
    groups = peopledom["groups"]
    print("===> Group number ",i,", domain = ",peopledom["domain"])
    people = people_initialization(dom, groups, dt,
        dmin_people=dmin_people, dmin_walls=dmin_walls, seed=seed,
        itermax=10, projection_method=projection_method, verbose=True)
    I, J, Vd = dom.people_desired_velocity(people["xyrv"],
        people["destinations"])
    people["Vd"] = Vd
    for ip,pid in enumerate(people["id"]):
        people["paths"][pid] = people["xyrv"][ip,:2]
    contacts = None
    if (with_graphes):
        colors = people["xyrv"][:,2]
        plot_people(100*i+20, dom, people, contacts, colors, time=t,
                    plot_people=plot_p, plot_contacts=plot_c,
                    plot_velocities=plot_v, plot_desired_velocities=plot_vd,
                    plot_sensors=plot_s, sensors=all_sensors[dom.name],
                    savefig=True, filename=prefix+dom.name+'_fig_'+ \
                    str(counter).zfill(6)+'.png')
    all_people[peopledom["domain"]] = people
#print("===> All people = ",all_people)

"""
    Main loop
"""

cc = 0
draw = True
safety = np.zeros((nn,Tf))

while (t<Tf):

    print("\n===> Time = "+str(t))

    ## Compute people desired velocity
    for idom,name in enumerate(domains):
        print("===> Compute desired velocity for domain ",name)
        dom = domains[name]
        people = all_people[name]
        I, J, Vd = dom.people_desired_velocity(people["xyrv"],
            people["destinations"])
        people["Vd"] = Vd
        people["I"] = I
        people["J"] = J

    ## Look at if there are people in the transit boxes
    print("===> Find people who have to be duplicated")
    virtual_people = find_duplicate_people(all_people, domains)
    #print("     virtual_people : ",virtual_people)

    ## Social forces
    for idom,name in enumerate(domains):
        print("===> Compute social forces for domain ",name)
        dom = domains[name]
        people = all_people[name]

        try:
            xyrv = np.concatenate((people["xyrv"],
                virtual_people[name]["xyrv"]))
            Vd = np.concatenate((people["Vd"],
                virtual_people[name]["Vd"]))
            Uold = np.concatenate((people["Uold"],
                virtual_people[name]["Uold"]))
        except:
            xyrv = people["xyrv"]
            Vd = people["Vd"]
            Uold = people["Uold"]

        if (xyrv.shape[0]>0):

            if (np.unique(xyrv, axis=0).shape[0] != xyrv.shape[0]):
                print("===> ERROR : There are two identical lines in the")
                print("             array xyrv used to determine the \
                    contacts between")
                print("             individuals and this is not normal.")
                sys.exit()

            contacts = compute_contacts(dom, xyrv, dmax)
            print("     Number of contacts: ",contacts.shape[0])
            Forces = compute_forces_dz(F, Fwall, xyrv, contacts, Uold, Vd,
                    lambda_, delta, kappa, eta, Fdz, dzy_down, dzy_up, awr)
            nn = people["xyrv"].shape[0]
            all_people[name]["U"] = dt*(Vd[:nn,:]-Uold[:nn,:])/tau + \
                          Uold[:nn,:] + \
                          dt*Forces[:nn,:]/mass
            ## only for the plot of virtual people :
            virtual_people[name]["U"] = dt*(Vd[nn:,:]-Uold[nn:,:])/tau + \
                          Uold[nn:,:] + \
                          dt*Forces[nn:,:]/mass


            all_people[name], all_sensors[name] = move_people(t, dt,
                                           all_people[name],
                                           all_sensors[name])

        if (draw and with_graphes):
            ## coloring people according to their radius
            colors =  all_people[name]["xyrv"][:,2]
            ## coloring people according to their destinations
            # colors = np.zeros(all_people[name]["xyrv"].shape[0])
            # for i,dest_name in enumerate(all_people[name]["destinations"]):
            #     ind = np.where(all_people[name]["destinations"]==dest_name)[0]
            #     colors[ind]=i
            plot_people(100*idom+20, dom, all_people[name], contacts,
                        colors, virtual_people=virtual_people[name], time=t,
                        plot_people=plot_p, plot_contacts=plot_c,
                        plot_paths=plot_pa, plot_velocities=plot_v,
                        plot_desired_velocities=plot_vd, plot_sensors=plot_s,
                        sensors=all_sensors[dom.name], savefig=True,
                        filename=prefix+dom.name+'_fig_'
                        + str(counter).zfill(6)+'.png')
            plt.pause(0.01)

    ## Safety stuff'
    t_int = int(t)
    for i in range(nn):
        if people["xyrv"][i,1] >= dzy_up:
            safety[i,t_int] = ((people["xyrv"][i,1]-dzy_up)/(53.9-dzy_up))**2
        elif people["xyrv"][i,1] <= dzy_down:
            safety[i,t_int] = ((-people["xyrv"][i,1]+dzy_down)/(dzy_down))**2
        else:
            safety[i,t_int] = 0

    ## Update people destinations
    all_people = people_update_destination(all_people,domains,dom.pixel_size)

    ## Update previous velocities
    for idom,name in enumerate(domains):
        all_people[name]["Uold"] = all_people[name]["U"]

    ## Print the number of persons for each domain
    for idom,name in enumerate(domains):
        print("===> Domain ",name," nb of persons = ",
            all_people[name]["xyrv"].shape[0])

    t += dt
    cc += 1
    counter += 1
    if (cc>=drawper):
        draw = True
        cc = 0
    else:
        draw = False


for idom,domain_name in enumerate(all_sensors):
    print("===> Plot sensors of domain ",domain_name)
    plot_sensors(100*idom+40, all_sensors[domain_name], t, savefig=True,
                filename=prefix+'sensor_'+str(i)+'_'+str(counter)+'.png')
    plt.pause(0.01)

plt.ioff()
plt.show()
sys.exit()
