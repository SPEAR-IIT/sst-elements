
import sst
from sst.merlin.base import *
from sst.merlin.endpoint import *
from sst.merlin.topology import *
from sst.merlin.interface import *
from sst.merlin.router import *

platdef = PlatformDefinition("platform_dragon")
PlatformDefinition.registerPlatformDefinition(platdef)


##------------------------------------------------

routing_alogrithm = 'q-adaptive' # VALn ugal-g ugal-n par minimal q-adaptive
num_msgs = 5000
msg_size = '128B'
msg_interval = '64ns'
pattern = 'Uniform'  # Uniform, Tornado, Stencil_3D, FFT3D_all2all, RandomNeighbors

dynamicload = False


##------------------------------------------------

#1056-node system
hpr = 4
rpg = 8
inter_link = 1
num_group = 33
size3dx = 4
size3dy = 8
size3dz = 33

#2550-node system:
# hpr = 5
# rpg = 10
# inter_link = 1
# num_group = 51
# size3dx = 5
# size3dy = 10
# size3dz = 51

platdef.addParamSet("topology",{
    "hosts_per_router" : hpr,
    "routers_per_group" : rpg,
    "intergroup_links" : inter_link,
    "num_groups" : num_group,
})

platdef.addParamSet("topology",{
    "algorithm" : routing_alogrithm,
    'link_lat_host' : '10ns',
    'link_lat_local' : '30ns',
    'link_lat_global' : '300ns',
})


platdef.addParamSet("topology",{
    "adaptive_threshold": 2,
})

platdef.addParamSet("topology",{
    "learning_rate": 0.2,
    "learning_rate2": 0.04,
    "epsilon": 0.001,
    "q_threshold1": 0.2,
    "q_threshold2": 0.35,

    "save_qtable": 'yes',
    "save_qtable_time": 1000, ##us

})

platdef.addClassType("topology","sst.merlin.topology.topoDragonFly")

platdef.addParamSet("router",{
    "link_bw" : '4GB/s',
    "link_bw:host" : '4GB/s',
    "xbar_bw" : '40GB/s',
    "flit_size" : '128B',
    "input_buf_size" : '2560B',
    "output_buf_size" : '2560B',
    "input_latency" : "10ns",
    "output_latency" : "10ns",
    "input_buf_size:host" : '2560B',
    "output_buf_size:host" : '2560B',
    "num_vns" : 1,
    "xbar_arb" : "merlin.xbar_arb_lru",
})

platdef.addClassType("router","sst.merlin.base.hr_router")

platdef.addParamSet("network_interface",{
    "link_bw" : '4GB/s',
    "input_buf_size" : '2560B',
    "output_buf_size" : '2560B',
})

#platdef.addClassType("network_interface","sst.merlin.base.ReorderLinkControl")
platdef.addClassType("network_interface","sst.merlin.interface.LinkControl")


PlatformDefinition.setCurrentPlatform("platform_dragon")

# Allocate the system. This will create the topology since it's
# set up in the platform file
system = System()

### set up the endpoint

#----------------------------------------

syssize = system.topology.getNumNodes()

ep = TrafficGenJob(0, syssize)
ep.packets_to_send = num_msgs
ep.packet_size = msg_size
ep.delay_between_packets = msg_interval
ep.message_rate = '1GHz'

ep.extargs["PacketDest:pattern"] = pattern

if dynamicload:
    ep.extargs["PacketDelay:pattern"] = "Step"
    ep.extargs["PacketDelay:packet_delay_list"] = ['160ns', '80ns']
    ep.extargs["PacketDelay:packet_num_list"] = [1000,1000]

if pattern == 'Tornado':
    ep.topology = 'dragonfly'
    ep.extargs["dragonfly:hosts_per_router"] = hpr
    ep.extargs["dragonfly:routers_per_group"] = rpg
    ep.extargs["dragonfly:num_groups"] = num_group
    ep.extargs["Tornado:shift"] = 1

if pattern == 'Stencil_3D':
    ep.extargs["PacketDest:Stencil_3D:3DSize"] = "{} {} {}".format(size3dx, size3dy, size3dz)  

if pattern == 'FFT3D_all2all':
    ep.extargs["PacketDest:FFT3D_all2all:3DSize"] = "{} {} {}".format(size3dx, size3dy, size3dz)  

if pattern == 'RandomNeighbors':
    ep.extargs["PacketDest:RandomNeighbors:range_min"] = 6
    ep.extargs["PacketDest:RandomNeighbors:range_max"] = 21

system.allocateNodes(ep,"linear") 

system.build()

sst.setStatisticLoadLevel(9)
    
sst.setStatisticOutput("sst.statOutputCSV")
sst.setStatisticOutputOptions({
    "filepath" : 'stats.csv',
    "separator" : ", "
})

sst.enableAllStatisticsForAllComponents()
