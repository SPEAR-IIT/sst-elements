#!/usr/bin/env python
#
# Copyright 2009-2021 NTESS. Under the terms
# of Contract DE-NA0003525 with NTESS, the U.S.
# Government retains certain rights in this software.
#
# Copyright (c) 2009-2021, NTESS
# All rights reserved.
#
# This file is part of the SST software package. For license
# information, see the LICENSE file in the top level directory of the
# distribution.

import sst
from sst.merlin.base import *
from sst.merlin.endpoint import *
from sst.merlin.interface import *
from sst.merlin.topology import *

from sst.ember import *

if __name__ == "__main__":

    PlatformDefinition.setCurrentPlatform("firefly-defaults")

    ### Setup the topology
    topo = topoDragonFly()

    topo.addParams({
        "hosts_per_router" : 3,
        "routers_per_group" : 6,
        "intergroup_links" : 1,
        "num_groups" : 19,
        "algorithm": ["ugal-4vc"], # "minimal","ugal-4vc" "q-routing1"
        'link_lat_host' : '10ns',
        'link_lat_local' : '30ns',
        'link_lat_global' : '300ns',
    })

    # Set up the routers
    router = hr_router()
    router.addParams({
        "link_bw" : "4GB/s",
        "link_bw:host" : "8GB/s",
        "xbar_bw" : "10GB/s",
        "flit_size" : "8B",
        "input_latency" : "10ns",
        "output_latency" : "10ns",
        "input_buf_size" : "4kB",
        "output_buf_size" : "4kB",
        # "input_buf_size:host" : "4kB",
        # "output_buf_size:host" : "4kB",
        "num_vns" : 1,
        # Optionally set up the QOS
        # "qos_settings" : [50,50],
        # Set up the arbitration type for the routers
        "xbar_arb" : "merlin.xbar_arb_lru",
        # "xbar_arb" : "merlin.xbar_arb_age",
        # "xbar_arb" : "merlin.xbar_arb_rr",
    })
    topo.router = router

    ### set up the endpoint
    networkif = ReorderLinkControl()
    networkif.link_bw = "4GB/s"
    networkif.input_buf_size = "1kB"
    networkif.output_buf_size = "1kB"

    # Set up VN remapping
    #networkif.vn_remap = [0]
    
    # ep = EmberMPIJob(0,topo.getNumNodes())
    ep = EmberMPIJob(0,114)
    ep.network_interface = networkif
    ep.addMotif("Init")
    # ep.addMotif("Sweep3D pex=2 pey=57 nx=16 ny=16 nz=10 kba=10 fields_per_cell=10 iterations=1 computetime=1")
    ep.addMotif("Allreduce iterations=3")
    ep.addMotif("Fini")
    ep.nic.nic2host_lat= "100ns"
        
    ep2 = EmberMPIJob(1,114)
    ep2.network_interface = networkif
    ep2.addMotif("Init")
    # ep2.addMotif("Sweep3D pex=2 pey=57 nx=16 ny=16 nz=10 kba=10 fields_per_cell=10 iterations=1 computetime=1")
    ep2.addMotif("Allreduce")
    ep2.addMotif("Fini")
    ep2.nic.nic2host_lat= "100ns"


    ep3 = EmberMPIJob(2,114)
    ep3.network_interface = networkif
    ep3.addMotif("Init")
    # ep3.addMotif("Sweep3D pex=2 pey=57 nx=16 ny=16 nz=10 kba=10 fields_per_cell=10 iterations=1 computetime=1")
    ep3.addMotif("Allreduce")
    ep3.addMotif("Fini")
    ep3.nic.nic2host_lat= "100ns"


    system = System()
    system.setTopology(topo)

    system.allocateNodes(ep,"linear")
    system.allocateNodes(ep2,"linear")
    # system.allocateNodes(ep3,"linear")

    # ep3 = EmberMPIJob(2,114)
    # ep3.network_interface = networkif
    # ep3.addMotif("Init")
    # ep3.addMotif("Sweep3D pex=2 pey=57 nx=16 ny=16 nz=10 kba=10 fields_per_cell=10 iterations=1 computetime=1")
    # ep3.addMotif("Allreduce")
    # ep3.addMotif("Fini")
    # ep3.nic.nic2host_lat= "100ns"

    # system.allocateNodes(ep3,"linear")
    
    system.build()

    # sst.setStatisticLoadLevel(9)

    # sst.setStatisticOutput("sst.statOutputCSV");
    # sst.setStatisticOutputOptions({
    #     "filepath" : "stats.csv",
    #     "separator" : ", "
    # })

