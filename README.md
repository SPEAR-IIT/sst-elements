![SST](http://sst-simulator.org/img/sst-logo-small.png)

# Structural Simulation Toolkit (SST)

#### This is a fork from the official [SST Github page](https://github.com/sstsimulator/sst-elements) with implementations maintained by [SPEAR group](http://www.cs.iit.edu/~zlan/) at IIT.  
Visit [sst-simulator.org](http://sst-simulator.org) to learn more about SST.

---
## Study of Workload Interference with Intelligent Routing on Dragonfly

This branch contains the SST enhancement for Dragonfly network interference study

### Installation
Follow [installation guide](http://sst-simulator.org/SSTPages/SSTBuildAndInstall_11dot1dot0_SeriesDetailedBuildInstructions/) to install SST.

### Run a test

```bash
cd tests/
./dragonfly_interference.sh
```

#### Note
**tests/dragonfly_interference.sh** will launch a CosmoFlow-UR pairwise simulation using Q-adaptive routing.

The **SST_SRC** variable needs to be updated to point the directory containing sst-element source code

---




##### [LICENSE](https://github.com/sstsimulator/sst-elements/blob/devel/LICENSE)
