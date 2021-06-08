![SST](http://sst-simulator.org/img/sst-logo-small.png)

# Structural Simulation Toolkit (SST)

#### This is a fork from the official [SST Github page](https://github.com/sstsimulator/sst-elements) with implementations maintained by [SPEAR group](http://www.cs.iit.edu/~zlan/) at IIT.  
Visit [sst-simulator.org](http://sst-simulator.org) to learn more about SST.

---
## Q-adaptive routing

The current master branch has the implementation of Q-adatpive routing in Merlin Dragonfly topology.

Note: if you use this SST implementation, please cite the following papers:

Yao Kang, Xin Wang, and Zhiling Lan. "Q-adaptive: A Multi-Agent Reinforcement Learning Based Routing on Dragonfly Network". In Proceedings of the 30th International Symposium on High-Performance Parallel and Distributed Computing (HPDC ’21).

### Installation
Clone [sst-core](https://github.com/SPEAR-IIT/sst-core) and [sst-elements](https://github.com/SPEAR-IIT/sst-elements) from this repository.
Use other SST versions may cause compatible issue.  

Follow [installation guide](http://sst-simulator.org/SSTPages/SSTBuildAndInstall10dot1dot0SeriesDetailedBuildInstructions/) to install SST.

### Run a test

```bash
mpirun -np 4 sst sst-elements/tests/dragonfly_q-adaptive.py 
```
---




##### [LICENSE](https://github.com/sstsimulator/sst-elements/blob/devel/LICENSE)
