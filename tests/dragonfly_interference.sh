#!/bin/bash

## Need to update this path
SST_SRC=$HOME/sst/src

CONFIG_DIR=$PWD/config

echo "SST source code path: " ${SST_SRC}
echo "Reading config files from " ${CONFIG_DIR}

export PYTHONPATH=${CONFIG_DIR}
export PYTHONPATH=${SST_SRC}/sst-element/src/sst/elements/ember/test:${PYTHONPATH}

mpirun -np 32 sst --model-options="--loadFile=${CONFIG_DIR}/df_1056_auto.load --platform=platform_df1056 --topo=dragonfly --algoRting=q-adaptive --stats_startat=0us --numVNs=1 --smallCollectiveVN=-1 --smallCollectiveSize=-1 --qlr1=0.2 --qlr2=0.04 --explore=0.001 --qthld1=0.2 --qthld2=0.35 --qtable_path= --io_level=2 --io_thld=30000000" ${CONFIG_DIR}/emberLoad_df1056.py