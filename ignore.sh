#!/bin/bash

target=(cacheTracer cassini CramSim GNA kingsley shogun thornhill ariel Messier miranda Opal prospero Samba scheduler serrano VaultSimC zodiac vanadis firefly hermes ember memHierarchy simpleElementExample simpleSimulation merlin firefly)

# firefly hermes merlin simpleElementExample simpleSimulation ember memHierarchy

for i in "${target[@]}"
do
   touch src/sst/elements/$i/.ignore
done
