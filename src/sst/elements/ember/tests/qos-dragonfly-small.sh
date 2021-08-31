LOADFILE="qos_small.load"
PLATFORM="platform_df"

export PYTHONPATH="../test"

# df342="3:6:1:19"
# --shape=$df342 \

sst \
--model-options=" \
--loadFile=$LOADFILE \
--platform=$PLATFORM \
--topo=dragonfly \
" \
../test/emberLoad_df.py

