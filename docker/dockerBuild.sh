#!/usr/bin/env bash 
docker build --progress=plain -t kjdurbin/$(basename ${PWD}) . 
#&& docker push kjdurbin/$(basename ${PWD})

# For debugging can try:
# docker build --no-cache --progress=plain -t kjdurbin/hic_breakfinder .