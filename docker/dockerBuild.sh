#!/usr/bin/env bash 
docker build -t kjdurbin/$(basename ${PWD}) . && docker push kjdurbin/$(basename ${PWD})
