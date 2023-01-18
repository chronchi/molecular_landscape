#!/usr/bin/env bash

port_image=8787
port_exposed=8000

docker run \
    -p $port_exposed:$port_image \
    --name molecular_landscape \
    -e PASSWORD=password \
    molecular_landscape
    # -e USERID=$(id -u) \
    # -e GROUPID=$(id -g) \
    # molecular_landscape
