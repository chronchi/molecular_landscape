#!/bin/bash

base_name="ember"

# first build the image that will generate the docs
sudo docker build --cpuset-cpus=0-30 -f Dockerfile.data -t ${base_name} .

image_name="chronchi/ember" 
tags=$(curl -sS "https://registry.hub.docker.com/v2/repositories/chronchi/ember/tags/" | jq -r '.results[] | .name' | grep -oE '[0-9]+$' | sort -n)

last_version=$(echo "$tags" | tail -n1)
version_number=$(echo "$last_version" | grep -oE '[0-9]+')

if [ -n "$version_number" ]; then
    new_version=$((version_number + 1))
    new_tag="${image_name}:v${new_version}"

    echo "The new version tag is $new_tag"

    sudo docker tag ${base_name} ${image_name} 
    sudo docker push ${image_name}

    # Optionally, tag the Docker image with the new version
    sudo docker tag "${image_name}:latest" "$new_tag"
    sudo docker push "$new_tag"
else

    echo "No version number found in the latest version of ${image_name}, pushing v1"
    new_tag="${image_name}:v1"

    echo "The new version tag is $new_tag"

    sudo docker tag ${base_name} ${image_name} 
    sudo docker push ${image_name}

    # Optionally, tag the Docker image with the new version
    sudo docker tag "${image_name}:latest" "$new_tag"
    sudo docker push "$new_tag"
fi
