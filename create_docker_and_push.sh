#!/bin/bash

image_name="chronchi/ember"  # Replace with your Docker image name
tags=$(curl -sS "https://registry.hub.docker.com/v2/repositories/chronchi/ember/tags/" | jq -r '.results[] | .name' | grep -oE '[0-9]+$' | sort -n)

last_version=$(echo "$tags" | tail -n1)
version_number=$(echo "$last_version" | grep -oE '[0-9]+')

if [ -n "$version_number" ]; then
    new_version=$((version_number + 1))
    new_tag="${image_name}:v${new_version}"

    echo "The new version tag is $new_tag"

    docker tag ember ${image_name} 
    docker push ${image_name}

    # Optionally, tag the Docker image with the new version
    docker tag "${image_name}:latest" "$new_tag"
    docker push "$new_tag"
else
    echo "No version number found in the latest version of $image_name"
fi
