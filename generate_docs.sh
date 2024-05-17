IMAGE_NAME="run_analysis"

# first build the image that will generate the docs at the same time
sudo docker build --cpuset-cpus=0-30 -f Dockerfile.report -t ${IMAGE_NAME} .

# now get the docs
docker create --name dummy ${IMAGE_NAME}
docker cp dummy:/home/rstudio/ember/docs/ ./docs
docker rm -f dummy

# remove the image since it was created only for the docs
sudo docker rmi ${IMAGE_NAME}
