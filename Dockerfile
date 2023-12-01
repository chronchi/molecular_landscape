# we follow the instructions by bruno rodrigues as described in the
# following chapters of his book:
# https://rap4mads.eu/09-docker.html
# https://rap4mads.eu/10-github_actions.html


# Here we use the chronchi/ember image to run the pipeline. Note that if 
# we start adding more packages the ember image needs to be updated to account
# for that. This Dockerfile here is created for the paper version

FROM chronchi/ember:v1

COPY scripts/ /home/rstudio/ember/scripts

CMD cd /home/rstudio/ember/scripts
CMD bash quarto_render.sh
