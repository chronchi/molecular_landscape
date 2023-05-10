## A tour on the molecular landscape

# basic packages to get a 3D visualization
library(ggplot2)
# remotes::install_github("tylermorganwall/rayrender")
library(rayrender)
library(rayshader)

# the packages for the roller coaster
library(rayimage)
library(raster)

ggdiamonds <- ggplot2::ggplot(diamonds, aes(x, depth)) +
    ggplot2::stat_density_2d(
        ggplot2::aes(fill = stat(nlevel)), 
        geom = "polygon", n = 200, bins = 100,contour = TRUE
    ) +
    ggplot2::facet_wrap(clarity~.) +
    ggplot2::scale_fill_viridis_c(option = "A")

ggdiamonds

plot_gg(
    ggdiamonds,
    width=7,
    height=7,
    scale=250,
    windowsize=c(1100,700),
    raytrace=FALSE, 
    zoom = 0.55, 
    phi = 30, 
    triangulate = TRUE, 
    max_error = 0,
    multicore = TRUE
)
render_snapshot()

# this is the command to save the object
save_obj("scripts/visual_ray/ggplot.obj")
rgl::rgl.close()

# we now proceed to perform the example shown in the
# webpage
disk(
    radius=1000,
    y=-1, 
    material = diffuse(
        checkerperiod = 6,
        checkercolor = "#0d401b", 
        color = "#496651"
    )
) %>% 
    add_object(
        obj_model(
            "scripts/visual_ray/ggplot.obj",
            y=-0.02, 
            texture=TRUE, 
            scale_obj = 1/100
        )
    ) %>%
    add_object(
        sphere(
            y = 30,
            z = 10,
            radius = 5,
            material = light(intensity=40)
        )
    ) %>% 
    render_scene(
        lookfrom = c(20,20,20),
        fov=0,
        ortho_dimensions = c(30,30), 
        width=800,
        height=800
    )
