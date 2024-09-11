include("intro.jl")

rad_gf = npzread("data/GF_rad_echam-001.npy")
sst_rcp85 = npzread("data/rcp85_tos.npy")
nt, ny, nx = size(sst_rcp85)

hasnan(sst_rcp85)
hasnan(rad_gf)

continent_mask = isnan.(sum(sst_rcp85, dims=1)[1, :, :])
heatmap(continent_mask')

replace!(sst_rcp85, NaN => 0)
replace!(rad_gf, NaN => 0)

hasnan(sst_rcp85)
hasnan(rad_gf)

sst_rcp85_base = mean(sst_rcp85[1:21, :, :], dims=1)[1, :, :]
heatmap(sst_rcp85_base')

sst_rcp85_anom = zeros(size(sst_rcp85))
for I in CartesianIndices(sst_rcp85)
    sst_rcp85_anom[I] = sst_rcp85[I] - sst_rcp85_base[I[2], I[3]]
end

area = ncread("data/area.nc", "cell_area")
earth_surface_area = 5.10e14 # area of Earth

nyrs = 250
rad_response_area_mean = fill(NaN, nyrs)

for k in 1:nyrs
    rad_response_area_mean[k] = sum((rad_gf * vec(sst_rcp85_anom[k, :, :])) .*
        vec(area)) / earth_surface_area
end

fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, rad_response_area_mean)
ax.xlabel = L"Years $\,$"
ax.ylabel = L"Radiative response (W/m^2) $\,$"
save("radresponse-timeseries.png", fig)


copts = (colormap = cgrad([:cornflowerblue, :white, :red]), colorrange = (-10, 10))
rad_response_obs = zeros(nx, ny)

fig = Figure(size=(800, 450))
ga = GeoAxis(
    fig[1, 1];                  # any cell of the figure's layout
    dest = "+proj=natearth",    # the CRS in which you want to plot
    xticklabelspace = 50.0,
    xticklabelsvisible = false,
    yticklabelsvisible = false,
)
xlims!(ga, (-180, 180))
ylims!(ga, (-90, 90))

lons = range(-180, stop=180, length=nx)
lats = range(-90, stop=90, length=ny)

perturb_coord = (20, -60)
perturb_width = 80
perturb_height = 5
perturb_lon = perturb_coord[1] .+ (-perturb_width, perturb_width)
perturb_lat = perturb_coord[2] .+ (-perturb_height, perturb_height)

in_perturb_lon(lon) = perturb_lon[1] < lon < perturb_lon[2]
in_perturb_lat(lat) = perturb_lat[1] < lat < perturb_lat[2]
in_perturb_region(lon, lat) = in_perturb_lon(lon) & in_perturb_lat(lat)
perturb_mask = fill(false, nx, ny)
for i in 1:nx, j in 1:ny
    perturb_mask[i, j] = in_perturb_region(lons[i], lats[j])
end

Colorbar(fig[1, 2], height = Relative(0.5),
    label = L"TOA radiative anomaly $\mathrm{W \, m^{-2}}$"; copts...)

record(fig, plotsdir("radresponse.mp4"), 1:12:nyrs, framerate=12) do i
    rad_response_obs .= reshape(rad_gf * vec(sst_rcp85_anom[i, :, :]), nx, ny)
    # rad_response_obs[continent_mask'] .= NaN
    surface!(ga, lons, lats, reverse(rad_response_obs, dims=2),
        nan_color = :black, shading = NoShading; copts...)
    contour!(ga, lons, lats, reverse(continent_mask', dims=2), levels=[0.5],
        linewidth=2, color=:gray)
end


sst_perturb = zeros(ny, nx) # sst_rcp85_anom[50, :, :] # zeros(nx, ny)
sst_perturb[perturb_mask'] .-= 2

fig_perturb = Figure()
ga_perturb = GeoAxis(
    fig_perturb[1, 1];                  # any cell of the figure's layout
    dest = "+proj=natearth",            # the CRS in which you want to plot
    xticklabelsvisible = false,
    yticklabelsvisible = false,
)

copts_small = (colormap = cgrad([:cornflowerblue, :white, :red]), colorrange = (-3, 3))
rad_response_obs .= reshape(rad_gf * vec(sst_perturb), nx, ny)
sf = surface!(ga_perturb, lons, lats, reverse(rad_response_obs, dims=2), shading = NoShading; copts_small...)
contour!(ga_perturb, lons, lats, perturb_mask, levels=[0.5], linewidth=2, color=:black)
contour!(ga_perturb, lons, lats, reverse(continent_mask', dims=2), levels=[0.5], linewidth=2, color=:gray)
Colorbar(fig_perturb[1, 2], sf, height = Relative(0.5),
    label = L"TOA Radiative anomaly ($\mathrm{W \, m^{-2}}$)";)
save(plotsdir("sst-perturb.png"), fig_perturb)