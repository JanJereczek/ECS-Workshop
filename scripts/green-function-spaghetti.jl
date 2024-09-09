using NPZ
using NetCDF
using CairoMakie
using GeoMakie
using StatsBase

hasnan(x) = any(isnan, x)

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

sst_rcp85_base = reshape(mean(sst_rcp85[1:21, :, :], dims=1), ny, nx)
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

kobs = Observable(1)
fig = Figure()
ax = Axis(fig[1, 1])
hidedecorations!(ax)
copts = (colormap = cgrad([:cornflowerblue, :white, :red]), colorrange = (-5, 5)) #, nan_color = :black)
rad_response_obs = zeros(nx, ny)

fig = Figure()
ga = GeoAxis(
    fig[1, 1]; # any cell of the figure's layout
    dest = "+proj=wintri", # the CRS in which you want to plot
)
xlims!(ga, (-180, 180))
ylims!(ga, (-90, 90))

lons = range(-180, stop=180, length=nx)
lats = range(-90, stop=90, length=ny)

record(fig, "radresponse.mp4", 1:10) do i
    kobs[] = i
    rad_response_obs .= reshape(rad_gf * vec(sst_rcp85_anom[i, :, :]), nx, ny)
    rad_response_obs[continent_mask'] .= NaN
    surface!(ga, lons, lats, reverse(rad_response_obs, dims=2); copts...)
end