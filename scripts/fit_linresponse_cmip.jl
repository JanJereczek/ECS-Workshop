include("intro.jl")

# Load files
file4x = datadir("data_martin/All_CMIP6_TOA_Tas_piC21Run_a4x.nc")
tas = ncread(file4x, "tas")
ecs4x = ncread(file4x, "ecs4x_t")
fns = readdir(datadir("data_martin/tas_global_average"), join = true)
models = ncread(file4x, "models")

# Define params
τ = [1, 10, 100]
β = 1.01
np = length(τ)
ne = length(fns)
t_ramp = 1/12:1/12:150
nt_ramp = length(t_ramp)

# Init arrays
response = zeros(nt_ramp)
c_mat = zeros(ne, np)
response_mat = zeros(nt_ramp, np)

# Compute matrix of c-coefficients
for i in axes(c_mat, 1)

    for j in axes(response_mat, 2)
        @. response_mat[:, j] = (1 / (log(β) + 1/τ[j])) *
            (exp(t_ramp * log(β)) - exp(-t_ramp / τ[j])) -
            τ[j] * (1 - exp(-t_ramp / τ[j]))
    end

    # for j in axes(response_mat, 2)
    #     @. response_mat[:, j] = (τ[j] / (log(β) + 1)) *
    #         (exp(t_ramp * log(β) / τ[j]) - exp(-t_ramp / τ[j])) # -
    #         # τ[j] * (1 - exp(-t_ramp / τ[j]))
    # end

    response .= view(ncread(fns[i], "tas"), 1:nt_ramp)
    response .-= mean(view(response, 1:60))
    view(c_mat, i, :) .= response_mat \ response
end

# Plot response of the last model
fig = Figure(size=(700, 500), fontsize = 20)
ax = Axis(
    fig[1, 1],
    xlabel = latexify("Time (yr)"),
    ylabel = latexify("Temperature anomaly (K)"),
    xticks = latexticks(0:50:150),
    yticks = latexticks(-2.5:2.5:7.5),
)
check = response_mat * c_mat[end, :]
lines!(ax, t_ramp, response, label = latexify("data"))
lines!(ax, t_ramp, check, label = latexify("fit"), linewidth = 3)
axislegend(ax, position = :rb)
save(plotsdir("cmip/rampresponse.png"), fig)

# Box plot C-matrix entries
fig = Figure(size=(700, 500), fontsize = 20)
ax = Axis(
    fig[1, 1],
    ylabel = L"$c_i \, \mathrm{(K \, yr^{-1})} $",
    xticks = latexticks(1:3),
    yticks = latexticks(-2:0.4:2),
)
ylims!(ax, (-2, 2))
cats = vec([j for j in axes(c_mat, 2), i in axes(c_mat, 1)]')
boxplot!(ax, cats, vec(c_mat))
save(plotsdir("cmip/cmatrix-box.png"), fig)

# responses = zeros(nt_ramp, ne)
# for i in axes(c_mat, 1)
#     for j in eachindex(τ)
#         @. responses[:, i] += (1 - exp(-t_ramp / τ[j])) * c[j] * τ[j]
#     end
# end
# lines(t_ramp, responses[:, 2])


# Init arrays
nt = 150
t = collect(1:nt)
response = zeros(nt)
c_mat = zeros(ne, np)
response_mat = zeros(nt, np)

# Compute matrix of c-coefficients
for i in axes(c_mat, 1)

    for j in axes(response_mat, 2)
        @. response_mat[:, j] = (1 - exp(-t / τ[j])) * τ[j]
    end

    response .= view(tas, :, i)
    response .-= mean(view(response, 1:60))
    view(c_mat, i, :) .= response_mat \ response
end

fig = Figure(size=(700, 500), fontsize = 20)
ax = Axis(
    fig[1, 1],
    xlabel = latexify("Time (yr)"),
    ylabel = latexify("Temperature anomaly (K)"),
    xticks = latexticks(0:50:150),
    yticks = latexticks(-2.5:2.5:7.5),
)
check = response_mat * c_mat[end, :]
lines!(ax, t_ramp, response, label = latexify("data"))
lines!(ax, t_ramp, check, label = latexify("fit"), linewidth = 3)
axislegend(ax, position = :rb)
save(plotsdir("cmip/rampresponse.png"), fig)

# Plot ECS for all models
ecs_analytic = c_mat * τ
fig = Figure(size=(800, 500), fontsize = 20)
ax = Axis(
    fig[1, 1],
    ylabel = latexify("ECS"),
    xticks = (1:ne, latexify.(models)),
    yticks = latexticks(-1.0:0.5:6.0),
    xticklabelrotation = π/2,
    xticklabelsize = 16,
)
idx = 1:2:ne
scatter!(ax, idx, ecs_analytic[idx], markersize = 20, label = latexify("Estimated"))
scatter!(ax, idx, ecs4x[idx] ./ 2, markersize = 20, label = latexify("Gregory"))
axislegend(ax, location = :rb)
save(plotsdir("cmip/ecs-gregory.png"), fig)
