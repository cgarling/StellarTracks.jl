using CairoMakie
CairoMakie.activate!(type="svg")

function plot_hr(iso)
    fig = Figure(size=(500, 500)) # hide
    ax = Axis(fig[1,1], xlabel="logTe", ylabel="Mbol", # hide
              xreversed=true, yreversed=true, # hide
              limits=(3.5, 3.85, nothing, nothing)) # hide
    lines!(ax, iso.logTe, iso.Mbol) # hide
    return fig # hide
end

function plot_cmd(iso; xfilters=(:F090W, :F150W), yfilter=:F090W, xlim=(0.4, 1.62))
    colors = getproperty(iso, xfilters[1]) .- getproperty(iso, xfilters[2]) # hide
    mags = getproperty(iso, yfilter) # hide
    fig = Figure(size=(500, 500)) # hide
    ax = Axis(fig[1,1], # hide
              xlabel=string(xfilters[1]) * " - " * string(xfilters[2]), # hide
              ylabel=string(yfilter), # hide
              yreversed=true, # hide
              limits=(xlim[1], xlim[2], nothing, nothing)) # hide
    lines!(ax, colors, mags) # hide
    return fig # hide
end

# Plot multiple isochrones on same HR diagram
# p=BaSTIv1Library(0, true)
# mh = -2.5:0.1:0.0
# isos = [isochrone(p, 9, i) for i in mh]
# plot_hrs(isos, mh)
function plot_hrs(isos, mh)
    fig = Figure(size=(600, 500))
    ax = Axis(fig[1,1], xlabel="logTe", ylabel="Mbol",
              xreversed=true, yreversed=true)
    mh_min, mh_max = minimum(mh), maximum(mh)
    for (i, iso) in enumerate(isos)
        t = (mh[i] - mh_min) / (mh_max - mh_min)
        lines!(ax, iso.logTe, iso.Mbol, color=t, colormap=:gnuplot2, colorrange=(0, 1))
    end
    Colorbar(fig[1,2], colormap=:gnuplot2, limits=(mh_min, mh_max))
    return fig
end
