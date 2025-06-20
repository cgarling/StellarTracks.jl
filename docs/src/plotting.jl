import PyPlot as plt
plt.ioff()
ENV["MPLBACKEND"] = "agg"
import PyPlot: @L_str # For LatexStrings
plt.rc("text", usetex=true)
plt.rc("font", family="serif", serif=["Computer Modern"], size=16)
# This gets close but not quite
# plt.matplotlib.rcParams["axes.formatter.use_mathtext"] = true
# plt.rc("font", family="serif", serif=["cmr10"], size=14)
plt.rc("figure", figsize=(5,5))
plt.rc("patch", linewidth=1, edgecolor="k", force_edgecolor=true)
# https://matplotlib.org/stable/gallery/images_contours_and_fields/interpolation_methods.html
plt.rc("image", interpolation="none")

function plot_hr(iso)
    fig,ax1 = plt.subplots() # hide
    ax1.plot(iso.logTe, iso.Mbol) # hide
    ax1.set_xlim([3.85, 3.5]) # hide
    ax1.set_ylim(reverse(ax1.get_ylim())) # hide
    ax1.set_xlabel("logTe") # hide
    ax1.set_ylabel("Mbol") # hide
    return fig # hide
end

function plot_cmd(iso; xfilters=(:F090W, :F150W), yfilter=:F090W, xlim=[0.4, 1.62])
    fig,ax1 = plt.subplots() # hide
    ax1.plot(getproperty(iso, xfilters[1]) .- getproperty(iso, xfilters[2]),
             getproperty(iso, yfilter)) # hide
    ax1.set_ylim(reverse(ax1.get_ylim())) # hide
    ax1.set_xlim(xlim) # hide
    ax1.set_xlabel(string(xfilters[1]) * " - " * string(xfilters[2])) # hide
    ax1.set_ylabel(string(yfilter)) # hide
    return fig # hide
end

# Plot mutiple isochrones on same HR diagram
# p=BaSTIv1Library(0, true)
# mh = -2.5:0.1:0.0
# isos = [isochrone(p, 9, i) for i in mh]
# plot_hrs(isos, mh)
function plot_hrs(isos, mh)
    fig,ax1 = plt.subplots()
    cmap = plt.cm.gnuplot2
    norm = plt.matplotlib.colors.Normalize(vmin=minimum(mh), vmax=maximum(mh))
    for (i, iso) in enumerate(isos)
        ax1.plot(iso.logTe, iso.Mbol, c=cmap(norm(mh[i])))
    end
    # ax1.set_xlim([3.85, 3.5])
    ax1.set_xlim(reverse(ax1.get_xlim()))
    ax1.set_ylim(reverse(ax1.get_ylim()))
    ax1.set_xlabel("logTe")
    ax1.set_ylabel("Mbol")
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), ax=ax1)
    return fig
end
