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
