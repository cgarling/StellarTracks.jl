module MIST

# imports from parent module
using ..StellarTracks: AbstractTrack, AbstractTrackSet, AbstractTrackLibrary, uniqueidx, Mbol
import ..StellarTracks: mass, post_rgb, isochrone
import ..StellarTracks: X, X_phot, Y, Y_phot, Z, Z_phot, MH, chemistry

# Imports for data reading / processing
using DataDeps: register, DataDep, @datadep_str, unpack
import CSV
# using DelimitedFiles: readdlm
using TypedTables: Table

# We will simply reuse the MISTChemistry defined in BolometricCorrections
using BolometricCorrections.MIST: MISTChemistry, unpack_txz

# Imports for core module code
using ArgCheck: @argcheck
using StaticArrays: SVector

# const eep_idxs = (PMS_BEG = 1,      # beginning of PMS; log(T_c) = 5
#                   MS_BEG = 202,     # beginning of MS; H-burning L > 99.9% of total L
#                   IAMS = 353,       # intermediate age MS; X_c = 0.3
#                   MS_TO = 454,      # terminal age MS; X_c = 1e-12
#                   RG_TIP = 605,     # max(L) or min(Teff) after core H-burning, before core-He burning
#                   HE_BEG = 631,     # beginning core He fusion to end core He fusion
#                   END_CHEB = 707,   # end of core He burning; Y_c = 1e-4
#                   TPAGB_BEG = 808,  # For low mass stars, this is the beginning of TP-AGB;
#                                     # Y_c < 1e-6, H and He shells comparable mass;
#                                     # for high mass stars, this is the end of core carbon burning
#                   END_TPAGB = 1409, # For stars that will become white dwarfs, end of TP-AGB phase
#                   WD = 1710)        # White dwarfs
# # Indices where the different phases begin

# const eep_lengths = NamedTuple{keys(eep_idxs)[begin:end-1]}(eep_idxs[i+1] - eep_idxs[i] for i in eachindex(values(eep_idxs))[begin:end-1])

const eep_lengths = (PMS_BEG = 201,   # beginning of PMS; log(T_c) = 5
                     MS_BEG = 151,    # beginning of MS; H-burning L > 99.9% of total L
                     IAMS = 101,      # intermediate age MS; X_c = 0.3
                     MS_TO = 151,     # terminal age MS; X_c = 1e-12
                     RG_TIP = 26,     # max(L) or min(Teff) after core H-burning, before core-He burning
                     HE_BEG = 76,     # beginning core He fusion to end core He fusion
                     END_CHEB = 101,  # end of core He burning; Y_c = 1e-4
                     TPAGB_BEG = 601, # For low mass stars, this is the beginning of TP-AGB;
                                      # Y_c < 1e-6, H and He shells comparable mass;
                                      # for high mass stars, this is the end of core carbon burning
                     END_TPAGB = 301) # For stars that will become white dwarfs, end of TP-AGB phase
                                      # White dwarfs at 1710, final EEP point
# Indices where the different phases begin
const eep_idxs = NamedTuple{keys(eep_lengths)}((1, (cumsum(values(eep_lengths)[begin:end-1]) .+ 1)...))
# Which columns to actually keep after reading track file; used in Track and track_table
const select_columns = SVector(:star_age, :star_mass, :log_L, :log_Teff, :log_g, :log_surf_cell_z)
const track_type = Float64 # Float type to use to represent values


##########################################################################

# Data download, organization, etc.
include("init.jl")

##########################################################################


end # module
