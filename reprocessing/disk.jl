module disk

using Unitful
using UnitfulAstro
import FITSIO
using ConfParser
using Plots
using DelimitedFiles
using Interpolations
using Statistics
using Pkg
using ProgressMeter
# include("struct/disk_struct.jl")
# include("struct/func_parameter.jl")
include("diskconfig.jl")
include("transconfig.jl")
include("reprocessing.jl")

export trans_curve, wavelength_filter, disk_parameter, disk_systerm, thin_disk, reprocessing
end # module
