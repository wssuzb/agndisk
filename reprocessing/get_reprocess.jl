using DelimitedFiles
using ProgressMeter
using Plots
using DataFrames
# using AccretionDisk
using ProgressMeter
using AstroLib
using FITSIO

include("disk.jl")

band = ["uw2", "um2", "uw1", "uuu", "ubb", "uvv"]

wave = [2000 ,3000, 4000, 5000, 6000.0, 7000, 8000]
t = disk.trans_curve("./transfile", band, wave; plot_trans_curve=true)

w_start = 100.0
w_end = 8000.0
wf = disk.wavelength_filter(w_start, w_end, t.trans, band)


par = disk.disk_parameter()

disk_par_update = disk.disk_systerm(par)

lc = FITS("./drw_UVW2.fits")[2]
mjd = read(lc, "mjd")
flux = read(lc, "flux")
err = read(lc, "err")
thindisk = disk.thin_disk(band, mjd, flux, err, wf, disk_par_update)
timerange = thindisk.idx

items = [[idx, item] for (idx, item) in zip(1:length(timerange),timerange)]
result = [Array{Any}(undef, 6) for i in 1:length(timerange)]

sixflux = Dict()
ill_sed = []



for x in items
    res = disk.reprocessing(thindisk, x, disk_par_update)
    idx = round(Int64, res.idx)
    result[idx] = res.ary
end

df = DataFrame()
df.idx = timerange
df[!, :time] .= 0.0
df[!, :uw2] .= 0.0
df[!, :um2] .= 0.0
df[!, :uw1] .= 0.0
df[!, :uuu] .= 0.0
df[!, :ubb] .= 0.0
df[!, :uvv] .= 0.0

for (i, r) in zip(1:length(result), eachrow(df))
    r.time  = result[i][1]
    r.uw2 = result[i][2]
    r.um2 = result[i][3]
    r.uw1 = result[i][4]
    r.uuu = result[i][5]
    r.ubb = result[i][6]
    r.uvv = result[i][7]

end

data = Dict(
    "time"=>df[!, :time],#.-2456600
    "uw2"=>df[!, :uw2], "um2"=>df[!, :um2], "uw1"=>df[!, :uw1],
    "uuu"=>df[!, :uuu], "ubb"=>df[!, :ubb], "uvv"=>df[!, :uvv]
            )

plot(data["time"], data["uw2"])
plot!(data["time"], data["um2"])
plot!(data["time"], data["uw1"])
plot!(data["time"], data["uuu"])
plot!(data["time"], data["ubb"])
plot!(data["time"], data["uvv"])

fid2 = FITS("./ngc5548.fits", "w")

myheader = FITSHeader(["Units"], ["erg cm-2 s-1 A-1"], ["FLUX units"])
sedheader = FITSHeader(["Units"], ["Flux erg s-1 A-1 s r-1"], ["SED units"])
write(fid2, data; name="simulation lightcurves", header = myheader)
close(fid2)