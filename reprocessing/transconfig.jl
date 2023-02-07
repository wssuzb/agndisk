
mutable struct trans_curve_info
    band::Vector{String}
    weff::Vector{Float64}
end


function Base.show(trans_info::trans_curve_info)
    println("Band used in reprocessing model: ", trans_info.band)
    println("Band W_eff: ", trans_info.weff)
end


function trans_curve(filepath::String, band::Vector{String}, weff::Vector{Float64}; plot_trans_curve::Bool=true)

    println("The transmission file located in: ", filepath)
    println("Transmission curve file name: ", readdir(filepath))
    
    trans_info = trans_curve_info(
        band, weff
    )

    show(trans_info)

    if plot_trans_curve
        p = plot(
            legend = true, 
            framestyle = :box,
            # aspect_ratio = 1, 
            size = (400, 200),
            format = :svg
        )        
        for i in 1:length(trans_info.band)
            t = readdlm(filepath * "/total_" * band[i] * ".txt")
            plot!(p, t[:, 1], t[:, 2], label=trans_info.band[i])
        end
        display(p)
        trans = Dict()
        for i in 1:length(trans_info.band)
            t = readdlm(filepath * "/total_" * band[i] * ".txt")
            trans[trans_info.band[i]] = (wavelength=t[:, 1], strength=t[:, 2])
        end
    else
        trans = Dict()
        for i in 1:length(trans_info.band)
            t = readdlm(filepath * "/total_" * band[i] * ".txt")
            trans[trans_info.band[i]] = (wavelength=t[:, 1], strength=t[:, 2])
        end
    end
    return (trans=trans, band=band, weff=weff)
end


function wavelength_filter(w_start::Float64, w_end::Float64,transmission, band::Vector{String})
    wavelength = range(w_start, stop=w_end,length=Int(w_end - w_start) + 1)
    six_filter = zeros((Int(w_end - w_start) + 1, length(transmission)+1))
    six_filter[:, 1] = wavelength

    for i=1:length(transmission)
        tran_w = transmission[band[i]].wavelength # wavelength
        resp = transmission[band[i]].strength # specresp
        itp = LinearInterpolation(tran_w, resp, extrapolation_bc=Flat())
        six_filter[:, i+1] = itp(wavelength)
    end
    return six_filter
end