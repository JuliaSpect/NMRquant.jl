import Markdown

using Colors: @colorant_str
using Measures: mm
using NMR: intrng_shifts, intrng_indices
using Plots: @layout, grid

@userplot MatchPlot
@userplot RangePlots
@userplot ReconPlot
@userplot LibraryPlot

@recipe function f(m::MatchPlot; short = true, dark = false)
    s, ref, re_ft, fit_scores, weights = m.args
    labels = ["Sample" "Reference"]
    Δ = chemical_shifts(s)

    framestyle --> :box
    grid --> nothing
    legend --> :inside
    linewidth --> 1.0
    xflip := true
    yticks --> []
    label --> ""

    if dark
        # plot overall match
        xticks --> plot_xticks(Δ, 1.25)
        Δ, [ref[:], re_ft]
    else    
        RangePlots([s, ref, re_ft, fit_scores, weights, short])
    end
end

@recipe function f(r::RangePlots)
    s, ref, re_ft, fit_scores, weights, short = r.args

    #v3.4.1:
    #=rngs = [ (δ, [s[j], re_ft[j]]) for (δ, j) in
             zip(intrng_shifts(ref, NMR.TOL), intrng_indices(ref, NMR.TOL)) ] =#
			 
	#v.3.4.2: Adjusted integral indices to be relative to s rather than ref, to account for slight differences in SW in s
    rngs = [ (δ, [s[j], re_ft[j]]) for (δ, j) in
             zip(intrng_shifts(ref, NMR.TOL), [ ppmtoindex(s,x) for x in intrng_shifts(ref, NMR.TOL) ]) ]
    n = length(rngs)
    
    if !short
        layout := @layout [ a
                        grid(1,n){0.2h}]
        size --> (750, 750)
        sps = 1
    else
        layout := (1, n)
        size --> (750, 150)
        sps = 0
    end
    
    ovmax, ovmin, maxspan = 0, 0, 0
    totweight = sum(weights)
    
    for (i, rng) in enumerate(rngs)
        δ, data = rng
        ymin = minimum(min.(data...))
        ymax = maximum(max.(data...))
        span = ymax - ymin
        
        @series begin
            left_margin := -3.0mm
            right_margin := -3.0mm
            top_margin := 2mm
            subplot := sps + i
            xticks := subplot_xticks(δ)
            ylims := (ymin - 0.1 * span, ymax + 0.1 * span)
            titlefontsize := 8
            
            if i == 1
                title := "Fit: $(score_format(fit_scores[i]))\nWeight: $(score_format(weights[i]/totweight))"
            else
                title := "$(score_format(fit_scores[i]))\n$(score_format(weights[i]/totweight))"
            end 
            
            δ, [data[1], data[2]]
        end
        
        ovmax = max(ymax, ovmax)
        ovmin = min(ymin, ovmin)
        maxspan = max(span, maxspan)
    end
    
    if !short
        @series begin
            Δ = NMR.chemical_shifts(s)
            xl = plot_xlimits(ref)
            
            subplot := 1
            xlims := xl
            xticks := plot_xticks(xl)
            ylims := ovmin - 0.1 * maxspan, ovmax + 0.5 * maxspan
            Δ, [s[:], re_ft]
        end
    end
end

@recipe function f(r::ReconPlot)
    s, recon, residue, dark_refs = r.args
    Δ = NMR.chemical_shifts(s)
    recon = copy(recon)
    orig = copy(s[:])
    for dr in dark_refs
        for rng in ppmtoindex.(Ref(s), intrng(dr))
            recon[rng] .= 0.0
            residue[rng] .= 0.0
            orig[rng] .= 0.0
        end
    end

    framestyle --> :box
    size --> (750, 600)
    xflip --> true
    yformatter --> :scientific
    yticks := []
    xl = plot_xlimits(s, 0.1)
    xlims := xl
    xticks := plot_xticks(xl)

    @series begin
        label --> ["Sample" "Reconstruction"]
        Δ, [orig, recon]
    end
    @series begin
        seriescolor --> colorant"forestgreen"
        label --> "Residue"
        Δ, residue
    end
end

@recipe function f(l::LibraryPlot)
    ref, = l.args
    rngs = [ (δ, ref[j]) for (δ, j) in
             zip(intrng_shifts(ref), intrng_indices(ref)) ]
    n = length(rngs)
    
    layout := (1,n)
    size := (600, 130)
    xflip := true
    grid --> nothing
    yticks --> []
    label --> ""
    left_margin := -2.5mm
    right_margin := -2.5mm
    top_margin := -3.0mm
    seriescolor := 2
    
    for (i, rng) in enumerate(rngs)
        δ, data = rng
        @series begin
            subplot := i
            xticks := subplot_xticks(δ)
            framestyle --> :box
            if n >= 10 xtickfontsize := 6 end
                     
            δ, data
        end
    end
end

function subplot_xticks(Δ::AbstractRange)
    l = Δ.len
    i1, i2 = Δ[Int(ceil(0.2l))], Δ[Int(ceil(0.8l))]
    ([i1, i2], [ppm_format(i1), ppm_format(i2)])
end

function plot_xticks(Δ::AbstractRange, δ::Float64)
    m, M = extrema(Δ)
    ticks = (δ * ceil(m/δ)):δ:(δ * floor(M/δ))
    (ticks, map(string, ticks))
end

function plot_xticks(Δ::Tuple, δ::Float64 = round((Δ[2]-Δ[1])/10; digits=2))
    ticks = (δ * ceil(Δ[1]/δ)):δ:(δ * floor(Δ[2]/δ))
    (ticks, map(ppm_format, ticks))
end

function plot_xlimits(s::Spectrum, margin = 0.25)
    Δrng = intrng(s)
    if isempty(Δrng)
        O1P = s.acqu["O1P"]
        SW = s.acqu["SW"]
        return (O1P - SW/2, O1P + SW/2)
    else
        return (minimum(Δrng)[2] - margin, maximum(Δrng)[1] + margin)
    end
end

function match_plots(s::NMR.Spectrum, lib::Array{NMR.Spectrum,1}, d::NMR.DecompositionResult, quants, short = true)
    res = Dict()
    x = d.matrix .* d.coefficients'
    for (i,refno) in enumerate(d.refnums)
        l = lib[refno]
        if haskey(quants, refno)
            res[refno] = matchplot(s, l, (@view x[:,i]), d.fit_scores[i], d.weights[i], short = short)
        end
    end
    res
end

function library_plots(lib::Array{NMR.Spectrum,1})
    res = Dict()
    for (i, ref) in enumerate(lib)
        if !isdark(ref)
            res[i] = libraryplot(ref)
        end
    end
    res
end
