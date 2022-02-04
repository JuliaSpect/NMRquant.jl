module NMRquant

using Dates
using LinearAlgebra: norm
using Mustache
using NMR
using RecipesBase
using Pkg
using Printf
using Sockets
using Statistics
using TimeZones

function write_config(path)
    open(joinpath(path, "config"), "w") do f
        write(f, match(r"inet ([\d\.]+)", read(`ip address show enp0s3`, String)).captures[1])
    end
end

function sorted_kv(dict,args...)
    k = sort(collect(keys(dict)), args...)
    ((x, dict[x]) for x in k)
end

function dark_signals(s::NMR.Spectrum, dark_areas)
    darks = NMR.extract.(s, [area[1] for area in dark_areas])
    for (d,name) in zip(darks, [area[2] for area in dark_areas])
        d[d.default_proc].title = "Name = $name"
    end
    darks
end

function remove_darks!(s::NMR.Spectrum, darks::Array{NMR.Spectrum,1})
    d = NMR.decompose(NMR.lsq_analyze(s, darks)) # find dark signals in s
    s[:] = d[3] # residue (non-dark part)
end

function remove_darks!(lib)
    darks = [l for l in lib if !haskey(sampleinfo(l), "Protons")]
    for l in lib
        if haskey(sampleinfo(l), "Protons")
            remove_darks!(l, darks)
        end
    end
end


"""
    analyze(s::NMR.Spectrum, lib::Array{NMR.Spectrum,1}, stan::Spectrum; out = stdout)
    
Perform analysis on the NMR spectrum `s` against the library group `lib`, using the internal standard spectrum `stan` for quantitation.

Returns the tuple (`decomposition results`, `quantitation results`) of type (DecompositionResult, Dict{Int, Array{Float64}}).
`Quantitation results` is a dictionary with library reference numbers as keys, and an array of % quantitations as values.
"""
function analyze(s::NMR.Spectrum, lib::Array{NMR.Spectrum,1}, stan::Spectrum; out = stdout, tol = NMR.TOL)
    callback(refs) = println(out, "Found $(join([name(lib[r]) for r in refs], ", ")).")
    d = NMR.lsq_analyze(s, lib; callback=callback, tol=tol)
    quants = quantitate(s, lib, d, stan)
    (d, quants)
end

"""
    quantitate(s::NMR.Spectrum, lib::Vector{NMR.Spectrum}, d::NMR.DecompositionResult, stan::Spectrum, minquant = quant_threshold(s))

Return a Dictionary of quantitation results based on DecompositionResult `d`, reported as percent by weight based internal standard spectrum `stan`,
library spectra `lib` and unknown spectrum `s`. An identification is only reported if the result is above the quantitation threshold `minquant`.
"""
function quantitate(s::NMR.Spectrum, lib::Vector{NMR.Spectrum}, d::NMR.DecompositionResult, stan::Spectrum, minquant = quant_threshold(s))
    results = Dict{Int,Array{Float64}}()
    sampinfo = sampleinfo(s)
    stanint = NMR.integrate(s, stan[stan.default_proc].intrng[1])
    stanconc = blank_conc(s) / MW(stan)
    stanprotons = protons(stan)[1]
    sampleweight = parse(Float64, sampinfo["Sample weight (mg)"])
    volume = parse(Float64, sampinfo["Volume (mL)"])
    for (refno, ref, coeff) in zip(d.refnums, lib[d.refnums], d.coefficients)
        if !isdark(ref)
            pr = protons(ref)
            mw = MW(ref)
            q = NMR.integrate(ref)./pr.*(coeff/stanint*stanprotons*mw*stanconc*volume/sampleweight*100.0)
            if mean(q) > max(minquant, quant_threshold(ref)) # only report more than meanquant %
                results[refno] = q
            end
        end
    end
    results
end


"""
    serve(port, home, libs)
    
Establish a listening port specified by `port` for receiving NMR spectra data located in a directory defined by `home`.

Valid NMR spectra will be analyzed against a library group `lib`, which must be a Dictionary with library names as keys,
and ([Library Spectrum], [Dark Spectrum], IS Spectrum) as values. A report will be generated at the location of the NMR spectra.
"""
function serve(port, home, libs)
    @async begin
        server = listen(ip"0.0.0.0", port)
        while true
            sock = accept(server)
            @info "Connection opened"
            l = strip(readline(sock))
            @info "Received $l"
            if strip(l) == "bye"
                @info "Shutting down server"
                close(sock)
                close(server)
                return 0
            end
            try
                s = NMR.Spectrum(joinpath(home,l),1)
                sampleinfovalid(s, libs)
                lib, darks, is = getlib(s, libs)
                
                #3.4
                if lowercase(get(sampleinfo(lib[1]), "Reference Deconvolution", "false")) ∈ ["true", "yes"]
                    tm1 = get(sampleinfo(lib[1]), "tm1", NMR.TM1)
                    tm2 = get(sampleinfo(lib[1]), "tm2", NMR.TM2)
                    s = refdecon(s, lib[3];tm1=tm1, tm2=tm2)
                    s[1].procno = 99
                    s[1].title *= "\r\ntm1 = $tm1\r\ntm2 = $tm2"
                    dump(joinpath(home,l,"pdata","99"),s[1])
                end
                
                d,q = analyze(s, lib, is; out = sock)
                outfile = report(s, lib, d, q, darks, is, joinpath(home,l))
                println(sock, "Done")
            catch err
                @error err
                @error stacktrace(catch_backtrace())
                println(sock, err)
            finally
                close(sock)
                @info "Connection closed, listening"
            end
        end
    end
end

include("plots.jl")
include("report.jl")
include("util.jl")

export serve, analyze, report, library_report, quantitate, sampleinfo, library

end
