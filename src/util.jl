name(s::NMR.Spectrum) = sampleinfo(s)["Name"]
MW(s::NMR.Spectrum) = parse(Float64, sampleinfo(s)["MW"])
blank_conc(s::Spectrum) = parse(Float64, sampleinfo(s)["Blank conc. (mg/mL)"])

function protons(s::NMR.Spectrum)
    pr = get(sampleinfo(s), "Protons", "")
    parsed = eval(Meta.parse(pr))
    parsed !== nothing ? float.(parsed) : nothing
end

function quant_threshold(s::NMR.Spectrum)
    th = get(sampleinfo(s), "Quantitation threshold (%)", "0.0")
    parse(Float64, th)
end

# (internal standard) in the Name field signifies the internal standard.
isstandard(s::NMR.Spectrum) = occursin("(internal standard)", name(s))

# A spectrum is dark if its title doesn't contain the Proton field
# or if it is the internal standard.
isdark(s::NMR.Spectrum) = !haskey(sampleinfo(s), "Protons") || isstandard(s)

"""
    library(path)
    
Returns a tuple from a valid directory `path`: `lib`, an array of library spectra; `darks`, an array of 'dark' spectra; and `is`, the internal standard spectrum.
"""
function library(path)
    try
        expnos = [expno for expno in readdir(path) if tryparse(Int, expno) !== nothing]
        sort!(expnos, by=expno->parse(Int,expno))
        lib = [NMR.Spectrum(joinpath(path, expno),1) for expno in expnos]
        for (i,l) in enumerate(lib)
            pr = protons(l)
            if pr === nothing
                continue
            end
            l1 = length(NMR.intrng_shifts(l))
            l2 = length(pr)
            if l1 != l2
                error("Reference $(name(l)) ($i) contains $l1 integration regions != $l2 proton groups.")
            end
        end
        is = filter(isstandard, lib)
        if length(is) == 0
            error("$path: No internal standard found.")
        end
        if length(is) != 1
            error("One internal standard expected; found $(length(is)): $(name.(is))")
        end
        darks = filter(isdark, lib)
        lib, darks, is[1]
    catch e
        println(sprint(showerror, e)*"\nIgnoring directory $path.")
        return nothing
    end
end

function libraries(path)
    result = Dict( p => library(joinpath(path, p))
        for p in readdir(path) if isdir(joinpath(path, p)) )
    return filter(x -> x.second != nothing, result)
end

function getlib(s::Spectrum, libs)
    libname = get(sampleinfo(s), "Library", "D2O")
    libs[libname]
end

"""
    sampleinfo(s::NMR.Spectrum)
    
Returns a Dictionary of sample information obtained from the title of NMR spectrum `s`.
"""
function sampleinfo(s::NMR.Spectrum)
    result = Dict{String,String}()
    for l in split(NMR.title(s), '\n')
        if findfirst("=", l) !== nothing
            k,v = split(l, '='; limit = 2)
            result[strip(k)] = strip(v)
        end
    end
    result
end

function sampleinfovalid(s::NMR.Spectrum, libs)
    sinfo = sampleinfo(s)
    
    for inputfield in ["Sample weight (mg)", "Volume (mL)", "Blank conc. (mg/mL)", "Library"]
        if !haskey(sinfo, inputfield)
            error("$inputfield missing from Spectrum Title")
        end
    end
    
    for inputfield in ["Sample weight (mg)", "Volume (mL)", "Blank conc. (mg/mL)"]
        if tryparse(Float64, sinfo[inputfield]) == nothing
            error("Spectrum Title: Invalid $(inputfield)!")
        end
    end
    
    if !(sinfo["Library"] ∈ collect(keys(libs)))
        error("Spectrum Title: Invalid 'Library'")
    end
    
    return true
end

function statistics(percentages)
    m = mean(percentages)
    pm = (maximum(percentages)-minimum(percentages))/2
    if m > 10.0
        @sprintf "%.1f" m
    else
        @sprintf "%.2f" m
    end
end

ppm_format(δ) = @sprintf "%.2f" δ

score_format(score) = @sprintf "%.3f" score

findlibversion(lib) = get(sampleinfo(lib[1]), "Library version", "Missing")

findversion(name::String) = string(Pkg.installed()[name])

zonedtime(time::DateTime=now(), f=dateformat"yyyy-mm-dd H\hM\mS\s ZZZZ") = Dates.format(ZonedDateTime(time,localzone()), f)

function ansitohtml(text::String)
    rd = Dict(r"\e\[0m" => "",
    r"\e\[1m" => "<b>",
    r"\e\[22m" => "</b>",
    r"\e\[37m" => "<font color=\"silver\">",
    r"\e\[39m" => "</font>",
    r"\e\[91m" => "<font color=\"red\">",
    r"\e\[32m" => "<font color=\"green\">",
    r"\e\[36m" => "<font color=\"teal\">")
    
    for i in rd
        text = replace(text, i)
    end
    text
end

function restorehtml(text::String)
    dict = Dict(r"&lt;" => "<",
    r"&#61;" => "=",
    r"&quot;" => "\"",
    r"&gt;" => ">",
    r"&#x2F;" => "/")
    
    for i in dict
        text = replace(text, i)
    end
    text
end
