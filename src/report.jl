const DATE_FORMAT = dateformat"yyyy-mm-dd H\hM\mS\s ZZZZ"
const DATE_FORMAT2 = dateformat"yyyy-mm-dd HH:MM:SS ZZZZ"

DIR = @__DIR__

PREAMBLE = """<head><meta charset="utf-8"></head>
<style>
$(read(joinpath(DIR,"style.css"), String))
</style>
"""

PART1 = mt"""# Analysis Results for Sample *{{{:samplename}}}*

Sample *{{:samplename}}* was found to contain the following.

| **Compound** | **Amount (w%)** |
|:------------ | -----------:|
{{#:matches}}| {{name}}          | {{percent}}      |
{{/:matches}}

## Sample Information
{{#:info}}
- **{{{key}}}:** {{{value}}}
{{/:info}}
"""

PLOT_HEADER = mt"""### {{{:component}}} component
- Library source: {{{:source}}}
- Fit Score: {{{:overall_fit}}}
"""

PART2 = mt"""## Signal Reconstruction and Residue

Portion of signal recovered by reconstruction: {{:recovery}}.

Solvent, internal standard, and chemical shift reference omitted for clarity.
"""

EMPTY = mt"""# Analysis Results for Sample *{{{:samplename}}}*

Sample *{{:samplename}}* did not match any of the spectra in the reference library.

## Sample Information
{{#:info}}
- **{{{key}}}:** {{{value}}}
{{/:info}}
"""

"""
    report(s::NMR.Spectrum, lib::Array{NMR.Spectrum,1}, d::NMR.DecompositionResult, quants, darks, is, path, short = true)

Generate an HTML report from decomposition results `d` and quantitation results `quants`, obtained from the analysis of NMR spectrum `s` analyzed against library group `lib`.

The report will be written to the directory `path`. `darks`, an array of NMR spectra containing 'dark' resonances, and `is`, the internal standard spectrum, are passed for reference.
If `short` is **true** a short report is generated.
"""
function report(s::NMR.Spectrum, lib::Array{NMR.Spectrum,1}, d::NMR.DecompositionResult, quants, darks, is,
                path, short = true)
    sampinf = sampleinfo(s)
    if haskey(sampinf, "Long report")
        short = !(lowercase(sampinf["Long report"]) ∈ ["true", "yes"])
    end
    mp = match_plots(s, lib, d, quants, short)
    plots = Dict( refno => IOBuffer() for refno=keys(mp) ) # Buffers to hold plot output
    for (k,v) in mp
        show(plots[k], MIME("text/html"), v)
    end
    matches = [ Dict( "name" => name(lib[k]),"percent" => statistics(v))
                for (k,v) in quants ]
        
    #augment sample info
    if haskey(s.acqu, "AUTOPOS")
        sampinf["Autosampler position"] = string(s["AUTOPOS"])
    end
    sampinf["NMRquant algorithm version"] = findversion("NMRquant")
    sampinf["Library version"] = findlibversion(lib)
    sampinf["Experiment number"] = string(s.expno)
    sampinf["Date of data collection"] = string(Dates.format(astimezone(ZonedDateTime(unix2datetime(s.acqu["DATE"]), tz"UTC"), localzone()), DATE_FORMAT2))
    sampinf["Analysis method"] = s.acqu["EXP"]
    
    info = [ Dict( "key" => k, "value" => v) for (k,v) in sorted_kv(sampinf) ]
    overall_fit = Dict{Int64, Float64}()
    for refno in keys(mp)
        refindex = findfirst(x->x==refno,d.refnums)
        overall_fit[refno] = *((fit_score for fit_score in d.fit_scores[refindex])...)^(1/length(d.fit_scores[refindex]))
    end
    _, recon, residue = NMR.decompose(d)
    recon_plt = reconplot(s, recon, residue, darks)
    recovery = @sprintf "%.1f" 100.0 * norm(recon,1) / norm(recon .+ residue,1)
    md1 = render(PART1, samplename=s.name, matches=matches, info=info)
    md2 = render(PART2, recovery=recovery)
    #testing
    #v3.3.2:
    #outfile = joinpath(path, "report " * zonedtime(now(), DATE_FORMAT) * ".html")
    #v3.4.0: Sample name in report file name
    outfile = joinpath(path, "$(get(sampinf, "Sample", "report")) " * zonedtime(now(), DATE_FORMAT) * ".html")
    open(outfile, "w") do f
        write(f, PREAMBLE)
        if !isempty(quants)
            write(f, Markdown.html(Markdown.parse(md1)))
        else
            write(f, Markdown.html(Markdown.parse(render(EMPTY, samplename=s.name, info=info))))
        end
        for (refno,p) in plots
            write(f, Markdown.html(Markdown.parse(render(PLOT_HEADER, component=name(lib[refno]), source=get(sampleinfo(lib[refno]),"Source",""), overall_fit=score_format(overall_fit[refno])))))
            write(f, p.data)
        end
        if !short
            write(f, Markdown.html(Markdown.parse(md2)))
            show(f, MIME("text/html"), recon_plt)
        end
    end
    outfile
end

LIBRARY_PREAMBLE = """<head><meta charset="utf-8"></head>
<style>
$(read(joinpath(DIR,"librarystyle.css"), String))
</style>
"""

LIBRARY_TITLE = mt"""# DAS NMRquant {{:name}} Library v{{:version}} Integrals
"""

COLUMN_WIDTHS = """
    <col width='16%'>
    <col width='16%'>
    <col width='12%'>
    <col width='56%'>
"""

TABLE_TITLE = """<tr>
        <th>Compound</th>
        <th><font size='-1'>Number of Integral(s)</font></th>
        <th><font size='-1'>Number <sup>1</sup>H(s) represented in each Integral Window</font></th>
        <th>Integral Windows</th>
    </tr>
</table>
"""

LIBRARY_TABLE="""
    <tr>
        <td>{{:name}}</td>
        <td>{{:numint}}</td>
        <td><font size='-1'>{{:intprotons}}</font></td>
"""

function library_report(lib::Array{NMR.Spectrum,1}, path)
    lp = library_plots(lib)
    plots = Dict( refno => IOBuffer() for refno=keys(lp) ) # Buffers to hold plot output
    for (k,v) in lp
        show(plots[k], MIME("text/html"), v)
    end
    
    outfile = joinpath(path, "DAS NMRquant $(lib[1].name) Library v$(findlibversion(lib)) " * zonedtime(now(), DATE_FORMAT) * ".html")
    title = render(LIBRARY_TITLE, name=lib[1].name, version=findlibversion(lib))
    open(outfile, "w") do f
        write(f, LIBRARY_PREAMBLE)
        write(f, Markdown.html(Markdown.parse(title)))
        write(f, "<table>", COLUMN_WIDTHS, TABLE_TITLE)
        write(f, "<table>", COLUMN_WIDTHS)
        
        libnum = Dict(name(lib[x]) => x for x = keys(plots))
        for n in sort(collect(keys(libnum)))
            ref = lib[libnum[n]]
            
            intprotons = ""
            intprotons = "[" * string(round(Int, protons(ref)[1]))
            if length(protons(ref)) > 1
                for i = 2:length(protons(ref))
                    intprotons *= ","*string(round(Int, protons(ref)[i]))
                end
            end
            intprotons *= "]"
            
            write(f, render(LIBRARY_TABLE, name=n, numint=length(protons(ref)), intprotons=intprotons), "<td>", plots[libnum[n]].data, "</td></tr>")
        end
        write(f, "</table>")
    end
    outfile
end


VALIDATION_PREAMBLE = """<head><meta charset="utf-8"></head>
<style>
$(read(joinpath(DIR,"validationstyle.css"), String))
</style>
"""

VRPART1 = """# Validation Report

- Date: {{{:reportdate}}}
- Program version: {{{:version}}}
- {{{:libname}}} Library version: {{{:libversion}}}

## Summary:

| **Test** | **Result** |
|:------------ | :-----------:|
{{#:r}}| {{test}}          | {{{result}}}      |
{{/:r}}
"""

VRPART2 = """### Quantitation Results:

| **Sample** | **Expected** | **Result** |
| :--------- | :---------- | :-------- |
{{#:r}}|  {{{sample}}} | {{#expected}} {{.}} <br> {{/expected}} | <font color=\"{{#passed}}green{{/passed}}{{^passed}}red{{/passed}}\"> {{#result}} {{.}} <br> {{/result}} </font> |
{{/:r}}

## Test log:
"""

function validation_report(lib::Array{NMR.Spectrum,1}, results::Dict{String, Bool}, quantlist::Dict, quantresults, testlogs::Dict{String,IOBuffer}, outpath)
    version = findversion("NMRquant")
    libname = lib[1].name
    libversion = findlibversion(lib)
    reportdate = now()
    outfile = joinpath(outpath, "Validation Report NMRquant v$version $libname v$libversion $(zonedtime(reportdate, DATE_FORMAT)).html")
    open(outfile, "w") do f
        write(f, VALIDATION_PREAMBLE)
        
        r = [ Dict( "test" => k,"result" => v ? "<font color=\"green\">Passed</font>" : "<font color=\"red\">Failed</font>")
            for (k,v) in sorted_kv(results) ]
        vrp1 = Markdown.html(Markdown.parse(render(VRPART1, reportdate=zonedtime(reportdate, DATE_FORMAT2),version=version, libname=libname, libversion=libversion,r=r)))
        write(f, restorehtml(vrp1))
        
        r = [ Dict( "sample" => k,
                    "expected" => haskey(quantlist[k],0) ? ["None"] : ["$(ev)% $(sampleinfo(lib[en])["Name"])" for (en,ev) in quantlist[k]],
                    "result" => haskey(quantresults[k][1],0) ? ["None"] : ["$(statistics(qv))% $(sampleinfo(lib[qn])["Name"])" for (qn, qv) in quantresults[k][1]],
                    "passed" => quantresults[k][2] ? "true" : "")
            for (k,_) in sorted_kv(quantlist) ]
        temp = Markdown.html(Markdown.parse(render(VRPART2, r=r)))
        write(f, restorehtml(temp))
        
        write(f, "<pre>")
        for k in sort(collect(keys(testlogs)))
            write(f, Markdown.html(Markdown.parse("### $(k)\n")))
            
            #Writes with ANSI codes converted to HTML
            write(f, ansitohtml(String(take!(copy(testlogs[k])))))
        end
        write(f, "</pre>")
        close(f)
    end
    outfile
end
