using Documenter

push!(LOAD_PATH,"../src/")
using HPSAnalysis

@info "Making documentation..."


makedocs(
    sitename = "HPSAnalysis",
    authors = "Yannick Witzky",
    format = Documenter.HTML(;
    size_threshold_ignore = ["listfunctions.md"]
    ), # easier local build
    modules = [HPSAnalysis],
	remotes=nothing,
    pages = [
        "Home" => "index.md",
        "Slab Analysis" => "SlabAnalysis.md",
        "List of functions" => "listfunctions.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
