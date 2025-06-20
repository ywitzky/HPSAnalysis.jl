using Documenter

push!(LOAD_PATH,"../src/")
using HPSAnalysis

@info "Making documentation..."


makedocs(
    sitename = "HPSAnalysis",
    authors = "Yannick Witzky",
    repo=Remotes.GitHub("ywitzky", "HPSAnalysis.jl"),
    format = Documenter.HTML(;edit_link=:commit, prettyurls = true), #get(ENV, "CI", nothing) == "true"), # easier local build
    #size_threshold_ignore = ["listfunctions.md"], # easier local build
    modules = [HPSAnalysis],
    pages = [
        "Home" => "index.md",
        "Setup" => "setup.md",
        "Slab Analysis" => "SlabAnalysis.md",
        "List of functions" => "listfunctions.md",
        "Installation" => "Requirements.md",
        "Deprecated" => "Deprecated.md"
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    push_preview=true,
)
