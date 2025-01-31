using Documenter

push!(LOAD_PATH,"../src/")
using HPSAnalysis



@info "Making documentation..."


makedocs(
    sitename = "HPSAnalysis",
    authors = "Yannick Witzky",
    format = Documenter.HTML(;edit_link=:commit,
    prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
    modules = [HPSAnalysis],
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    pages = [
        "Home" => "index.md",
        "Setup" => "setup.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    devbranch="LF_Docs",
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    push_preview=true,
)=#
