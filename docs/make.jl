using Documenter

push!(LOAD_PATH,"../src/")
using HPSAnalysis



@info "Making documentation..."


makedocs(
    sitename = "HPSAnalysis",
    authors = "Yannick Witzky",
    format = Documenter.HTML(;
    edit_link="LF_Docs",
    size_threshold_ignore = ["man/listfunctions.md"]), # easier local build
    modules = [HPSAnalysis],
	remotes=nothing,
    pages = [
        "Home" => "index.md",
        "Setup" => "setup.md",
        "List of functions" => "listfunctions.md",
        #"Installation" => "man/installation.md",
        #"List of functions" => "man/listfunctions.md",
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    branch="LF_Docs"
)
