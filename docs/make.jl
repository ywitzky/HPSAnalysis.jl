using Documenter

push!(LOAD_PATH,"../src/")
using HPSAnalysis



@info "Making documentation..."


makedocs(
    sitename = "HPSAnalysis",
    authors = "Yannick Witzky",
<<<<<<< HEAD
    format = Documenter.HTML(;
    edit_link="LF_Docs",
    size_threshold_ignore = ["man/listfunctions.md"]), # easier local build
=======
    format = Documenter.HTML(;edit_link=:commit,
    prettyurls = get(ENV, "CI", nothing) == "true"), # easier local build
>>>>>>> 930e2ad244d91538085bd75b946dcd24142ff99c
    modules = [HPSAnalysis],
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    pages = [
        "Home" => "index.md",
        "Setup" => "setup.md",
<<<<<<< HEAD
        "List of functions" => "listfunctions.md",
        #"Installation" => "man/installation.md",
        #"List of functions" => "man/listfunctions.md",
=======
>>>>>>> 930e2ad244d91538085bd75b946dcd24142ff99c
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
<<<<<<< HEAD
deploydocs(
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    branch="LF_Docs"
)
=======
#=deploydocs(
    devbranch="LF_Docs",
    repo = "https://gitlab.rlp.net/ywitzky/HPSAnalysis.jl",
    push_preview=true,
)=#
>>>>>>> 930e2ad244d91538085bd75b946dcd24142ff99c
