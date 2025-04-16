function getPythonEnvironment(PkgSourcePath)
    EnvironmentPath=""
    file="$(PkgSourcePath)/../data/EnvironmentPath.txt" 
    println(file)
    if isfile(file)
        f = open(file,"r")
        EnvironmentPath = readline(f)
        println(EnvironmentPath)
    end
    if EnvironmentPath==""
        @warn("Python Environment has not been specified. All functionalty regarding Polyply and HOOMD may be breaking or omitted.")
    end
    return EnvironmentPath
end