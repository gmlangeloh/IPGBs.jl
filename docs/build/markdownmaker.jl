path = "../../src/"
files = readdir(path)

exceptions = ["Globals"]

for filename in files
    if filename in exceptions
        continue
    end
    parts = split(filename, ".")
    name = parts[1]
    mdfilename = name * ".md"
    module_name = ""
    if name == "IPGBs"
        module_name = "IPGBs"
    else
        module_name = "IPGBs." * name
    end
    f = open(mdfilename, "w")
    println(f, "# ", filename)
    println(f, "Documentation for ", filename)
    println(f)
    println(f, """
```@autodocs
Modules = [$module_name]
```
    """)
    close(f)
end
