
module StructDispatch
    SD_ChoosenFloatType=Float64
    SD_ChoosenIntType=Int64
    SD_StructNameToSymbols = Dict(Symbol("SimData") => Symbol("Sim"), Symbol("WorkData") => Symbol("Work"))
    SD_SymbolsTypesToReplace = [] #[Symbol("Sim"), Symbol("Work")]
    SD_PrecompileTypes  = Dict(Symbol("Sim")=>"SimData{$(SD_ChoosenFloatType), $(SD_ChoosenIntType)}", Symbol("Work")=>"WordData{$(SD_ChoosenIntType), $(SD_ChoosenIntType)}")
    SD_DataTypesToReplace = Dict{Tuple{Symbol, Symbol}, Bool}()

    macro GetType(Struct)
        ### get symbol of the Struct
        local struct_symbol = Struct.args[2].args[1]
        push!(SD_SymbolsTypesToReplace, SD_StructNameToSymbols[struct_symbol])

        ### get block expression
        local block = Struct.args[3].args

        for expression in block
            ### skip lineNumberNodes
            if typeof(expression) == Expr
                ### type definition or declaration?
                if expression.head == Symbol("::")
                    ### detectect a mutable type
                    if ~(typeof(expression.args[2]) == Symbol)
                        SD_DataTypesToReplace[(SD_StructNameToSymbols[struct_symbol], expression.args[1])] = true
                    end 
                end
            end
        end
        #println(Struct)
        #println(SD_DataTypesToReplace)
        return quote
            $(esc(Struct))
        end
    end

    macro apply(ex::Expr)
        local funcname  = ex.args[1].args[1]
        local DataTypesToReplace = SD_SymbolsTypesToReplace 
        local DataTypes = Vector{Symbol}()
        local PrecompileType = SD_PrecompileTypes 
        local ToUnroll = Vector{Set{Symbol}}()
        local ToCreate = Vector{Expr}() ### Why not set?
        local CreateFuncs = [Symbol("zeros"), Symbol("ones"), Symbol("OffsetArray"), Symbol("Dict")]
        local OtherFuncVars = Vector{Expr}()
        local Keywords = Vector{Any}()
        local Arguments = Vector{Symbol}()
        local ArgumentTypes = Vector{Any}()
        local AllDataTypesKnown = true
    
        if length(ex.args[1].args)>2
            for args in ex.args[1].args
                ### 
                if args in DataTypesToReplace
                    continue
                end
                if typeof(args)==Symbol && args != funcname
                    AllDataTypesKnown = false
                end
                if (typeof(args)==Expr)
                    if(args.head==Symbol("parameters"))
                        Keywords= args.args
                    elseif (args.head==Symbol("::"))
                        push!(Arguments, args.args[1])
                        push!(ArgumentTypes, args.args[2])
                    else
                        if length(args.args)<2
                            AllDataTypesKnown=false
                        end
                        #dump(args)
                        push!(OtherFuncVars, args)
                    end
                end
            end
        end
    
        for (ind,DataType) in enumerate(DataTypesToReplace)
            if DataType in ex.args[1].args
                push!(DataTypes, DataType)
            end
        end
    
    
    
        ### Replace Struct inside function body
        ### OPTIMIZE: replace stuff in one substitute iteration
        for (ind, DataType) in enumerate(DataTypes)
            push!(ToUnroll, Set{Symbol}())
            substitute(ex.args[2], DataType, ToUnroll[ind], ToCreate, CreateFuncs)
        end
    
        ### Write new function head
        local str="function $(funcname)("
        for (ind, DataType) in enumerate(DataTypes)
            str *= "$(DataType),"
        end
        
        for var in OtherFuncVars
            ### catch non mutable data types
            if length(var.args)>=2 && typeof(var.args[2])!=Symbol
                str *="$(var.args[1])::$(typeof(var.args[2]))=$(var.args[2]),"
            else
                str*="$(var),"
            end
        end
        
        for (ind, argument) in enumerate(Arguments)
            if ArgumentTypes[ind]!= nothing
                str *= "$argument::$(ArgumentTypes[ind]), "
            else
                precompile *= "$(argument)"
            end
        end

        str=Base.chopsuffix(str,",")
    
        if length(Keywords)>0
            str*=";"
            for keyword in Keywords
                if length(keyword.args)>=2
                    str *="$(keyword.args[1])::$(typeof(keyword.args[2]))=$(keyword.args[2]),"
                else
                    str *="$(keyword.args[1]),"
                    AllDataTypesKnown = false
                end
            end
            str=Base.chopsuffix(str,",")

        end
        str*=")\n"
    
    
        ### include allocations to function body
        for create in ToCreate
            str*= chop(repr(create);head=2,tail=1)*"\n"
        end
        
        ### call the the new function and return it
        str*="return $(funcname)_SD("
        for DataType in DataTypes
            str*="$(DataType), "
        end

        for (ind, DataType) in enumerate(DataTypes)
            ### collect converts Set to vector so that it can be iterated in reverse, needed later to push elements of Unroll to func definition
            for field in reverse(collect(ToUnroll[ind])) 
                str*="$(DataType).$(field),"
            end
        end

        for var in OtherFuncVars
            ### catch non mutable data types
            if length(var.args)>=2 && typeof(var.args[2])!=Symbol
                str *="$(var.args[1]),"#::$(typeof(var.args[2])),"
            else
                str*="$(var.args[1]),"
            end    
        end
        for (ind, argument) in enumerate(Arguments)
            if ArgumentTypes[ind]!= nothing
                str *= "$argument::$(ArgumentTypes[ind]), "
            else
                precompile *= "$(argument)"
            end
        end
    
        ### pass keyword arguments by name to new function
        if length(Keywords)>0
            str=Base.chopsuffix(str,",")
            str*=";"
            for keyword in Keywords
                str *="$(keyword.args[1])=$(keyword.args[1]),"
            end
        end
        str=Base.chopsuffix(str,",")
        str*=")\n end \n"
    
        new_func = Meta.parse(str)
    
        ### Rename the function in the initial code 
        ex.args[1].args[1]=Symbol("$(funcname)_SD")
        local NDataTypes = length(DataTypes)+2
        for (ind, DataType) in enumerate(DataTypes)
            for field in ToUnroll[ind]
                insert!(ex.args[1].args, NDataTypes, Symbol("$(DataType)_$(field)"))
            end
        end
    
        ### create precompile statement
        precompile=""
        if AllDataTypesKnown
            precompile= "precompile($(funcname), ("
            for DataType in DataTypes
                precompile*= "$(PrecompileType[DataType]),"
            end 
            for var in OtherFuncVars
                if length(var.args)>=2 && typeof(var.args[2])!=Symbol
                    precompile *="$(typeof(var.args[2])),"
                else
                    precompile *="$(var.args[2]),"
                end
            end
            for (ind, argument) in enumerate(Arguments)
                if ArgumentTypes[ind]!= nothing
                    precompile *= "$(ArgumentTypes[ind]),"
                #else
                #    precompile *= "$(argument)"
                end
            end
            for keyword in Keywords
                precompile *="$(typeof(keyword.args[2])),"
            end
            precompile*= " ) )"
        #else
        #    println("Precompile failed: $(funcname)")
        end
        
        if funcname==Symbol("computeSlabHistogram") && true
            #println(str)
            println(new_func)
            println(ex)
            #end
            println(precompile)
           # println(OtherFuncVars)
           # println(Arguments)
          #  println(String.(ArgumentTypes))
        end
    
        return quote 
            $(esc(new_func))
            $(esc(ex))
            $(esc(Meta.parse(precompile)))
        end
    end

    function substitute(expr::Expr, ToSub , ToUnroll, ToCreate, createFuncs)
        local todelete = []
        for i in 1:length(expr.args)
            arg = expr.args[i]
            if typeof(arg)== Expr
                if (arg.head)===Symbol("=")
                    if length(arg.args)>=2 
                        if(typeof(arg.args[1])==Expr && typeof(arg.args[2])==Expr)
                        ### dont parse further if the assignment is not to array type, those types would not be changed in the function scope if they are not part of a mutable struct
                        if(arg.args[1].args[1]==ToSub)
                            if(~haskey(SD_DataTypesToReplace, (arg.args[1].args[1], arg.args[1].args[2].value)))
                                continue
                            end
                        end

                        ### check if its struct field and a call assignment
                        if arg.args[1].head===Symbol(".") && arg.args[2].head===Symbol("call")
                            ### Check whether to be replaced field is part of the assignment
                            if(arg.args[1].args[1]==ToSub && typeof(arg.args[1].args[2])==QuoteNode)
                               
                                ### search for functions that allocate memory
                                if(arg.args[2].args[1] in createFuncs)
                                    push!(ToCreate, copy(arg))
                                    push!(todelete, i)
                                end
                            end
                        end 
                        end
                    end
                end
                if length(arg.args)>=2 
                    ### Check whether symbol matches and if the member variable is stored as a QuoteNode
                    if arg.args[1]===ToSub && typeof(arg.args[2])==QuoteNode
                        #####################
                        fieldSymbol = arg.args[2].value ### Storage of the Symbol for a Quote node
                        push!(ToUnroll,fieldSymbol) 
                        expr.args[i] = Symbol("$(ToSub)_$(fieldSymbol)") ### Replace Expression by Symbol, so that no fields are accessed
                        continue
                    end
                end
                substitute(arg, ToSub, ToUnroll, ToCreate,createFuncs)
            else
                continue
            end
        end
        for i in sort(todelete, rev=true)
            deleteat!(expr.args,i)
        end
        return 
    end
end