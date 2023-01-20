using JLD2
using Plots
using DelimitedFiles

function ReadData(fileName)
    fileNamesWithVector = ["halfLN", "LN", "smallLN", "corrlN"]
    if ContainsMeanAndVariance(fileName) && !(ContainsProperFluc(fileName))
        data = load("/data/htolppan/Data for the paper/" * fileName * "/data.jld2", "x", "y")
        return data[1], SeparateMeanAndVar(data[2][1]), SeparateMeanAndVar(data[2][2]), SeparateMeanAndVar(data[2][3]), SeparateMeanAndVar(data[2][4])
    else
        data = load("/data/htolppan/Data for the paper/" * fileName * "/data.jld2", "x", "y")
        return data[1], data[2][1], data[2][2], data[2][3], data[2][4]
    end
end
function ContainsProperFluc(filename)
    return chop(filename, head = 3, tail = 0) == "Proper_fluctuation"
end
function ContainsMeanAndVariance(fileName)
    fileNamesWithVector = ["halfLN", "LN", "smallLN", "corrlN"]
    return !(chop(fileName, head = 3, tail = 0) in fileNamesWithVector)
end

function SeparateMeanAndVar(data)
    return [res[1] for res in data], [res[2] for res in data]
end

function SaveForMatlab(fileName)
    x, y4, y6, y8, y10 = ReadData(fileName)
    SaveVectorMatlab(x, fileName, "probabilities")
    if ContainsProperFluc(fileName)
        SaveVectorMatlab(y4, fileName, "L4")
        SaveVectorMatlab(y6, fileName, "L6")
        SaveVectorMatlab(y8, fileName, "L8")
        SaveVectorMatlab(y10, fileName, "L10")
    elseif ContainsMeanAndVariance(fileName)
        SaveVectorMatlab(y4[1], fileName, "L4_Mean")
        SaveVectorMatlab(y4[2], fileName, "L4_Var")
        SaveVectorMatlab(y6[1], fileName, "L6_Mean")
        SaveVectorMatlab(y6[2], fileName, "L6_Var")
        SaveVectorMatlab(y8[1], fileName, "L8_Mean")
        SaveVectorMatlab(y8[2], fileName, "L8_Var")
        SaveVectorMatlab(y10[1], fileName, "L10_Mean")
        SaveVectorMatlab(y10[2], fileName, "L10_Var")
    else
        SaveVectorMatlab(reduce(hcat, y4), fileName, "L4")
        SaveVectorMatlab(reduce(hcat, y6), fileName, "L6")
        SaveVectorMatlab(reduce(hcat, y8), fileName, "L8")
        SaveVectorMatlab(reduce(hcat, y10), fileName, "L10")
    end
end

function SaveVectorMatlab(data, observableName, dataTypeName)
    path = pwd() * "/Data Matlab/" * observableName
	mkpath(path)
    writedlm(path * "/" * dataTypeName * ".csv", data)
end

function f()
    x, y4, y6, y8, y10 = ReadData("M2_Entanglement")
    display(length(y8[1]))
    plot(x, y8[1])
end

f()