using DataFrames, CSV
U = [1;2;3;4;;5;6;7;8;;;
9;10;11;12;;13;14;15;16;;;
17;18;19;20;;21;22;23;24]
X, Y, Z = size(U)
println("Sizes: ", X, ", ", Y, ", ", Z)

df = DataFrame(
    i = repeat(1:X, outer=Y*Z),
    j = repeat(repeat(1:Y, inner=X), outer=Z),
    k = repeat(1:Z, inner=X*Y),
    concentration = vec(U)
)
CSV.write("test.csv", df)