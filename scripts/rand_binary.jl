using Random

print("Enter n: ")
n = parse(Int, readline())

print("Enter d: ")
d = parse(Int, readline())

# Generate random dataset in [0, 10]
point = rand(n, d) .* 10.0

mkpath("datasets")

open("datasets/rand.bin", "w") do file
    write(file, Int64(d))
    write(file, Int64(n))

    for i in 1:n
        for j in 1:d
            write(file, point[i, j])
        end
    end
end

println("Dataset written to datasets/rand.bin")