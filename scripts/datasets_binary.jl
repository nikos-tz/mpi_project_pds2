"""
Write a matrix into a binary file (to read from C).
With this simple script we download and read the MNIST dataset, 
you can find other datasets here: https://juliaml.github.io/MLDatasets.jl/stable/
If you want to read other datasets you will have to change the values of the variables:
train_x, test_x, file_path and sometimes the dimension.
"""

using LinearAlgebra, MLDatasets

function write_matrix_binary(A::AbstractMatrix{Float64}, filename::String)

	open(filename, "w") do file
    	write(file, size(A,1)) # write an Int64
    	write(file, size(A,2)) # write an Int64
		for i in eachindex(A)
			write( file, A[i] ) # write Float64 (column-major by default)
		end
    end
end

file_path = "/tmp/mnist.bin";
dimension = 28;

train_x, _ = MNIST.traindata()
test_x,  _ = MNIST.testdata()
X = Float64.( [reshape(train_x,dimension*dimension,:) reshape(test_x,dimension*dimension,:)] )
write_matrix_binary(X, file_path)