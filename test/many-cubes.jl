using Plasm,Random


function randcubes(n)
	randomcubes = []
	for _ = 1:n
		f = rand(Float64,6)
		i1 = rand(1:3,1)
		i2 = rand(1:3,1)
		i3 = rand(1:3,1)
		i4 = rand(1:3,1)
		rcube = STRUCT(T(i1)(f[1]), R(1,2)(f[2]),  R(2,3)(f[3]),  R(1,3)(f[4]), CUBE(f[5]))
		push!(randomcubes, rcube)
	end
	return BOOL(UNION(randomcubes))
end

obj = randcubes(20)

VIEWCOMPLEX(obj, explode=[2,2,2])
VIEWCOMPLEX(BOUNDARY(obj), explode=[2,2,2])
VIEWCOMPLEX(BOUNDARY(obj))



global TEST_DIR = mk_output_dir()
cd(TEST_DIR)
