using Plasm

# //////////////////////////////////////////////////////
function MyMain()
	lar=ToLARForm(STRUCT([
		CUBOID([1,1,1]),
		T([1])([3]),
		CUBOID([1,1,1])
	]))
	println("points ",lar.childs[1].points)
	println("hulls  ",lar.childs[1].hulls)
	println("facets ",lar.childs[1].facets)

end

MyMain()