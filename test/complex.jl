using Test, Plasm, LinearAlgebra

import Base.+  
+(f::Function, g::Function) = (x...) -> f(x...) + g(x...)  

# //////////////////////////////////////////////////////////////////////////////
@testset "IsSimplex" begin
   @test IsPolytope(CUBE(100))
   @test IsSimplex(CUBOID([1]))  == true
   @test IsSimplex(SIMPLEX(1))  == true
   @test IsSimplex(CUBOID([1,1])) == false
   @test IsPolytope(CUBOID([1,1]))  == true
   @test IsPolytope(CUBOID([3,3]))  == true
   @test IsPolytope(CUBE(100))
end

# //////////////////////////////////////////////////////////////////////////////
@testset "SIMPLEX" begin
   @test SIMPLEX(1).T == MatrixNd(2)
   @test SIMPLEX(4).T == MatrixNd(5)
   @test SIMPLEX(20).childs[1].hulls == [1:21]
   @test SIMPLEX(20).childs[1].points[1] == zeros(20)
end

simplexfacets

# //////////////////////////////////////////////////////////////////////////////
@testset "CUBOID" begin
   @test typeof(CUBOID([1]).childs[1].childs[1].points) == Vector{Vector{Float64}}
   @test typeof(CUBOID([1]).childs[1].childs[1].hulls) == Vector{Vector{Int64}}
   @test CUBOID([1]).childs[1].childs[1].points == [[0.0],[1.0]]
   @test CUBOID([1]).childs[1].childs[1].hulls == [[1,2]]
   @test typeof(CUBOID([1,1,1,1]).childs[1].childs[1].points[1]) == Vector{Float64}
   @test typeof(CUBOID([1,1,1,1]).childs[1].childs[1].hulls[1]) == Vector{Int64}
   @test CUBOID(DIESIS(4)(1)).childs[1].childs[1].points[1] == [0.,0.,0.,0.]
   @test CUBOID(DIESIS(4)(1)).childs[1].childs[1].hulls[1] == 1:2^4
end;
@testset "CUBOID" "$shape" for shape in [[1],[3,2,1],[3,2],[1,1,1,1]]
   @test CUBOID(shape).childs[1].childs[1].points[1] == zeros(LEN(shape))
end; 

# //////////////////////////////////////////////////////////////////////////////
# scrgiorgio 20240117 this fails on Windows, so commented
#@testset "CHULL1" begin
# points = rand(6, 3)
#  obj = CHULL(points)
#  VIEW(Hpc(obj.V, obj.C[:EV]))
#end

# 20240117 scrgiorgio this fails on Windows, so commented
#@testset "CHULL2" begin
#  points = [0.7847031769013115 0.2583051472567749 0.2807219719987384; 0.5873015490214822 0.27654940664002525 0.09814339221099866; 0.312439556286509 0.3509835618138125 0.8606532158160981; 0.7999175319797596 0.48136552088275075 0.7250168308949106; 0.13440321955041745 0.9288317959331052 0.9330767543972907; 0.6934214563208536 0.4598639813793165 0.4985288983965731; 0.7957846859216022 0.1818042912907808 0.827158332120763; 0.8049974871790887 0.052020575464559626 0.3787205333768807; 0.11906445712905878 0.5758185746980895 0.9436673428232172; 0.742006419372256 0.14654659755790478 0.12596494694361882; 0.6660716740819993 0.6113128167395886 0.3929073845909179; 0.9135424489502044 0.31765286438760365 0.43453539121025553; 0.7781306720402988 0.2573100295147508 0.563954615698295; 0.05528182093812595 0.29015402634285503 0.298969626011406; 0.3267602613807892 0.854950337141562 0.06678416929022801; 0.022879218606704943 0.2777397179009111 0.46421071837052685; 0.8371228951517155 0.05514197476407556 0.8270188321860387; 0.3658455142208985 0.3126874269600419 0.7491840955501928; 0.7108748785929014 0.4086649443194168 0.7776370327406069; 0.7559114508004247 0.258350210798748 0.43571974025035654; 0.9116925963537305 0.2212404823266687 0.810970314757372; 0.5813429803777005 0.976162078802727 0.8498145982227862; 0.322627031463959 0.2526932738886688 0.7366817698037228; 0.18603910250206146 0.8539683831859929 0.14047668896152443; 0.5370794040559644 0.6389189387470219 0.9601323882002971; 0.43922613689329293 0.4268379347298201 0.6815530406735171; 0.011770873888801825 0.13329169190086454 0.6283773138712718; 0.27512751472348806 0.4351195985505101 0.8629721808125018; 0.4787794115149303 0.10190899229160122 0.3667944870403044; 0.13232502529863666 0.49398862907201657 0.38432554140905373; 0.1895733150898662 0.8942646502716064 0.7802325605123039; 0.5989601520031026 0.38655195808759935 0.08264133193402468; 0.7342284420078041 0.5261343060599059 0.99317412532261; 0.4672197356509594 0.18167514692420927 0.8991414610494608; 0.3065975848591662 0.42092197047670943 0.2673380426213542; 0.6237711099002246 0.3740701837032857 0.729743218114104; 0.32163567216171374 0.7189394838739409 0.39482269379038326; 0.47251455932447406 0.9519905755093248 0.860237213579729; 0.28061632557070604 0.5613764375614576 0.373183792703193; 0.15942853298022808 0.5826594476501387 0.5831925356433824; 0.3175057543064439 0.41491477439649616 0.8594538296869134; 0.3651769483519681 0.04550702057370326 0.04949923653032384; 0.25186155691437584 0.7950934998787013 0.5120664555074361; 0.9436629965666474 0.0502161197810822 0.6034181095111704; 0.36326824406066693 0.9695824865819831 0.5927792897939886; 0.5478951988518589 0.8505705866157713 0.6890423473828612; 0.8375075231283257 0.15171621123909096 0.337639705506534; 0.4767630629330035 0.31075103448631813 0.6774852422923894; 0.2803529656767 0.37315054566284 0.842027959128257; 0.9193874522838082 0.6562944639235494 0.9187606880834863]
#  obj = CHULL(points)
#  @test Set(union(vcat(obj.C[:EV]...))) == Set(union(vcat(obj.C[:FV]...)))
#  VIEW(Hpc(obj.V, obj.C[:EV]))
#end

# //////////////////////////////////////////////////////////////////////////////
@testset "VIEWCOMPLEX"  begin  # View 2D and 3D Lar
   mesh2D = CUBOIDGRID([5,5])
   VIEWCOMPLEX(mesh2D, properties=Dict())
   mesh3D = CUBOIDGRID([2,3,4])
   VIEWCOMPLEX(mesh3D, properties=Dict("background-color"=>Point4d(1,1,1,1)))
end;


# //////////////////////////////////////////////////////////////////////////////
# extracted from alberto's code

"""
SIMPLEX(1)
SIMPLEX(2)
SIMPLEX(3)
SIMPLEX(4)

CUBOID([1,1])
CUBOID([1,1,1])
CUBOID([1,1,1,1])

simplex6 = simplex(6,complex=true)
simplex6.V
simplex6.C
Hpc(simplex6.V, simplex6.C[:C5V])

Simplex6 = Lar(Hpc(simplex6.V, simplex6.C[:C5V]))

obj3 = Hpc(Lar(CUBOID([1,1,1]))); # OK
VIEW(Hpc(Lar(obj3).V, Lar(obj3).C[:EV])) # OK
VIEW(Hpc(Lar(obj3).V, Lar(obj3).C[:FV])) # OK

LAR(STRUCT(SQUARE(1), T(1,2)(0.5,0.25), SQUARE(1))) # 2D
LAR(STRUCT(CUBE(1), T(1,2)(0.5,0.25), CUBE(1))) # 3D |V|=16,|FV|=12,|EV|=24 ... OK!
LAR(STRUCT(CUBE(1), T(1)(1), CUBE(1))) # |V|=12,|FV|=11,|EV|=20 ... OK!  correct by using LAR ... uncorrect using Lar !!!

cube3 = (LAR∘STRUCT∘CAT∘DIESIS(3))([T(1)(1), CUBE(1)]); # NB:  LAR  !!!
VIEW(Hpc(cube3.V, cube3.C[:EV])) # OK
VIEW(Hpc(cube3.V, cube3.C[:FV])) # OK

grid1D = INTERVALS(10.0)(10) # equivalent to the next
grid1D = QUOTE(DIESIS(10)(1.0)) # equivalent to the former

mesh = CUBOIDGRID([5,5,5])
VIEW(Hpc(mesh.V, mesh.C[:EV]))

GRID1(5)

mesh = CUBOIDGRID([10,10])
VIEW(Hpc(mesh.V, mesh.C[:EV])) # OK ai bordi
VIEW(Hpc(mesh.V, mesh.C[:FV])) # OK ai bordi
"""







