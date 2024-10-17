using Test


@testset "Integration Tests" begin
	@testset "M" begin
		@test M(0,0)==0.5
		@test M(1,0)==0.16666666666666666
		@test M(2,0)==0.08333333333333333
		@test M(3,0)==0.05
		@test M(1,1)==0.041666666666666685
		@test M(2,0)==0.08333333333333333
		@test M(2,1)==0.016666666666666663
		@test M(3,0)==0.05
		@test M(3,1)==0.008333333333333338
		@test M(3,2)==0.0023809523809523586
		@test M(3,3)==0.0008928571428571397
	end


	@testset "TTT" begin
		tau=[0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		@test TTT(tau, 0,0,0)==0.5
		@test TTT(tau, 1,0,0)==0.16666666666666666
		@test TTT(tau, 1,1,0)==0.041666666666666685
		@test TTT(tau, 1,1,1)==0.0
		@test TTT(tau, 2,0,0)==0.08333333333333333
		@test TTT(tau, 2,1,0)==0.016666666666666663
		@test TTT(tau, 2,2,0)==0.005555555555555545
		@test TTT(tau, 2,2,1)==0.0
		@test TTT(tau, 2,2,2)==0.0
	end


	@testset "II" begin
		V = [0.0 1.0 0.0; 0.0 0.0 1.0; 0.0 0.0 0.0];
		FV = [[1,2,3]];
		P = V,FV;
		@test II(P, 0,0,0)==0.5
		@test II(P, 1,0,0)==0.16666666666666666
		@test II(P, 0,1,0)>=0.1666666666666666
		@test II(P, 0,0,1)==0.0
		@test II(P, 1,1,1)==0.0
	end


	@testset "III" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test III(P, 0,0,0)>0.166666666
		@test III(P, 0,0,0)<0.166666888
		@test III(P, 1,0,0)>0.041666666
		@test III(P, 1,0,0)<0.041666888
		@test III(P, 0,1,0)>0.041666666
		@test III(P, 0,1,0)<0.041666888
		@test III(P, 0,0,1)>0.041666666
		@test III(P, 0,0,1)<0.041666888
		@test III(P, 10,10,10)>1.3377e-11
		@test III(P, 10,10,10)<1.3388e-11
	end


	@testset "SURFACE" begin
		V,FV = SIMPLEXGRID([1,1]);
		P = [V;[0 0 0 0]], FV
		@test SURFACE(P)==1.0
		p = LAR(STRUCT( T(1,2)(0.5,0.5), R(2,3,)(pi/4), MKPOL(P)));
      q = p.V, p.C[:FV]; 
      @test SURFACE(q)>1.0000000
		@test SURFACE(q)<1.0000222
	end


	@testset "VOLUME" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test VOLUME(P)>0.166666666
		@test VOLUME(P)<0.166668888
	end


	@testset "FIRSTMOMENT" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test FIRSTMOMENT(P)[1]<0.0416667
		@test FIRSTMOMENT(P)[1]>0.0416665

		@test FIRSTMOMENT(P)[2]<0.0416667
		@test FIRSTMOMENT(P)[2]>0.0416665

		@test FIRSTMOMENT(P)[3]<0.0416667
		@test FIRSTMOMENT(P)[3]>0.0416665

		@test abs(FIRSTMOMENT(P)[1]-FIRSTMOMENT(P)[2])<0.00001
		@test abs(FIRSTMOMENT(P)[2]-FIRSTMOMENT(P)[3])<0.00001
		@test abs(FIRSTMOMENT(P)[3]-FIRSTMOMENT(P)[1])<0.00001
	end


	@testset "SECONDMOMENT" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test SECONDMOMENT(P)[1]<0.0166666669
		@test SECONDMOMENT(P)[1]>0.0166666664

		@test SECONDMOMENT(P)[2]<0.0166666669
		@test SECONDMOMENT(P)[2]>0.0166666664

		@test SECONDMOMENT(P)[3]<0.0166666669
		@test SECONDMOMENT(P)[3]>0.0166666664

		@test abs(SECONDMOMENT(P)[1]-SECONDMOMENT(P)[2])<0.00001
		@test abs(SECONDMOMENT(P)[2]-SECONDMOMENT(P)[3])<0.00001
		@test abs(SECONDMOMENT(P)[3]-SECONDMOMENT(P)[1])<0.00001
	end


	@testset "INERTIAPRODUCT" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test INERTIAPRODUCT(P)[1]<0.00833666
		@test INERTIAPRODUCT(P)[1]>0.00833000

		@test INERTIAPRODUCT(P)[2]<0.00833666
		@test INERTIAPRODUCT(P)[2]>0.00833000

		@test INERTIAPRODUCT(P)[3]<0.00833666
		@test INERTIAPRODUCT(P)[3]>0.00833000
	end


	@testset "CENTROID" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test CENTROID(P)[1]<0.26
		@test CENTROID(P)[1]>0.24

		@test CENTROID(P)[2]<0.26
		@test CENTROID(P)[2]>0.24

		@test CENTROID(P)[3]<0.26
		@test CENTROID(P)[3]>0.24
	end


	@testset "INERTIAMOMENT" begin
		V = [0.0 1.0 0.0 0.0; 0.0 0.0 1.0 0.0; 0.0 0.0 0.0 1.0];
		FV = [[1, 2, 4], [1, 3, 2], [4, 3, 1], [2, 3, 4]];
		P = V,FV;
		@test INERTIAMOMENT(P)[1]<0.0333555
		@test INERTIAMOMENT(P)[1]>0.0333111

		@test INERTIAMOMENT(P)[2]<0.0333555
		@test INERTIAMOMENT(P)[2]>0.0333111

		@test INERTIAMOMENT(P)[3]<0.0333555
		@test INERTIAMOMENT(P)[3]>0.0333111
	end

end

