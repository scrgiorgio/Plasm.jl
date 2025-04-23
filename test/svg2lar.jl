using Plasm

function myshow(filename)
	pol = SVG(filename)
	VIEWCOMPLEX(LAR(pol))
	return pol
end

pol = myshow("./svg/new.svg");
pol = myshow("svg/curved.svg");
pol = myshow("svg/twopaths.svg");
pol = myshow("svg/paths.svg");
pol = myshow("svg/boundarytest2.svg");
pol = myshow("svg/tile.svg");
pol = myshow("svg/interior.svg");
pol = myshow("svg/holes.svg");
VIEWCOMPLEX(ARRANGE2D(LAR(pol)));
pol = myshow("svg/Lar.svg");
