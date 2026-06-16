module face_vertex_mismatch

function face_vertex_mismatch_report(lar)
    EV = lar.C[:EV]
    FE = lar.C[:FE]
    FV = lar.C[:FV]

    println("Checking faces...")
    for f in eachindex(FE)
        fe = FE[f]
        fv = FV[f]

        usedverts = unique(vcat([EV[e] for e in fe]...))

        if length(fv) != length(unique(fv)) ||
           Set(fv) != Set(usedverts) ||
           length(fv) != length(fe)

            println("\nface ", f)
            println("  FE = ", fe)
            println("  |FE| = ", length(fe))
            println("  FV = ", fv)
            println("  |FV| = ", length(fv))
            println("  unique(FV) = ", unique(fv))
            println("  |unique(FV)| = ", length(unique(fv)))
            println("  vertices from FE = ", sort(usedverts))
            println("  |vertices from FE| = ", length(usedverts))
        end
    end
end

end # module