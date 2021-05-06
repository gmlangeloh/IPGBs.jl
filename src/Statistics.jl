module Statistics

export GBStats

abstract type GBStats end

function Base.show(
    io :: IO,
    stats :: T
) where {T <: GBStats}
    if hasfield(T, :title)
        println(io, stats.title)
    else
        println(io, "Algorithm statistics:")
    end
    i = 1
    fields = fieldnames(typeof(stats))
    num_fields = length(fields)
    for field in fields
        if field == :title
            i += 1
            continue
        elseif i < num_fields
            println(io, field, " => ", getfield(stats, field))
        else
            print(io, field, " => ", getfield(stats, field))
        end
        i += 1
    end
end

end
