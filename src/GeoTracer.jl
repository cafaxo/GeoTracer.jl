module GeoTracer

using LinearAlgebra
using StaticArrays

export
    Construction,
    new_free_point!,
    new_distance!,
    new_join!,
    new_intersect_line_line!,
    new_intersect_line_circle!,
    new_point_on_line!,
    new_circle_mp!,
    new_circle_mr!,
    new_intersect_circle_circle!,
    new_select_point!,
    new_point_on_circle!,
    move_to!

include("core.jl")

end # module
