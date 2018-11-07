import Base: getindex

include("point.jl")

abstract type GeoObject end

struct Construction{T <: Number}
    objects::Vector{GeoObject}
end

Construction{T}() where T = Construction{T}(Vector{GeoObject}())

struct GeoPointer{T}
    index::Int
end

function getindex(construction::Construction, pointer::GeoPointer{T}) where T
    return state(construction.objects[pointer.index])::T
end

struct FreePoint{T} <: GeoObject
    position::Point{T}
end

state(free::FreePoint) = free.position

function new_free_point!(construction::Construction{T}, position::Point) where T
    free = FreePoint(convert(Point{T}, position))
    push!(construction.objects, free)

    return GeoPointer{Point{T}}(length(construction.objects))
end

update(::Construction, p::FreePoint) = p

update_param(::Construction, ::FreePoint, p::Point) = FreePoint(p)

function interpolate(::Construction, t::Real, point::FreePoint, target::Point)
    t2 = t * t
    dt = 0.5 / (1 + t2)
    tc = complex((2 * t) * dt + 0.5, (1 - t2) * dt)

    source = point.position

    return source + tc * (target - source)
end

struct Distance{T,S <: Real} <: GeoObject
    a::GeoPointer{Point{T}}
    b::GeoPointer{Point{T}}

    d::S
end

state(dist::Distance) = dist.d

function new_distance!(construction::Construction{T}, a::GeoPointer{Point{T}}, b::GeoPointer{Point{T}}) where T
    dist = Distance(a, b, distance(construction[a], construction[b]))

    push!(construction.objects, dist)

    return GeoPointer{real(T)}(length(construction.objects))
end

update(construction::Construction, dist::Distance) = Distance(dist.a, dist.b, distance(construction[dist.a], construction[dist.b]))

struct Join{T} <: GeoObject
    a::GeoPointer{Point{T}}
    b::GeoPointer{Point{T}}

    line::Point{T}
end

state(j::Join) = j.line

function new_join!(construction::Construction{T}, a::GeoPointer{Point{T}}, b::GeoPointer{Point{T}}) where T
    join = Join(a, b, cross(construction[a], construction[b]))
    push!(construction.objects, join)

    return GeoPointer{Point{T}}(length(construction.objects))
end

update(construction::Construction, j::Join) = Join(j.a, j.b, cross(construction[j.a], construction[j.b]))

struct IntersectionLineLine{T} <: GeoObject
    l1::GeoPointer{Point{T}}
    l2::GeoPointer{Point{T}}

    point::Point{T}
end

state(ll::IntersectionLineLine) = ll.point

function new_intersect_line_line!(construction::Construction{T}, l1::GeoPointer{Point{T}}, l2::GeoPointer{Point{T}}) where T
    ll = IntersectionLineLine(l1, l2, cross(construction[l1], construction[l2]))

    push!(construction.objects, ll)

    return GeoPointer{Point{T}}(length(construction.objects))
end

function update(construction::Construction, ll::IntersectionLineLine)
    return IntersectionLineLine(ll.l1, ll.l2, cross(construction[ll.l1], construction[ll.l2]))
end

struct IntersectionLineCircle{T} <: GeoObject
    line::GeoPointer{Point{T}}
    circle::GeoPointer{CircleMatrix{T}}

    p1::Point{T}
    p2::Point{T}
end

state(lc::IntersectionLineCircle) = (lc.p1, lc.p2)

function new_intersect_line_circle!(construction::Construction{T}, line::GeoPointer{Point{T}}, circle::GeoPointer{CircleMatrix{T}}) where T
    lc = IntersectionLineCircle(line, circle, intersect_line_circle(construction[line], construction[circle])...)

    push!(construction.objects, lc)

    return GeoPointer{Tuple{Point{T},Point{T}}}(length(construction.objects))
end

function update(construction::Construction, lc::IntersectionLineCircle)
    result = tracing2((lc.p1, lc.p2), intersect_line_circle(construction[lc.line], construction[lc.circle]))

    if result === nothing
        return nothing
    end

    return IntersectionLineCircle(lc.line, lc.circle, result...)
end

struct PointOnLine{T} <: GeoObject
    line::GeoPointer{Point{T}}

    point::Point{T}
end

state(pol::PointOnLine) = pol.point

function new_point_on_line!(construction::Construction{T}, line::GeoPointer{Point{T}}, nearest::Point) where T
    nearest = convert(Point{T}, nearest)

    pol = PointOnLine(line, project_point_to_line(nearest, construction[line]))
    push!(construction.objects, pol)

    return GeoPointer{Point{T}}(length(construction.objects))
end

function update(construction::Construction, pol::PointOnLine, orig_construction::Construction, orig_pol::PointOnLine)
    new_line = construction[pol.line]
    old_point = pol.point

    real_line = orig_construction[orig_pol.line]
    real_point = orig_pol.point

    center = cross(real_line, new_line)

    # TODO: handle line stay almost the same
    @assert abs(center) > 3eps(Float64)

    circle = circle_from_mp(center, real_point)

    candidates = intersect_line_circle(new_line, circle)

    old_antipode = point_reflection(center, old_point)

    result = tracing2((old_point, old_antipode), candidates)

    if result === nothing
        return nothing
    end

    return PointOnLine(pol.line, result[1])
end

function update_param(construction::Construction, pol::PointOnLine, param::Point)
    return PointOnLine(pol.line, param)
end

function interpolate(construction::Construction, t::Real, pol::PointOnLine, target::Point)
    t2 = t * t
    dt = 0.5 / (1 + t2)
    tc = complex((2 * t) * dt + 0.5, (1 - t2) * dt)

    target = project_point_to_line(target, construction[pol.line])
    source = pol.point

    return source + tc * (target - source)
end

struct CircleMP{T} <: GeoObject
    m::GeoPointer{Point{T}}
    p::GeoPointer{Point{T}}

    matrix::CircleMatrix{T}
end

state(circle::CircleMP) = circle.matrix

function new_circle_mp!(construction::Construction{T}, m::GeoPointer{Point{T}}, p::GeoPointer{Point{T}}) where T
    circle = CircleMP(m, p, circle_from_mp(construction[m], construction[p]))
    push!(construction.objects, circle)

    return GeoPointer{CircleMatrix{T}}(length(construction.objects))
end

update(construction::Construction, c::CircleMP) = CircleMP(c.m, c.p, circle_from_mp(construction[c.m], construction[c.p]))

struct CircleMr{T,S <: Real} <: GeoObject
    m::GeoPointer{Point{T}}
    r::GeoPointer{S}

    matrix::CircleMatrix{T}
end

state(circle::CircleMr) = circle.matrix

function new_circle_mr!(construction::Construction{T}, m::GeoPointer{Point{T}}, r::GeoPointer{S}) where {T,S <: Real}
    @assert S == real(T)

    circle = CircleMr(m, r, circle_from_mr(construction[m], construction[r]))
    push!(construction.objects, circle)

    return GeoPointer{CircleMatrix{T}}(length(construction.objects))
end

update(construction::Construction, c::CircleMr) = CircleMr(c.m, c.r, circle_from_mr(construction[c.m], construction[c.r]))

struct IntersectionCircleCircle{T} <: GeoObject
    c1::GeoPointer{CircleMatrix{T}}
    c2::GeoPointer{CircleMatrix{T}}

    p1::Point{T}
    p2::Point{T}
end

state(cc::IntersectionCircleCircle) = (cc.p1, cc.p2)

function new_intersect_circle_circle!(construction::Construction{T}, c1::GeoPointer{CircleMatrix{T}}, c2::GeoPointer{CircleMatrix{T}}) where T
    cc = IntersectionCircleCircle(c1, c2, intersect_circle_circle(construction[c1], construction[c2])...)
    push!(construction.objects, cc)

    return GeoPointer{Tuple{Point{T},Point{T}}}(length(construction.objects))
end

function update(construction::Construction, cc::IntersectionCircleCircle)
    result = tracing2((cc.p1, cc.p2), intersect_circle_circle(construction[cc.c1], construction[cc.c2]))

    if result === nothing
        return nothing
    end

    return IntersectionCircleCircle(cc.c1, cc.c2, result...)
end

struct SelectPoint{T} <: GeoObject
    collection::GeoPointer{Tuple{Point{T},Point{T}}}

    index::Int
    p::Point{T}
end

state(select::SelectPoint) = select.p

function new_select_point!(construction::Construction{T}, collection::GeoPointer{Tuple{Point{T},Point{T}}}, index::Int) where T
    select = SelectPoint(collection, index, construction[collection][index])
    push!(construction.objects, select)

    return GeoPointer{Point{T}}(length(construction.objects))
end

function update(construction::Construction, select::SelectPoint)
    return SelectPoint(select.collection, select.index, construction[select.collection][select.index])
end

struct PointOnCircle{T} <: GeoObject
    circle::GeoPointer{CircleMatrix{T}}

    param::Point{T} # far point polar to the diameter through the point
    point::Point{T}
    antipodal::Point{T}
end

state(poc::PointOnCircle) = poc.point

function new_point_on_circle!(construction::Construction{T}, circle::GeoPointer{CircleMatrix{T}}, nearest::Point) where T
    nearest = convert(Point{T}, nearest)
    
    circle_mat = construction[circle]

    pos = normalize_z(nearest)
    mid = normalize_z(center_of_circle(circle_mat))

    dir = pos - mid

    param = Point(dir.y, -dir.x, 0)

    diameter = cross(pos, mid)

    candidates = intersect_line_circle(diameter, circle_mat)

    d1 = projective_dist_min_scal(pos, candidates[1])
    d2 = projective_dist_min_scal(pos, candidates[2])

    if d2 < d1
        point = candidates[2]
        antipodal = candidates[1]
    else
        point = candidates[1]
        antipodal = candidates[2]
    end
    
    poc = PointOnCircle(circle, param, point, antipodal)
    push!(construction.objects, poc)

    return GeoPointer{Point{T}}(length(construction.objects))
end

function update_param(construction::Construction, poc::PointOnCircle, param::Point)
    circle_mat = construction[poc.circle]

    diameter = circle_mat * param

    candidates = intersect_line_circle(diameter, circle_mat)

    result = tracing2((poc.point, poc.antipodal), candidates)
    
    if result === nothing
        return nothing
    end

    return PointOnCircle(poc.circle, param, result...)
end

update(construction::Construction, poc::PointOnCircle) = update_param(construction, poc, poc.param)

function poc_param_from_nearest(poc::PointOnCircle, circle::CircleMatrix, nearest::Point)
    mid = normalize_z(center_of_circle(circle))
    dir = nearest - mid

    oldparam = poc.param
    oldpos = normalize_z(poc.point)
    olddir = oldpos - mid

    oldsign = oldparam.x * olddir.y - oldparam.y * olddir.x

    if real(oldsign) < 0
        dir = -dir
    end

    return Point(dir.y, -dir.x, 0)
end

function interpolate(construction::Construction, tr::Real, poc::PointOnCircle, target::Point)
    src = poc.param
    dst = poc_param_from_nearest(poc, construction[poc.circle], target)

    sp = dot(src, dst)

    if real(sp) >= 0
        t2 = tr * tr
        dt = 0.5 / (1 + t2)
        tc = complex((2 * tr) * dt + 0.5, (1 - t2) * dt)

        return src + tc * (dst - src)
    end

    mid = Point(src.y - dst.y, dst.x - src.x, 0)

    sp = dot(src, mid)

    if real(sp) < 0
        mid = -mid
    end

    if tr < 0
        tr = 2 * tr + 1
        t2 = tr * tr
        dt = 0.25 / (1 + t2)
        tc = complex((2 * tr) * dt + 0.25, (1 - t2) * dt)
    else
        tr = 2 * tr - 1
        t2 = tr * tr
        dt = 0.25 / (1 + t2)
        tc = complex((2 * tr) * dt + 0.75, (1 - t2) * dt)
    end

    uc = 1 - tc

    tc2 = tc * tc

    uc2 = uc * uc

    tuc = tc * uc

    return uc2 * src + tuc * mid + tc2 * dst
end

function move_to!(construction::Construction{T}, pointer::GeoPointer{Point{T}}, target::Point) where T
    target = convert(Point{T}, target)

    # we are trying to reach a new good state
    # thus we start by working on a copy of our Construction
    # and apply the changes only after we have successfully reached a new good state
    sandbox = Construction{T}(copy(construction.objects))

    # also remember last good state
    lastgood = Construction{T}(copy(construction.objects))

    last = -1.0
    step = 0.9

    t = last + step

    mover_index = pointer.index
    mover = construction.objects[mover_index]

    while last != t
        refine = false

        result = update_param(sandbox, sandbox.objects[mover_index], interpolate(construction, t, mover, target))

        if result !== nothing
            sandbox.objects[mover_index] = result
        
            # update all dependent objects
            for i = mover_index+1:length(sandbox.objects)
                if sandbox.objects[i] isa PointOnLine # special snowflake...
                    result = update(sandbox, sandbox.objects[i], construction, construction.objects[i])
                else
                    result = update(sandbox, sandbox.objects[i])
                end

                if result === nothing
                    refine = true
                    break
                end

                sandbox.objects[i] = result
            end
        else
            refine = true
        end

        if refine
            # this is bad. we could not update position for this object
            # throw current sandbox away, create new sandbox from lastgood and refine
            sandbox = Construction{T}(copy(lastgood.objects))

            step *= 0.5 # reduce step size
            t = last + step

            continue
        end

        last = t # successfully traced up to t
        step *= 1.25
        t += step

        if t >= 1
            t = 1.0
        end

        lastgood = sandbox
        sandbox = Construction{T}(copy(lastgood.objects))
    end

    if t != 1
        error("could not trace to end")
    end

    copyto!(construction.objects, lastgood.objects)
    return nothing
end
