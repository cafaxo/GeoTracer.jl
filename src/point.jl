import Base: abs, *, /, +, -, convert
import LinearAlgebra: dot, cross, transpose

struct Point{T}
    x::T
    y::T
    z::T

    function Point(x::Number, y::Number, z::Number)
        x_p, y_p, z_p = promote(x, y, z)
        return new{typeof(x_p)}(x_p, y_p, z_p)
    end
end

convert(::Type{Point{T}}, p::Point) where T = Point(convert(T, p.x), convert(T, p.y), convert(T, p.z))

Point(x::SVector) = Point(x[1], x[2], x[3])

svector(p::Point) = SVector((p.x, p.y, p.z))

cross(a, b) = Point(cross(svector(a), svector(b)))

+(a::Point, b::Point) = Point(a.x + b.x, a.y + b.y, a.z + b.z)
-(a::Point, b::Point) = Point(a.x - b.x, a.y - b.y, a.z - b.z)

-(p::Point) = Point(-p.x, -p.y, -p.z)

abs(p::Point) = norm(svector(p))

*(s::Number, p::Point) = Point(s * p.x, s * p.y, s * p.z)

/(p::Point, s::Number) = Point(p.x / s, p.y / s, p.z / s)

@inline *(A::SMatrix{3,3}, p::Point) = Point(A * svector(p))
@inline *(A::SDiagonal{3}, p::Point) = Point(A * svector(p))

conjugate(p::Point) = Point(p.x', p.y', p.z')

dot(a::Point, b::Point) = a.x * b.x + a.y * b.y + a.z * b.z

normalize_abs(p::Point) = inv(abs(p)) * p

normalize_z(p::Point) = inv(p.z) * p

distance(a::Point, b::Point) = abs(normalize_z(a) - normalize_z(b))

function cross_operator(p::Point)
    return @SMatrix [ 0   -p.z  p.y;
                      p.z  0   -p.x;
                     -p.y  p.x  0]
end

const CircleMatrix{T} = SMatrix{3, 3, T, 9}

center_of_circle(c::CircleMatrix) = Point(c[3,1], c[3,2], -c[1,1])

function circle_from_mp(m::Point, p::Point)
    x = m.x
    y = m.y

    mz = -m.z

    tang = @SMatrix [
        mz 0  x;
        0  mz y;
        x  y  0
    ]

    fundmat = SDiagonal(SVector(0, 0, 1))

    mu = dot(tang * p, p)
    
    la = dot(fundmat * p, p)

    m1 = mu * fundmat
    m2 = la * tang

    return m1 - m2
end

function circle_from_mr(m::Point, r::Real)
    sr = m.z * r
    rr = sr * sr

    return @SMatrix [
        m.z*m.z   0       -m.z*m.x;
        0         m.z*m.z -m.z*m.y;
        -m.z*m.x -m.z*m.y  m.x*m.x + m.y*m.y - rr
    ]
end

function intersect_line_circle(l::Point, c::CircleMatrix)
    l1 = cross_operator(l)

    s = transpose(l1) * (c * l1)

    maxidx = argmax((abs(l.x), abs(l.y), abs(l.z)))

    if maxidx == 1
        a11 = s[2,2]
        a12 = s[2,3]
        a21 = s[3,2]
        a22 = s[3,3]
        b = l.x
    elseif maxidx == 2
        a11 = s[1,1]
        a12 = s[1,3]
        a21 = s[3,1]
        a22 = s[3,3]
        b = l.y
    else
        a11 = s[1,1]
        a12 = s[1,2]
        a21 = s[2,1]
        a22 = s[2,2]
        b = l.z
    end

    alp = sqrt(a12 * a21 - a11 * a22) / b

    erg = s + alp * l1

    erg1 = Point(erg[argmax((norm(erg[1,:]), norm(erg[2,:]), norm(erg[3,:]))),:])
    erg2 = Point(erg[:,argmax((norm(erg[:,1]), norm(erg[:,2]), norm(erg[:,3])))])

    return erg1 / abs(erg1), erg2 / abs(erg2)
end

function intersect_circle_circle(c1::CircleMatrix, c2::CircleMatrix)
    ct1 = c2[1,1]
    line1 = ct1 * c1[3,:]

    ct2 = c1[1,1]
    line2 = ct2 * c2[3,:]

    ll = line1 - line2

    return intersect_line_circle(Point(ll[1], ll[2], ll[3] / 2), c1)
end

function projective_dist_min_scal(a::Point, b::Point)
    sa = abs(a)
    sb = abs(b)

    if sa == 0 || sb == 0
        return zero(sa)
    end

    cb = conjugate(b)

    p = dot(a, cb)

    if abs(p) < 1e-13
        #error("p too small")
        p = 1.0
    end

    np = p / abs(p)

    na = a / sa
    nb = np * (b / sb)

    d1 = abs(na + nb)
    d2 = abs(na - nb)

    return min(d1, d2)
end

function tracing2(old::Tuple{Point{T},Point{T}}, new::Tuple{Point{T},Point{T}}) where T
    o1, o2 = old
    n1, n2 = new

    safety = 3

    do1n1 = projective_dist_min_scal(o1, n1)
    do1n2 = projective_dist_min_scal(o1, n2)
    do2n1 = projective_dist_min_scal(o2, n1)
    do2n2 = projective_dist_min_scal(o2, n2)
    do1o2 = projective_dist_min_scal(o1, o2)
    dn1n2 = projective_dist_min_scal(n1, n2)
    cost1 = do1n1 + do2n2
    cost2 = do1n2 + do2n1

    # Always sort output: we don't know yet whether it's correct, but
    # it's our best bet.

    if cost1 > cost2
        res = (n2, n1)
        cost = cost2
    else
        res = (n1, n2)
        cost = cost1
    end

    # is nan test?

    if do1o2 > cost * safety && dn1n2 > cost * safety
        # Distance within matching considerably smaller than distance
        # across matching, so we could probably match correctly.
        return res
    elseif dn1n2 < 1e-5
        error("New points too close: we presumably are inside a singularity.")
    elseif do1o2 < 1e-5
        error("Moved out of singularity.")
    end

    # refine
    return nothing
end

function project_point_to_line(point::Point, line::Point)
    return cross(cross(Point(line.x, line.y, 0), point), line)
end

function point_reflection(center::Point, point::Point)
    return 2 * point.z * center - center.z * point
end
