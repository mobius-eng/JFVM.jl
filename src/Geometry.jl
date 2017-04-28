__precompile__()
module Geometry

import Base: +, -, *, /

export Point. tetrahedronVolume, triangleArea, tetrahedronCentroid, triangleCentroid, pyramidVolumeCentroid

# %% Point: point in 3D space

"""
Point in 3D space

Constructor

```julia
Point(x, y, [z=0])
Point(T)
```

where `x`, `y` and `z` are of some type `T<:Real`. `Point(T)` (e.g. `Point(Float64)`)
constructs an origin point of type `T`.
"""
immutable Point{T<:Real}
    x :: T
    y :: T
    z :: T
end

Point{T<:Real}(x::T, y::T) = Point(x,y, zero(T))
Point{T<:Real}(::Type{T}) = Point(zero(T), zero(T), zero(T))
"""
```
Point(x,y,z)
```
Construct an array of points given the array of `x`, `y` (and `z`) coordinates.
If `z` coordinates are not provided, they are assumed to be zero.
"""
Point{T<:Real}(x :: AbstractArray{T}, y :: AbstractArray{T}, z :: AbstractArray{T}) = map(Point, x, y, z)

Point{T<:Real}(x :: AbstractArray{T}, y :: AbstractArray{T}) = map(Point, x, y)

function Base.broadcast{T<:Real}(f, p :: Point{T})
    Point(f(p.x), f(p.y), f(p.z))
end

function Base.broadcast{T1<:Real, T2<:Real}(f, p1 :: Point{T1}, p2 :: Point{T2})
    Point(f(p1.x, p2.x), f(p1.y, p2.y), f(p1.z, p2.z))
end

function Base.broadcast{T1<:Real, T2}(f, x :: T2, p :: Point{T1})
    Point(f(x,p.x), f(x,p.y), f(x,p.z))
end

function Base.broadcast{T1<:Real,T2}(f, p :: Point{T1}, x :: T2)
    Point(f(p.x,x), f(p.y,x), f(p.z,x))
end

Base.norm(p :: Point) = norm([p.x, p.y, p.z])

function Base.cross{T<:Real}(p :: Point{T}, q :: Point{T})
    c = [p.x, p.y, p.z] × [q.x, q.y, q.z]
    Point(c...)
end

function Base.dot{T<:Real}(p :: Point{T}, q :: Point{T})
    p.x*q.x + p.y*q.y + p.z*q.z
end

function +{T<:Real}(p :: Point{T}, q :: Point{T})
    Point(p.x+q.x, p.y+q.y, p.z+q.z)
end

function -{T<:Real}(p :: Point{T}, q :: Point{T})
    Point(p.x-q.x, p.y-q.y, p.z-q.z)
end

function /{T<:Real}(p :: Point{T}, f :: Real)
    Point(p.x/f, p.y/f, p.z/f)
end

function *{T<:Real}(f :: Real, p :: Point{T})
    Point(f*p.x, f*p.y, f*p.z)
end

function Base.isapprox{T<:Real}(p :: Point{T}, q :: Point{T}; kwargs...)
    isapprox(p.x, q.x; kwargs...) && isapprox(p.y, q.y; kwargs...) && isapprox(p.z, q.z; kwargs...)
end

function tetrahedronVolume{T<:Real}(a :: Point{T}, b :: Point{T}, c :: Point{T}, d :: Point{T})
    p = b - a
    q = c - a
    r = d - a
    abs(p.x*(q.y*r.z-q.z*r.y) - p.y*(q.x*r.z-q.z*r.x) + p.z*(q.x*r.y-q.y*r.x)) / 6
end

function pyramidVolumeCentroid{T<:Real}(a :: Point{T}, b, c, d, p)
    # Tetrahedron abcp
    v1 = tetrahedronVolume(a,b,c,p)
    c1 = tetrahedronCentroid(a,b,c,p)
    # Tetrahedron acdp
    v2 = tetrahedronVolume(a,c,d,p)
    c2 = tetrahedronCentroid(a,c,d,p)
    # Comdine
    v = v1 + v2
    (v, (v1/v)*c1 + (v2/v)*c2)

end

function tetrahedronCentroid{T<:Real}(a :: Point{T}, b :: Point{T}, c :: Point{T}, d :: Point{T})
    (a+b+c+d) / (4one(T))
end

function triangleArea(a :: Point{T}, b :: Point{T}, c ::Point{T})
    Sf = (b - a) × (c - a)
    area = norm(Sf)
    (area, Sf)
end

function triangleCentroid(a :: Point{T}, b :: Point{T}, c ::Point{T})
    (a+b+c) / (3one(T))
end

end
