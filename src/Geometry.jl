__precompile__()
module Geometry

using StaticArrays

export Point2D, triangleArea, triangleCentroid

# %% Point: point in 2D space

"""
Point in 2D space

Constructor

```julia
Point2D(x, y)
```

where `x` and `y` are `Float64` values, constructs a 2D point (or otherwise viewed as
the vector from origin to this point).
"""
immutable Point2D <: FieldVector{Float64}
    x :: Float64
    y :: Float64
end

"""
```julia
Point2D(x,y)
```
Construct an array of points given the array of `x`, `y` coordinates.
"""
Point2D{U <: Real, V <: Real}(x :: AbstractArray{U}, y :: AbstractArray{V}) = map(Point2D, x, y)

StaticArrays.similar_type(::Type{Point2D}, t :: Type{Float64}, sz :: Size{(2,)}) = Point2D

"""
```julia
cross( p :: Point2D, q :: Point2D)
```
Cross product of two 2D vectors is the (signed) area of the parallelegram constructed
on these vectors
"""
Base.cross(p :: Point2D, q :: Point2D) = p.x*q.y - p.y*q.x

"""
Computes the area of the triangle defined by its vertices
"""
function triangleArea(a :: Point2D, b :: Point2D, c ::Point2D)
    Sf = ((b - a) × (c - a)) / 2
    abs(Sf)
end

"""
Computes the centroid of the triangle defined by its vertices
"""
function triangleCentroid(a :: Point2D, b :: Point2D, c ::Point2D)
    (a+b+c) / 3
end

end
