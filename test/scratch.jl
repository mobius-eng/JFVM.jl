module JFVMMScratch

import JFVMM

FV = JFVMM

m = FV.createMesh2D(3.0, 2.0, 3, 4)


m3 = FV.createMesh3D(3.0, 2.0, 4.0, 3, 4, 2)

m3.vertices

end


import JFVMM

Fv = JFVMM

m = Fv.createMesh2D(3.0, 2.0, 3, 4)


type Foo
    x :: Int64
    y :: Int64
end

f = Foo(1,2)

f.x = 5

f
