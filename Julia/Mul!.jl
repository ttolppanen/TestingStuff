using LinearAlgebra

function f()
   d = 1000
   v = rand(d)
   m = rand(d, d) 
   @time mul!(v, m, v)
   @time v .= m * v
   @time v = m * v
end

f();