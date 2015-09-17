# solve 1.2.1 in the textbook

function driver()
  xmin = 0
  xmax = 1
  tmax = 0.901
  N = 11
#  r = 0.5
  delta_t = 0.02
  nu = 1/6
  ICFunc = ICSin
  BCL = BCZero
  BCR = BCZero

  u, tmax_ret = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
#  u = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
  vals = [xmin, xmax, tmax_ret, nu, delta_t]
  writedlm("counts.dat", vals)
  writedlm("u.dat", u)
end

function solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc::Function, BCL::Function, BCR::Function)
# xmin = minimum x coordinate
# xmax = maximum x coordinate
# tmax = maximum time value
# N = number of x points
# r = nu*delta_t/delta_x^2
# ICFunc = initial condition function with signature val = ICFunc(x)
# BCL = left boundary condition function with signature val = BCL(t)
# BCR = right boundary condition function

delta_x = (xmax - xmin)/(N-1)
#delta_t = (r*delta_x^2)/nu  # nu*delta_t
r = nu*delta_t/(delta_x^2)
nStep = convert(Int, div(tmax, delta_t)) + 1

println("delta_x = ", delta_x)
println("delta_t = ", delta_t)
println("r = ", r)
println("nStep = ", nStep)


# allocate storage
u = Array(Float64, N, nStep)

# apply IC
# Not applying BC at initial condition
# hopefully IC and BC are consistent
for i=1:N
  u[i, 1] = ICFunc(xmin + (i-1)*delta_x)
end


# use low order method for first step
# calculate interior points
tstep = 2
for i=2:(N-1)
  u_k = u[i, tstep - 1]
  u_k_1 = u[i-1, tstep - 1 ]
  u_k_p1 = u[i+1, tstep - 1]
  u[i, tstep] = u[i, tstep - 1] + r*(u_k_p1 - 2*u_k + u_k_1)
  
end




flops = 0

time = @elapsed for tstep=3:nStep  # loop over timesteps

#  println("tstep = ", tstep)

  # apply BC
  u[1, tstep] = BCL((tstep - 1)*delta_t)
  u[N, tstep] = BCR((tstep - 1)*delta_t)
#   u[1, tstep] = 0.0
#   u[N, tstep] = 0.0


  # calculate interior points
  for i=2:(N-1)
    u_k = u[i, tstep - 1]
    u_k_1 = u[i-1, tstep - 1 ]
    u_k_p1 = u[i+1, tstep - 1]
    u[i, tstep] = u[i, tstep - 2] + r*(u_k_p1 - 2*u_k + u_k_1)
    
  end

  flops += 5*(N-2)

end

println("time = ", time)
flop_rate = 1e-6*flops/time  # MFlops/seconds
println("flop rate = ", flop_rate, " MFlops/sec")

return u, delta_t*(nStep - 1)

end



function ICSin(x)
  return sin(2*pi*x)
end

function BCZero(x)
  return 0.0
end



# run
driver()
