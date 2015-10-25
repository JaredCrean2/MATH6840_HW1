# solve 1.2.1 in the textbook

function driver()
  xmin = 0
  xmax = 1
  tmax = .901
  N = 11
#  r = 0.5
  delta_t = 0.02
  nu = 1/6
  ICFunc = ICSin
  BCL = BCD2
  BCR = BCD2

  u, tmax_ret = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
#  u = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
  vals = [xmin, xmax, tmax_ret, nu, delta_t]
  writedlm("counts.dat", vals)
  writedlm("u.dat", u[2:end-1])
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
n = N+2  # number of points in u
u_i = Array(Float64, n) # current timestep
u_i_1 = Array(Float64, n)  # previous timestpe

# apply IC
# Not applying BC at initial condition
# hopefully IC and BC are consistent
for i=2:(n-1)
  println("x = ", (i-2)*delta_x)
  u_i_1[i] = ICFunc(xmin + (i-2)*delta_x)
end

flops = 0

time = @elapsed for tstep=2:nStep  # loop over timesteps

#  println("tstep = ", tstep)

  # apply BC to u_i_1
  BCL(u_i_1)
  BCR(u_i_1, isleft=false)

  println("u_i_1 =\n", u_i_1)

  # calculate interior points
  for i=3:(n-2)
    u_k = u_i_1[i]
    u_k_1 = u_i_1[i-1]
    u_k_2 = u_i_1[i-2]
    u_k_p1 = u_i_1[i+1]
    u_k_p2 = u_i_1[i+2]

#    u_i[i] = u_k + r*(u_k_p1 - 2*u_k + u_k_1)
    u_i[i] = u_k + r*(-u_k_2/6 + 5*u_k_1/3 - 3*u_k + 5*u_k_p1/3 - u_k_p2/6)
  end

  flops += 12*(n-4)

  # rebind names, using u_i as u_i_1, and reusing the array that was
  # u_i_1 as u_i for the next timestep
  tmp = u_i_1
  u_i_1 = u_i
  u_i = tmp



end

#apply BC to final time
BCL(u_i_1)
BCR(u_i_1, isleft=false)



# apply BC to fina
println("time = ", time)
flop_rate = 1e-6*flops/time  # MFlops/seconds
println("flop rate = ", flop_rate, " MFlops/sec")

return u_i_1, delta_t*(nStep - 1)  # return u_i_1

end


function BCZero(u; isleft=true)
# apply direchlet BC to boundary point, numerical BC =0 to ghost point

  if isleft
    u[2] = 0  # direchlet BC
    u[1] = 0  # numerical BC
  else
    n = length(u)
    u[n-1] = 0 # dirchlet BC
    u[n] = 0  # numerical BC
  end

  return nothing
end

function BCD2(u; isleft=true)
# apply direchlet BC to boundary pointer, numerical u''=0 to ghost point

  if isleft
    u[2] = 0  # dirchlet BC
    u[1] = 2*u[2] - u[3]  # numerical BC
  else
    n = length(u)
    u[n-1] = 0  # dirchlet BC
    u[n] = 2*u[n-1] - u[n-2]
  end

  return nothing
end


function ICSin(x)
  return sin(2*pi*x)
end



# run
driver()
