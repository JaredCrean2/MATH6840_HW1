# solve 1.2.1 in the textbook

function driver()
  xmin = 0
  xmax = 1
  tmax = 0.10
  N = 161
  r = 0.4
#  delta_t = 0.02
  nu = 1
  ICFunc = ICSin
  BCL = BCZero
  BCR = BCZero

  # for accurate performance benchmarking, JIT the function first
  u, tmax_ret = solve(xmin, xmax, tmax, N, r, nu, ICFunc, BCL, BCR, write_files=false)

  # now that the function is compiled, we wil get an accurate timing
  u, tmax_ret = solve(xmin, xmax, tmax, N, r, nu, ICFunc, BCL, BCR, write_files=true)

  delta_x = (xmax - xmin)/(N-1)
  delta_t = (r*delta_x^2)/nu  # nu*delta_t

  println("tmax = ", tmax, ", tmax_ret = ", tmax_ret)
  vals = [xmin, xmax, tmax_ret, nu, delta_t]
  writedlm("counts.dat", vals)
  writedlm("u.dat", u[2:(end-1), end])
end

function solve(xmin, xmax, tmax, N, r, nu, ICFunc::Function, BCL::Function, BCR::Function; write_files=true)
# xmin = minimum x coordinate
# xmax = maximum x coordinate
# tmax = maximum time value
# N = number of x points
# r = nu*delta_t/delta_x^2
# ICFunc = initial condition function with signature val = ICFunc(x)
# BCL = left boundary condition function with signature val = BCL(t)
# BCR = right boundary condition function
# write_files = whether or not to write to convergence.dat

delta_x = (xmax - xmin)/(N-1)
delta_t = (r*delta_x^2)/nu  # nu*delta_t
#r = nu*delta_t/(delta_x^2)
nStep = convert(Int, div(tmax, delta_t)) + 1

println("delta_x = ", delta_x)
println("delta_t = ", delta_t)
println("r = ", r)
println("nStep = ", nStep)


# allocate storage
u = Array(Float64, N+2, nStep)

# apply IC
# Not applying BC at initial condition
# hopefully IC and BC are consistent
for i=2:(N+1)
  u[i, 1] = ICFunc(xmin + (i-2)*delta_x)
end


flops = 0 # count number of FLOPs

time = @elapsed for tstep=2:nStep  # loop over timesteps

#  println("tstep = ", tstep)
  tstep_1 = tstep-1

  # apply BC
  t_i = (tstep - 1)*delta_t
#  u[1, tstep] = BCL(t_i)
  applyCompatBC(BCL, 2, t_i, tstep, delta_x, nu, u)
  applyNeumannBC(BCR, N+1, t_i, tstep, delta_x, u; is_left=false)

#  println("u at tstep $tstep_1 = ", u[:, tstep_1])
  # calculate non-ghost points
  for i=2:(N+1)
    u_k = u[i, tstep - 1]
    u_k_1 = u[i-1, tstep - 1 ]
    u_k_p1 = u[i+1, tstep - 1]
    u[i, tstep] = u[i, tstep - 1] + r*(u_k_p1 - 2*u_k + u_k_1) 
  end

  flops += 5*(N-2)

end

println("time = ", time)
flop_rate = 1e-6*flops/time  # MFlops/seconds
println("flop rate = ", flop_rate, " MFlops/sec")


# calculate maximum error
err_max = typemin(Float64)
for i=2:(N+1)
  x_i = xmin + (i-2)*delta_x
  ue = uExact(x_i, delta_t*(nStep - 1), nu)
  err_i = abs(u[i, nStep] - ue)
  if err_i > err_max
    err_max = err_i
  end
end

println("final error = ", err_max)

# write to file
if write_files
  f = open("convergence.dat", "a")
  @printf(f, "%d %16.15f %16.15f %d\n", N, err_max, time, nStep)
  close(f)
end
return u, delta_t*(nStep - 1)

end


function applyNeumannBC(f::Function, i::Integer, t, tstep, delta_x, u; is_left=true)
# apply Neumann BC to point i in the array u using ghost point
# there must be an element u[i+1] to use as a ghost point
# t is the current time
#  println("applying Neumann BC")
  f_val = f(t)

  if is_left
    ghost_val = u[i+1, tstep - 1] - 2*delta_x*f_val
    u[i-1, tstep - 1] = ghost_val
  else
    ghost_val = u[i-1, tstep - 1] + 2*delta_x*f_val
    idx = i+1
    idt = tstep-1
    u[i+1, tstep - 1] = ghost_val
  end

  return nothing
end


function applyCompatBC(f::Function, i::Integer, t, tstep, delta_x, nu, u; is_left=true)
# apply a Direchlet compatable BC to the array u using a ghost point
# f is the Direchlet function

#  println("applying compatability BC")
  # use complex step to get the df/dt
  epsilon = 1e-20
  val_pert = t + complex(0, epsilon)
  dfdt = imag( f(val_pert) )/epsilon

  if is_left
    ghost_val = delta_x*delta_x*dfdt/nu - u[i+1, tstep - 1] + 2*u[i, tstep - 1]
    u[i-1, tstep - 1] = ghost_val
  else
    ghost_val = delta_x*delta_x*dfdt/nu - u[i-1, tstep - 1] + 2*u[i, tstep - 1]
    u[i+1, tstep - 1] = ghost_val
  end

  return nothing
end



function ICSin(x)
  return sin(5*pi*x/2)
end

function BCZero(x)
  return 0.0
end

function uExact(x, t, nu)
# calculate exact solution for an IC of sin(5*pi/2x)
  k = 5*pi/2
  return exp(-nu*k*k*t)*sin(k*x)
end

function contains(A, val)
  for i=1:length(A)
    if A[i] == val
      return true
    end
  end

  return false
end

# run
driver()
