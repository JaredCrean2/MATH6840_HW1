# solve Problem Set 7 Problem 1
using ArrayViews  # non-copying subarray package

function driver()
  xmin = 0
  xmax = 1
  N = 10  # N = # of grid points (not including ghosts)
  delta_x = (xmax - xmin)/N
#  r = 0.5
#  sigma = 0.75
  delta_t = delta_x
#  nu = 1/6


  tmax = delta_t + eps()
  ICFunc = IC1
  BCL = BC1
  BCR = BC2
  src = SRC

  u, tmax_ret = solve(xmin, xmax, tmax, N, delta_t, ICFunc, BCL, BCR, src)
#  u = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
  println("tmax_ret = ", tmax_ret)
  max_err, l2_err = calcError(u, xmin, xmax, tmax_ret, N)
  vals = [xmin, xmax, tmax_ret, delta_t]
  writedlm("counts.dat", vals)
  writedlm("u.dat", u)

  f = open("convergence.dat", "a+")
  @printf(f, "%d %16.15f %16.15f\n", N, max_err, tmax_ret)
  close(f)
end

function calcError(u, xmin, xmax, tmax, N)
# calculate the error
# u should not include the ghost point

  print("\n")
  delta_x = (xmax - xmin)/N
  err = zeros(size(u))
  u_ex = zeros(size(u))
  for i=1:length(u)
    x_i = xmin + (i-1/2)*delta_x
    u_exact_i = uExact(x_i, tmax)
    u_ex[i] = u_exact_i
    err[i] = u_exact_i - u[i]
  end

  err_max = maximum(abs(err))
  l2_err = norm(err)

  println("err = \n", err)
  println("u = \n", u)
  println("u_ex = \n", u_ex)



  println("max error = ", err_max)
  println("L2 error = ", l2_err)

  writedlm("u.dat", u)
  writedlm("uexact.dat", u_ex)
  return err_max, l2_err
end



function solve(xmin, xmax, tmax, N, delta_t,ICFunc::Function, BCL::Function, BCR::Function, source::Function)
# xmin = minimum x coordinate
# xmax = maximum x coordinate
# tmax = maximum time value
# N : N+1 =  number of x points
# r = nu*delta_t/delta_x^2
# sigma = delta_t/delta_x
# ICFunc = initial condition function with signature val = ICFunc(x)
# BCL = left boundary condition function with signature val = BCL(t)
# BCR = right boundary condition function
# source = source function, signature val = source(x, t)
# this function assumes no ghost point on the left, ghost point on the right
delta_x = (xmax - xmin)/N
#delta_t = (r*delta_x^2)/nu  # nu*delta_t
r = delta_t/(2*delta_x^2)
nStep = convert(Int, div(tmax, delta_t))
mat_size = N+2

println("tmax = ", tmax)
println("delta_x = ", delta_x)
println("delta_t = ", delta_t)
println("r = ", r)
println("nStep = ", nStep)


# allocate storage
A = zeros(Float64, mat_size, mat_size)  # this could be a SparseMatrixCSC
rhs = Array(Float64, mat_size) # right hand side
u_i = Array(Float64, mat_size)  # current timestep solution values

# apply IC
# Not applying BCL at initial condition
for i=2:(mat_size-1)
  u_i[i] = ICFunc(xmin + (i-3/2)*delta_x)
end

# set ghost point value at IC
u_i[1] = -delta_x*BCL(0) + u_i[2]
u_i[mat_size] = delta_x*BCR(0) + u_i[mat_size-1]


println("\nu_initial = \n", u_i)

# apply BC to IC

# construct the matrix A
# do left BC
A[1, 1] = -1/delta_x
A[1, 2] = 1/delta_x

# do right BC
A[mat_size, mat_size] = 1/(delta_x)
A[mat_size, mat_size-1] = -1/(delta_x) 


for i=2:(mat_size-1)  # loop over interior of matrix
  # evaluate coefficient D
  x_j = (i-3/2)*delta_x  # x coordinate at j
  x_jp = (i-1/2)*delta_x  # x coordinate at j + 1
  x_jm = (i-5/2)*delta_x  # x coordinate at j - 1
  D_j = calcD(x_j)
  D_jp = calcD(x_jp)
  D_jm = calcD(x_jm)

  D_jph = 0.5*(D_j + D_jp)  # Dj+1/2
  D_jmh = 0.5*(D_j + D_jm)  # Dj-1/2

  stencil_l = -r*D_jmh
  stencil_c = 1 + r*(D_jph + D_jmh)
  stencil_r = -r*D_jph

  A[i, i-1] = stencil_l
  A[i, i] = stencil_c
  A[i, i+1] = stencil_r
end


println("A = \n", A)

Af = lufact(A)
#println("Af = ", Af)
println("typeof(Af) = ", typeof(Af))


# set up stencil for rhs
print("\n")

time = @elapsed for tstep=1:nStep  # loop over timesteps
# advance from timestep tstep to tstep + 1

  println("\ntstep = ", tstep)
#  print("\n")


  # calculate right hand size interior points
  for i=2:(mat_size-1)
    # get solution values
    u_k = u_i[i]
    u_k_1 = u_i[i-1]
    u_k_p1 = u_i[i+1]
    src_val = source( (i-3/2)*delta_x, (tstep - 0.5)*delta_t) 

    # get coefficient values
    x_j = (i-3/2)*delta_x
    x_jp = (i-1/2)*delta_x  # x coordinate at j + 1
    x_jm = (i-5/2)*delta_x  # x coordinate at j - 1
    D_j = calcD(x_j)
    D_jp = calcD(x_jp)
    D_jm = calcD(x_jm)

    D_jph = 0.5*(D_j + D_jp)  # Dj+1/2
    D_jmh = 0.5*(D_j + D_jm)  # Dj-1/2

    stencil_l = r*D_jmh
    stencil_c = 1 - r*(D_jph + D_jmh)
    stencil_r = r*D_jph


    rhs[i] = stencil_l*u_k_1 + stencil_c*u_k  + stencil_r*u_k_p1 + delta_t*src_val
    
  end

  # apply BC terms to rhs
  rhs[1] = BCL(tstep*delta_t)
  rhs[mat_size] = BCR(tstep*delta_t)

  println("rhs = \n", rhs)
  # solve for next time step
#  A_ldiv_B!(Af, rhs)  # rhs gets overwritten with new solution values
  u_i = A\rhs
#  println("secondard u_next = \n", u_next)

  println("u_next = \n", u_i)



end

# apply BCs to final time
#u_i_1[1] = BCL((nStep - 1)*delta_t)
#u_i_1[N] = BCR((nStep - 1)*delta_t)


println("time = ", time)

return u_i[2:end-1], delta_t*(nStep)  # plus 1 because we are at the beginning of the next timestep

end


#=
function IC1(x)
  return cos(x)
end

#function BCZero(x)
#  return 0.0
#end

function BC1(t)
  return 0
end

function BC2(t)
  return -sin(1)*cos(t)
end

function SRC(x, t)
#  return -cos(x)*sin(t) + cos(x)*cos(t) + sin(x)*sin(x)*cos(t) - cos(x)*cos(x)*cos(t)
  return cos(x)*cos(t) + 2*cos(x)*sin(t)
end

function uExact(x, t)
  return cos(x)*cos(t)
end
=#

function calcD(x)
  return 1e6
end


function IC1(x)
  return x*x + x + 1
end


function BC1(t)
  return 1
end

function BC2(t)
  return 3
end

function SRC(x, t)
#  return -0.1*sin(x)
   return 2*t - 1e6*2 + 1
end

function uExact(x, t)
  return x*x + x + t*t + t + 1
end



# run
driver()
