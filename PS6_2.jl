# solve Problem Set 6 Problem 2
using ArrayViews  # non-copying subarray package

function driver()
  xmin = 0
  xmax = pi
  ymin = 0
  ymax = pi
  N = 1000  # N+1 = # of grid points in each direction
  delta_x = (xmax - xmin)/N
#  r = 0.5
#  sigma = 0.75
  delta_t = delta_x
#  nu = 1/6


  tmax = 1.0
  ICFunc = IC1
  BCxl = BC1
  BCxr = BC2
  BCyl = BC3
  BCyr = BC4
  src = SRC

  u, tmax_ret = solve(xmin, xmax, ymin, ymax, tmax, N, delta_t, ICFunc, BCxl, BCxr, BCyl, BCyr, src)
#  u = solve(xmin, xmax, tmax, N, delta_t, nu, ICFunc, BCL, BCR)
  println("tmax_ret = ", tmax_ret)
  u_vals = u[ 2:(end-1), :]
  max_err, l2_err = calcError(u_vals, xmin, xmax, ymin, ymax, tmax_ret, N)
  vals = [xmin, xmax, tmax_ret, delta_t]
  writedlm("counts.dat", vals)
#  writedlm("u.dat", u)

  f = open("convergence.dat", "a+")
  @printf(f, "%d %16.15f %16.15f\n", N, max_err, tmax_ret)
  close(f)
end

function calcError(u, xmin, xmax, ymin, ymax, tmax, N)
# calculate the error
# u should not include the ghost point

  println("size(u) = ", size(u))
  print("\n")
  delta_x = (xmax - xmin)/N
  delta_y = (ymax - ymin)/N
  err = zeros(size(u))
  u_ex = zeros(size(u))

  for i=1:size(u, 1)  # loop over y
    y_i = ymin + (i-1)*delta_y
    for j=1:size(u, 2)  # loop over x
      x_j = xmin + (j-1)*delta_x
      u_exact_i = uExact(x_j, y_i, tmax)
      u_ex[i, j] = u_exact_i
      err[i, j] = u_exact_i - u[i, j]
    end
  end

  err_max = maximum(abs(err))
  l2_err = norm(err)

#  println("err = \n", err)
#  println("u = \n", u)
#  println("u_ex = \n", u_ex)



  println("max error = ", err_max)
  println("L2 error = ", l2_err)

#  println("size(u) = ", size(u))
  writedlm("u.dat", u)
  writedlm("uexact.dat", u_ex)
  return err_max, l2_err
end



function solve(xmin, xmax, ymin, ymax, tmax, N, delta_t,ICFunc::Function, BCxl::Function, BCxr::Function, BCyl::Function, BCyr::Function, source::Function)
# xmin = minimum x coordinate
# xmax = maximum x coordinate
# tmax = maximum time value
# N : N+1 =  number of grid points in each direction
# r = nu*delta_t/delta_x^2
# sigma = delta_t/delta_x
# ICFunc = initial condition function with signature val = ICFunc(x, y)
# BCxl = x direction left boundary condition function with signature val = BCL(y, t)
# BCxr = x direction right boundary condition function
# BCyl = y direction left (bottom) BC function with signature val = BCyl(x, t)
# BCyr = y direction right (top) BC function
# source = source function, signature val = source(x, y, t)
# this function assumes no ghost point on the left, ghost point on the right
delta_x = (xmax - xmin)/N
delta_y = (ymax - ymin)/N
#delta_t = (r*delta_x^2)/nu  # nu*delta_t
rx = delta_t/(delta_x^2)
ry = delta_t/(delta_y^2)
nStep = convert(Int, div(tmax, delta_t))
mat_sizex = N+1  # x_i = (index-1)*delta_x
mat_sizey = N+3  # y_j = (index-2)*delta_y because of ghost points

println("tmax = ", tmax)
println("delta_x = ", delta_x)
println("delta_y = ", delta_y)
println("delta_t = ", delta_t)
println("rx = ", rx)
println("ry = ", ry)
println("nStep = ", nStep)


# allocate storage
Ax = zeros(Float64, mat_sizex, mat_sizex)  # this could be a SparseMatrixCSC
Ay = zeros(Float64, mat_sizey, mat_sizey)
rhs_x = Array(Float64, mat_sizex) # right hand side for x matrix
rhs_y = Array(Float64, mat_sizey)  # right hand side for y matrix
u_i = Array(Float64, mat_sizey, mat_sizex)  # current timestep solution values
u_i_1 = Array(Float64, mat_sizey, mat_sizex)  # next timestep solution values

# apply IC
# Not Dirichlet BCs at initial condition
# hopefully IC and BC are consistent
for i=1:mat_sizex  # loop over x points
  x_i = xmin + (i-1)*delta_x
  for j=2:(mat_sizey-1)
    y_j = ymin + (j-2)*delta_y
    u_i[j, i] = ICFunc(x_i, y_j)
  end
end

# set ghost point value at IC
for i=1:mat_sizex  # loop along x direction
  x_j = xmin + (i-1)*delta_x
  u_i[1, i] = -2*delta_x*BCyl(x_j, 0.0) + u_i[3, i]  # bottom BC
  u_i[mat_sizey, i] = 2*delta_x*BCyr(x_j, 0.0) + u_i[mat_sizey-2, i]
end


#println("\nu_initial = \n", u_i)

# apply BC to IC

# construct Ax, Ay
# Ax BCs
Ax[1, 1] = 1
Ax[mat_sizex, mat_sizex] = 1

# Ay BCs
# do right BC
Ay[1, 1] = -1/(2*delta_y)
Ay[1, 3] = 1/(2*delta_y)

Ay[mat_sizey, mat_sizey] = 1/(2*delta_y)
Ay[mat_sizey, mat_sizey-2] = -1/(2*delta_y)


# assemble Ax
stencil_l = -rx/2
stencil_c = 1 + rx
stencil_r = -rx/2

for i=2:(mat_sizex-1)  # loop over interior Ax
  Ax[i, i-1] = stencil_l
  Ax[i, i] = stencil_c
  Ax[i, i+1] = stencil_r
end


# assemble Ay
stencil_l = -ry/2
stencil_c = 1 + ry
stencil_r = -ry/2

for i=2:(mat_sizey-1)  # loop over interior of Ay
  Ay[i, i-1] = stencil_l
  Ay[i, i] = stencil_c
  Ay[i, i+1] = stencil_r
end


Axf = lufact(Ax)
Ayf = lufact(Ay)

#println("Ax = \n", Ax)
#println("Ay = \n", Ay)

#println("A = \n", A)

#Af = lufact(A)
#println("Af = ", Af)
#println("typeof(Af) = ", typeof(Af))


# set up stencils for rhs
# stencil for Ax rhs
stencilx_l = ry/2
stencilx_c = 1 - ry
stencilx_r = ry/2

# stencil for Ay rhs
stencily_l = rx/2
stencily_c = 1 - rx
stencily_r = rx/2


println("\nstencilx_l = ", stencilx_l)
println("stencilx_c = ", stencilx_c)
println("stencilx_r = ", stencilx_r)
println("\nstencily_l = ", stencily_l)
println("stencily_c = ", stencily_c)
println("stencily_r = ", stencily_r)
print("\n")

time = @elapsed for tstep=1:nStep  # loop over timesteps
# advance from timestep tstep to tstep + 1

  println("\ntstep = ", tstep)
#  print("\n")
  # it shouldn't be necessary to apply the BC at subsequent time steps
  # because the BC is built into the matrix A
  # print verification here
#  uL = BCL( (tstep-1)*delta_t )
#  println("uL = ", uL, "u[1] = ", u_i[1])
#  ghost_val = 2*delta_x*BCR( (tstep-1)*delta_t ) + u_i[mat_size - 2]
#  println("ghost_val = ", ghost_val, " u[mat_size] = ", u_i[mat_size])


  # x step
  println("\nPerforming x step")
  t_val = (tstep-0.5)*delta_t  # time at n+1/2
  t_qval = (tstep-0.75)*delta_t # time at n+1/4

  for j=2:(mat_sizey-1)  # loop over x grid lines, not including ghost lines
#    println("x grid line ", j-1)
    y_j = (j-2)*delta_y
#    println("y_j = ", y_j)

    # apply BC to rhs
    rhs_x[1] = BCxl(y_j, t_val)
    rhs_x[mat_sizex] = BCxr(y_j, t_val)

    # calculate rhs interior points
    for i=2:(mat_sizex-1)
      x_i = (i-1)*delta_x
      u_k = u_i[j, i]
      u_k_1 = u_i[j - 1, i]
      u_k_p1 = u_i[j + 1, i]
      src_val = source( x_i, y_j, t_qval) 
      rhs_x[i] = stencilx_l*u_k_1 + stencilx_c*u_k  + stencilx_r*u_k_p1 + delta_t*src_val/2      
    end

#    println("u[(j-1):(j+1), :] = \n", u_i[(j-1):(j+1), :])
#    println("rhs_x= \n", rhs_x)
    # solve for next time step
  #  A_ldiv_B!(Af, rhs)  # rhs gets overwritten with new solution values
    u_i_1[j, :] = Axf\rhs_x  # check that this actually works

#    println("underlying array = \n", u_i_1[j, :])


  end  # end loop over x grid lines

  # apply the ghost cell BC at ymin, ymax
  for i=1:mat_sizex
    x_i = (i-1)*delta_x
    u_i_1[1, i] = -2*delta_x*BCyl(x_i, t_val) + u_i_1[3, i]
    u_i_1[mat_sizey, i] = 2*delta_x*BCyr(x_i, t_val) + u_i_1[mat_sizex - 2, i]
  end

#  u_tmp = u_i_1[ 2:(end-1), :]

#  println("after x step, u_i_1 = \n", u_i_1)
#  calcError(u_tmp, xmin, xmax, ymin, ymax, t_val, N)

  # update variables
  tmp = u_i
  u_i = u_i_1
  u_i_1 = tmp

  println("\nPerforming y step")

  # loop over y grid lines
  # skip first, last gridlines because they are determined by Dirchelet BC
  # manually apply ghost point BC (although this might be unnecessary)

  t_val = (tstep)*delta_t  # time at n+ 1
  t_qval = (tstep-0.25)*delta_t # time at n+3/4


  for j=2:(mat_sizex-1)
#    println("y grid line ", j)
    x_j = (j-1)*delta_x

    # apply BC to rhs
    rhs_y[1] = BCyl(x_j, t_val)
    rhs_y[mat_sizey] = BCyr(x_j, t_val)

    # calculate rhs interior points (not including ghost points)
    for i=2:(mat_sizey-1)
      y_i = (i-2)*delta_y
      u_k = u_i[i, j]
      u_k_1 = u_i[i, j - 1]
      u_k_p1 = u_i[i, j + 1]
      src_val = source( x_j, y_i, t_qval) 
      rhs_y[i] = stencily_l*u_k_1 + stencily_c*u_k  + stencily_r*u_k_p1 + delta_t*src_val/2      
    end

#    println("(j-1):(j+1) = ", (j-1):(j+1))
#    println("u[:, (j-1):(j+1)] = \n", u_i[:, (j-1):(j+1)])
#    println("rhs_y= \n", rhs_y)
 
    u_i_1[:, j] = Ayf\rhs_y

#    println("underlying array = \n", u_i_1[:, j])
  end  # end loop over y grid lines


  # apply the dirchilet BC at xmin, xmax
  for i=2:(mat_sizey - 1)
    y_j = (i-2)*delta_y
    u_i_1[i, 1] = BCxl(y_j, t_val)
    u_i_1[i, mat_sizex] = BCxr(y_j, t_val)
  end

  # apply the ghost cell BC at ymin, ymax
  for i=1:mat_sizex
    x_i = (i-1)*delta_x
    u_i_1[1, i] = -2*delta_x*BCyl(x_i, t_val) + u_i_1[3, i]
    u_i_1[mat_sizey, i] = 2*delta_x*BCyr(x_i, t_val) + u_i_1[mat_sizey - 2, i]
  end


#  u_tmp = u_i_1[ 2:(end-1), :]

#  println("after y step, u_i_1 = \n", u_i_1)
#  calcError(u_tmp, xmin, xmax, ymin, ymax, t_val, N)



  # update variables
  tmp = u_i
  u_i = u_i_1
  u_i_1 = tmp



end  # end time step loop

# apply BCs to final time
#u_i_1[1] = BCL((nStep - 1)*delta_t)
#u_i_1[N] = BCR((nStep - 1)*delta_t)


println("time = ", time)

return u_i, delta_t*(nStep)  # plus 1 because we are at the beginning of the next timestep

end


function IC1(x, y)
  return sin(x)*cos(y) - 3*sin(x)*cos(2*y)
end


function BC1(y, t)
  return 0
end

function BC2(y, t)
  return 0
end

function BC3(x, t)
  return 0
end

function BC4(x, t)
  return 0
end


function SRC(x, y, t)
  return 0
end

function uExact(x, y, t)
  return sin(x)*cos(y)*e^(-2*t) - 3*sin(x)*cos(2*y)e^(-5*t)
end



#=
function IC1(x, y)
  return x*x + x + y*y +y + 1
end


function BC1(y, t)
  return y*y + y + t*t + t + 1
end

function BC2(y, t)
  return 2 + y*y + y + t*t + t + 1
end

function BC3(x, t)
  return 1
end

function BC4(x, t)
  return 3
end


function SRC(x, y, t)
  return 2*t -3
end

function uExact(x, y, t)
  return x*x + x + y*y + y + t*t + t + 1
end
=#


# run
driver()
