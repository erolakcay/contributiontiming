### 1- Initialize parameters ...
#println("Here is your queue number: ", queue)
grid = collect(Base.product(
# lines with value = 99 are currently depracted
[0.4 0.45 0.5 0.55 0.6], #1 investment per capita threshold α
[0.1 0.2 0.3 0.4 0.5], #2 delay δ
[0.0001], #3 diffusion rate
))[queue]; # additional repetitions (global; to spread jobs over more cores)

# # Record parameters for export
# setup = Dict(
#   # Setup for simulation
#   "nRound"            => grid[1], # ...
#   );
#
# # Make parameters globally accessable
# nRound            = setup["nRound"];


## Set local parameters
global α = grid[1] # investment per capita threshold
global δ = grid[2] # delay
global D = grid[3] # diffusion rate

global maxtime = 10.0
global numberofsims = 2
global sidelength = 4
global carryingcapacity = 1000.0
