-- This file demonstrates how to use the GRN library.

local G = require("GRN")

local n = 5

local N = 1 << n
local a = -G.FD.tables.factorial_logs[N]
local b = a + N * G.FD.tables.logs[N]

math.randomseed(os.time())

local function out(str, ...)
	print(string.format(str, ...))
end

-- Create a random n-gene update function.
local update = G.FD.U.List.new(function() return math.random(1, N) end, N)

print("Update function:")

for i, v in ipairs(update) do
	out("%d -> %d", i, v)
end

-- Initialize a new GRN.
local g = G.new(update)

-- Compute essential properties.
g:compute()

-- Print out some properties.
out("This GRN has %d connected component(s).", #g.fd.components)

for i, v in ipairs(g.fd.components) do
	out("Component %d contains states: %s.", i, table.concat(v, ", "))
end

out("It belongs to the isomorphism class represented by str(C) = %s.", g)

out("The approximate entropy is %.3f.", g.entropy)

out("This means that ln |C| = %.3f.", a+g.entropy)

out("For f chosen uniformly at random, ln(P(f in C)) = %.3f.", b+g.entropy)

out("It has a stability of %.3f (epsilon = 1/%d).", g:getStability(1/n), n)
