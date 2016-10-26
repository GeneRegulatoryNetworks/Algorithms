local GRN = {
	tables = {},
	FD = require("FD"),
	functions = {},
	metatable = {}
}

local L = GRN.FD.U.List
local M = GRN.FD.U.Matrix

-- Initialize lookup tables.
do
	local weights = L.new(0, 1)

	for i = 1, 5 do
		local length = #weights

		for j = 1, length do
			weights[length+j] = weights[j] + 1
		end
	end

	GRN.tables.weights = weights

	local genes = L.new()

	for i = 1, 5 do
		genes[1<<i] = i
	end

	GRN.tables.genes = genes
end

function GRN.functions:equal(other)
	return self.fd == other.fd
end

function GRN.functions:getString()
	return self.fd:getString()
end

GRN.metatable = {
	__index = GRN.functions,
	__tostring = GRN.functions.getString,
	__eq = GRN.functions.equal
}

-- Returns a GRN (with no properties computed) given the underlying update
-- function. The number of states must arise from 1-5 genes.
function GRN.new(update)
	local length = #update

	assert(length >= 2 and length <= 32, "Number of states out of range.")

	assert(length & (length - 1) == 0, "Number of states not power of 2.")

	return setmetatable({
		fd = GRN.FD.new(update),
		genes = GRN.tables.genes[length]
		--[[
			Other fields:
			stabilities
			entropy
		]]
	}, GRN.metatable)
end

-- Computes T values for each component of the network.
function GRN.functions:computeStabilities()
	local stabilities = {}

	for index, component in ipairs(self.fd.components) do
		local length = #component
		local spectrum = L.new(0, self.genes+1)
		spectrum[1] = length

		for i = 1, length-1 do
			for j = i+1, length do
				local bin =
					GRN.tables.weights[
						((component[i] - 1) ~
						(component[j] - 1)) +
						1
					] + 1
				spectrum[bin] = spectrum[bin] + 2
			end
		end

		stabilities[index] = spectrum
	end

	self.stabilities = stabilities
end

-- Computes properties of the GRN.
function GRN.functions:compute()
	self.fd:compute()

	self.entropy = self.fd.entropy

	self:computeStabilities()
end

-- Returns the stability of the network given error rate e.
function GRN.functions:getStability(e)
	local stability = 0.0
	local n = self.genes
	local nu = #self.fd.components
	local q = 1 - e

	for k = 0, n do
		local combinations = 0

		for i = 1, nu do
			combinations = combinations + self.stabilities[i][k+1]
		end

		stability = stability + e ^ k * q ^ (n - k) * combinations
	end

	return stability / self.fd.N
end

return GRN
