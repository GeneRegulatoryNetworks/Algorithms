local FD = {
	tables = {},
	U = require("Utility"),
	functions = {},
	metatable = {}
}

local insert, concat, sort = table.insert, table.concat, table.sort
local create, yield = coroutine.create, coroutine.yield

local L = FD.U.List
local M = FD.U.Matrix

-- Initialize logarithm lookup tables.
FD.tables.logs = L.new(function(x) return -math.log(x) end, 32)

FD.tables.factorial_logs = L.new(create(function(length)
	local acc = 0.0

	yield()

	for i = 1, length do
		acc = acc - math.log(i)

		yield(acc)
	end
end), 32)

-- Returns a functional digraph (with no properties computed) assuming
-- update is a list representing the images of [1, N] under an endomorphism
-- for some N <= 32.
function FD.new(update)
	local N = #update

	assert(N <= 32, "Length of update list must be <= 32.")

	return setmetatable({
		N = N,
		update = L.new(update),
		complete = false
		--[[
			Other fields:
				preimages
				adjacency
				transition (optional)
				cycles (optional)
				components
				trees
				classes
				entropy = -ln |Aut(G)|
				string (optional)
		]]
	}, FD.metatable)
end

-- Computes the preimage table of the update function.
function FD.functions:computePreimages()
	local preimages = {}

	for i, v in ipairs(self.update) do
		local preimage = preimages[v]

		if preimage then
			insert(preimage, i)
		else
			preimages[v] = { i }
		end
	end

	self.preimages = preimages
end

-- Utility function to build a recursive tree structure via the preimage
-- table. A tree vertex looks like { children = array }, where array is a
-- (possibly empty) list of vertices. The identity field is used for
-- bookkeeping.
function FD.functions:buildTree(root)
	local children = self.preimages[root.identity]

	if not children then
		return root
	end

	for i, v in ipairs(children) do
		root.children[i] = self:buildTree({ children = {}, identity = v })
	end

	return root
end

-- Computes the canonical number and entropy of a tree simultaneously. See
-- paper for algorithm description.
function FD.functions:canonicalNumber(root)
	local length = #root.children

	-- Assign 10 to leaves.
	if length == 0 then
		return {
			number = 2,
			entropy = 0.0,
			length = 2 -- Number of bits in canonical number.
		}
	end

	local data = {}

	-- Apply computation to children. Recursion is neat :)
	for i, v in ipairs(root.children) do
		data[i] = self:canonicalNumber(v)
	end

	-- Sort labels.
	sort(data, function(d1, d2) return d1.number < d2.number end)

	-- Compute entropy.
	local entropy = 0.0
	local count = 1
	local current = data[1]

	for i = 2, length+1 do
		if i == length+1 or current.number ~= data[i].number then
			entropy =
				entropy +
				count * current.entropy +
				FD.tables.factorial_logs[count]
			current = data[i]
			count = 1
			i = i + 1
		else
			count = count + 1
		end
	end

	local number = 0
	local pos = 0

	-- Concatenate labels in descending order.
	for i, v in ipairs(data) do
		number = number | (v.number << pos)
		pos = pos + v.length
	end

	assert(pos <= 62, "Tree must have <= 32 nodes.")

	return {
		number = (1 << (pos + 1)) | (number << 1),
		entropy = entropy,
		length = pos + 2
	}
end

-- Computes the adjacency matrix of the digraph.
function FD.functions:computeAdjacency()
	self.adjacency = M.new(create(function(N)
		yield()

		for i = 1, N do
			local pos = self.update[i]

			for j = 1, N do
				yield(j == pos and 1 or 0)
			end
		end
	end), self.N)
end

-- Computes the transition matrix of the digraph.
function FD.functions:computeTransition()
	self.transition = self.adjacency:transpose()
end

-- Computes the connected components of the digraph. The adjacency matrix
-- must already be computed.
function FD.functions:computeComponents()
	self.components = self.adjacency:oneEigenbasis()
end

-- Computes the cycles of the digraph. The transition matrix must already be
-- computed.
function FD.functions:computeCycles()
	self.cycles = self.transition:oneEigenbasis()
end

-- Computes tree and component isomorphisms and automorphisms. The
-- components must already be computed.
function FD.functions:computeTrees()
	local trees = {}

	for index, component in ipairs(self.components) do
		-- Find cycle states.
		local length = #component
		-- Bijection from states in component to [1, length].
		local inverse = {}

		for i, v in ipairs(component) do
			inverse[v] = i
		end

		local transition = M.new(create(function(N)
			yield()

			for i = 1, N do
				for j = 1, N do
					yield(
						i == inverse[self.update[component[j]]] and
						1 or
						0
					)
				end
			end
		end), length)

		local cycle = transition:oneEigenbasis()[1]
		local c_length = #cycle

		-- Recreate cycle in full space with correct order.
		cycle[1] = component[cycle[1]]

		for i = 2, c_length do
			cycle[i] = self.update[cycle[i-1]]
		end

		-- Extract tree structures and run canonical number algorithm.
		local numbers = L.new()
		local entropies = {}

		for i, v in ipairs(cycle) do
			local predecessor = cycle[i>1 and i-1 or c_length]
			local root = { children = {}, identity = v }

			for _, w in ipairs(self.preimages[v]) do
				if w ~= predecessor then
					insert(
						root.children,
						self:buildTree({ children = {}, identity = w })
					)
				end
			end

			local result = self:canonicalNumber(root)
			numbers[i] = result.number
			entropies[i] = result.entropy
		end

		numbers = L.new(numbers, numbers:blockLength())
		numbers = numbers:rotate(numbers:lmr())
		local block_length = #numbers
		local num_blocks = c_length // block_length

		-- Compute canonical string and entropy of component.
		for i, v in ipairs(numbers) do
			numbers[i] = string.format("%x", v)
		end

		insert(numbers, 1, num_blocks)

		local entry = { string = concat(numbers, ",") }
		local entropy = 0.0

		for i = 1, block_length do
			entropy = entropy + entropies[i]
		end

		entry.entropy = num_blocks * entropy + FD.tables.logs[num_blocks]

		trees[index] = entry
	end

	self.trees = trees
end

-- Finds component isomorphisms using the canonical strings.
function FD.functions:computeClasses()
	local classes = {}
	local entropy = 0.0

	for _, v in ipairs(self.trees) do
		local class = classes[v.string]

		if class then
			class.count = class.count + 1
		else
			classes[v.string] = {
				count = 1,
				entropy = v.entropy
			}
		end
	end

	for _, v in pairs(classes) do
		entropy =
			entropy +
			v.count * v.entropy +
			FD.tables.factorial_logs[v.count]
	end

	self.classes = classes
	self.entropy = entropy
end

-- Returns the canonical string of the digraph. Computes and stores the
-- result if currently unknown. Two digraphs have the same string iff
-- they're isomorphic.
function FD.functions:getString()
	if self.string then
		return self.string
	end

	assert(self.complete, "Must call compute() before tostring().")

	local strings = {}

	for k, v in pairs(self.classes) do
		insert(strings, concat({ v.count, k }, ","))
	end

	sort(strings)

	self.string = concat(strings, "+")

	return self.string
end

-- Computes the selected functions on the digraph.
function FD.functions:selected(list)
	for _, v in ipairs(list) do
		FD.functions[v](self)
	end
end

-- Computes non-optional properties of the digraph.
function FD.functions:compute()
	self:selected({
		"computePreimages",
		"computeAdjacency",
		"computeComponents",
		"computeTrees",
		"computeClasses"
	})

	self.complete = true
end

-- Determines whether two digraphs are isomorphic.
function FD.functions:equal(other)
	return self.N == other.N and tostring(self) == tostring(other)
end

FD.metatable = {
	__index = FD.functions,
	__tostring = FD.functions.getString,
	__eq = FD.functions.equal
}

return FD
