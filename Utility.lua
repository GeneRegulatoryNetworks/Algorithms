-- Disable this assertion at your own risk :P
assert(_VERSION == "Lua 5.3", "This library is written for Lua 5.3.")

assert(1 << 63 ~= 0, "This library requires an integer width of >= 64.")

local Utility = {
	List = {
		functions = {},
		metatable = {}
	},
	Matrix = {
		functions = {},
		metatable = {}
	}
}

local insert, concat = table.insert, table.concat
local resume = coroutine.resume

-- Returns a string representation of the list.
function Utility.List.functions:display()
	return concat(self, ", ")
end

Utility.List.metatable = {
	__index = Utility.List.functions,
	__tostring = Utility.List.functions.display
}

-- Creates a list from initializer according to its type:
---- nil: empty list.
---- number, boolean, string: *length* copies of initializer.
---- table: duplicates entries of initializer up to *length*.
------ if length == nil or length > #initializer, it is entirely duplicated.
---- function: list[i] = initializer(i).
---- coroutine: resumes coroutine with *length* argument and takes
------ subsequent return values as list elements, up to the first error.
function Utility.List.new(initializer, length)
	local t = type(initializer)
	local new = {}

	if t == "nil" then
		-- Do nothing.
	elseif t == "number" or t == "boolean" or t == "string" then
		for i = 1, length do
			new[i] = initializer
		end
	elseif t == "table" then
		length = math.min(#initializer, length or math.huge)

		for i = 1, length do
			new[i] = initializer[i]
		end
	elseif t == "function" then
		for i = 1, length do
			new[i] = initializer(i)
		end
	elseif t == "thread" then
		resume(initializer, length)

		for i = 1, length do
			local status, value = resume(initializer)

			if status then
				new[i] = value
			else
				break
			end
		end
	else
		error("Invalid list initializer.")
	end

	return setmetatable(new, Utility.List.metatable)
end

-- Duplicates the list.
function Utility.List.functions:dup()
	return Utility.List.new(self)
end

-- Returns a rotated list beginning at index.
function Utility.List.functions:rotate(index)
	local length = #self
	local new = Utility.List.new()

	for i = 1, length do
		new[i] = self[(i+index-2)%length+1]
	end

	return new
end

-- Returns the length of the smallest repeating sublist (prime block).
-- This block is unique and aperiodic (i.e., every rotation is unique).
-- When C_L acts on lists of length L by rotation, the stabilizer of a list
-- l has order L / blockLength(l).
-- Uses ideas from the KMP algorithm.
function Utility.List.functions:blockLength()
	local length = #self
	local nxt = Utility.List.new(0, length)

	for i = 1, length-1 do
		local k = nxt[i]

		while true do
			if self[i+1] == self[k+1] then
				nxt[i+1] = k + 1

				break
			elseif k == 0 then
				nxt[i+1] = 0

				break
			else
				k = nxt[k]
			end
		end
	end

	local block_length = length - nxt[#nxt]

	return length % block_length == 0 and block_length or length
end

-- Returns the Lyndon factorization of the list.
function Utility.List.functions:factor()
	local factorization = Utility.List.new()
	local k = 0
	local length = #self

	while k < length do
		local i = k + 1
		local j = k + 2

		while j <= length do
			local ai = self[i]
			local aj = self[j]

			if ai > aj then
				break
			end

			i = ai == aj and i + 1 or k + 1
			j = j + 1
		end

		local difference = j - i

		repeat
			k = k + difference

			insert(factorization, k)
		until k >= i
	end

	return factorization
end

-- Outputs the starting point of the lexicographically minimal rotation of
-- the list. Two lists are circularly equivalent iff they have the same LMR.
-- Note that rotating a list is equivalent to rotating its prime block by
-- the same distance, so it's advisable to use blockLength() first.
function Utility.List.functions:lmr()
	local factorization = self:factor()
	local f_length = #factorization

	if f_length == 1 then
		return 1
	end

	-- offset is the starting position of the final Lyndon word.
	local offset = factorization[f_length-1] + 1
	local difference = #self - offset

	-- This loop goes backwards until it finds the first consecutive
	-- instance of the final Lyndon word. It doesn't check the first word,
	-- which isn't a problem; if all the words are equal, then it doesn't
	-- matter which one you pick.
	for i = f_length-2, 1, -1 do
		local f_begin = factorization[i] + 1
		local f_end = factorization[i+1]

		if f_end - f_begin ~= difference then
			return offset
		end

		for j = 0, difference do
			if self[f_begin] ~= self[offset+j] then
				return offset
			end

			f_begin = f_begin + 1
		end

		offset = f_end - difference
	end

	return offset
end

-- Adds two matrices of the same size.
function Utility.Matrix.functions:add(other)
	local N = #self

	assert(N == #other, "Mismatched matrices.")

	local new = {}

	for i = 1, N do
		local row = {}
		local my_row = self[i]
		local other_row = other[i]

		for j = 1, N do
			row[j] = my_row[j] + other_row[j]
		end

		new[i] = row
	end

	return setmetatable(new, Utility.Matrix.metatable)
end

-- Returns a string representation of the matrix.
function Utility.Matrix.functions:display()
	local representation = {}

	for i, row in ipairs(self) do
		representation[i] = concat(row, "\t")
	end

	return concat(representation, "\n")
end

Utility.Matrix.metatable = {
	__index = Utility.Matrix.functions,
	__add = Utility.Matrix.functions.add,
	__tostring = Utility.Matrix.functions.display
}

-- Creates a matrix from initializer from its type.
---- number: the scalar matrix initializer*I.
---- matrix: duplicates initializer. N is ignored.
---- function: matrix[i][j] = initializer(i, j).
---- coroutine: passes N and uses return values in row-major order. Errors
------ are propagated.
function Utility.Matrix.new(initializer, N)
	local t = type(initializer)
	local new = {}

	if t == "number" then
		for i = 1, N do
			local row = {}

			for j = 1, N do
				row[j] = i == j and initializer or 0
			end

			new[i] = row
		end
	elseif t == "table" then
		for i, old_row in ipairs(initializer) do
			local row = {}

			for j, v in ipairs(old_row) do
				row[i] = v
			end

			new[i] = row
		end
	elseif t == "function" then
		for i = 1, N do
			local row = {}

			for j = 1, N do
				row[j] = initializer(i, j)
			end

			new[i] = row
		end
	elseif t == "thread" then
		resume(initializer, N)

		for i = 1, N do
			local row = {}

			for j = 1, N do
				local status, value = resume(initializer)

				assert(status, value)

				row[j] = value
			end

			new[i] = row
		end
	else
		error("Invalid matrix initializer.")
	end

	return setmetatable(new, Utility.Matrix.metatable)
end

-- Duplicates the matrix.
function Utility.Matrix.functions:dup()
	return Utility.Matrix.new(self)
end

-- Returns the matrix transpose.
function Utility.Matrix.functions:transpose()
	local N = #self
	local new = {}

	for i = 1, N do
		local row = {}

		for j = 1, N do
			row[j] = self[j][i]
		end

		new[i] = row
	end

	return setmetatable(new, Utility.Matrix.metatable)
end

-- Computes the rref of the matrix in place and returns the matrix.
-- IMPORTANT: only use on matrices of the form A - I where A's columns are
-- selected from the standard basis vectors (with replacement) or the
-- transpose of such a matrix.
function Utility.Matrix.functions:rref()
	local N = #self
	local j = 1

	for i = 1, N do
		-- Find a nonzero column and produce a pivot.
		for k = j, N do
			j = k

			-- if A_ij is zero, swap in a nonzero pivot.
			if self[i][k] == 0 then
				for l = i+1, N do
					if self[l][k] ~= 0 then
						self[i], self[l] = self[l], self[i]

						break
					end
				end

				if self[i][k] ~= 0 then
					break
				end
			else
				break
			end
		end

		-- If the pivot is -1, multiply row i by -1.
		if self[i][j] == -1 then
			local row = self[i]

			for k = j, N do
				row[k] = -row[k]
			end
		end

		-- Eliminate other nonzero entries in the column.
		for k = 1, N do
			local entry = self[k][j]

			if entry ~= 0 and k ~= i then
				local row_i = self[i]
				local row_k = self[k]

				if entry == -1 then
					-- If A_kj = -1, add row i to row k.
					for l = j, N do
						row_k[l] = row_k[l] + row_i[l]
					end

					-- If A_kj = 1, subtract row i from row k.
				else
					for l = j, N do
						row_k[l] = row_k[l] - row_i[l]
					end
				end
			end
		end

		j = j + 1

		if j > N then
			break
		end
	end

	return self
end

-- Returns a list of lists that represents a basis for the 1-eigenspace of
-- the matrix. Each basis vector is the sum of some subset of the standard
-- basis vectors. Each list in the return value contains the corresponding
-- indices.
-- IMPORTANT: only works on the A matrices described above.
function Utility.Matrix.functions:oneEigenbasis()
	local N = #self
	local matrix = self + Utility.Matrix.new(-1, N)
	local basis = {}
	local ones = {}

	matrix:rref()

	for j = 1, N do
		for i = 1, j do
			local entry = matrix[i][j]

			if entry == 1 then
				-- If we see a 1 in row i, keep track of the column number.
				ones[i] = j

				break
			elseif entry == -1 then
				-- If we see a -1 column, recall the positions of the 1's.
				local vector = { ones[i] }

				for k = i+1, j do
					if matrix[k][j] == -1 then
						insert(vector, ones[k])
					end
				end

				insert(vector, j)

				insert(basis, vector)

				break
			elseif i == j and entry == 0 then
				-- A zero column gives us a basis vector by itself.
				insert(basis, { j })

				break
			end
		end
	end

	return basis
end

return Utility
