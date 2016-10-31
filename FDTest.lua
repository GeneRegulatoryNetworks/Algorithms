local FD = require("FD")

local epsilon = 0.01 -- 1%.

local function ent(x)
	return -math.log(x)
end

local function msg(...)
	return table.concat(table.pack(...))
end

local test_cases = {
	{
		update = { 1, 2 },
		entropy = ent(2),
		string = "2,1,2",
	},
	{
		update = { 1, 1 },
		entropy = ent(1),
		string = "1,1,c",
	},
	{
		update = { 2, 1 },
		entropy = ent(2),
		string = "1,2,2",
	},
	{
		update = { 1, 2, 3, 4 },
		entropy = ent(24),
		string = "4,1,2"
	},
	{
		update = { 1, 1, 3, 4 },
		entropy = ent(2),
		string = "1,1,c+2,1,2"
	},
	{
		update = { 1, 1, 1, 4 },
		entropy = ent(2),
		string = "1,1,2+1,1,34"
	},
	{
		update = { 1, 1, 3, 3 },
		entropy = ent(2),
		string = "2,1,c"
	},
	{
		update = { 3, 3, 4, 4 },
		entropy = ent(2),
		string = "1,1,e8"
	},
	{
		update = { 2, 3, 3, 4 },
		entropy = ent(1),
		string = "1,1,2+1,1,38"
	},
	{
		update = { 2, 3, 4, 4 },
		entropy = ent(1),
		string = "1,1,f0"
	},
	{
		update = { 2, 4, 4, 4 },
		entropy = ent(1),
		string = "1,1,e4"
	},
	{
		update = { 1, 1, 1, 1 },
		entropy = ent(6),
		string = "1,1,d4"
	},
	{
		update = { 2, 1, 3, 4 },
		entropy = ent(4),
		string = "1,2,2+2,1,2"
	},
	{
		update = { 2, 1, 4, 3 },
		entropy = ent(8),
		string = "2,2,2"
	},
	{
		update = { 3, 3, 4, 3 },
		entropy = ent(2),
		string = "1,1,2,34"
	},
	{
		update = { 2, 3, 2, 3 },
		entropy = ent(2),
		string = "1,2,c"
	},
	{
		update = { 2, 3, 4, 3 },
		entropy = ent(1),
		string = "1,1,2,38"
	},
	{
		update = { 2, 3, 2, 4 },
		entropy = ent(1),
		string = "1,1,2+1,1,2,c"
	},
	{
		update = { 2, 2, 4, 3 },
		entropy = ent(2),
		string = "1,1,c+1,2,2"
	},
	{
		update = { 2, 3, 1, 4 },
		entropy = ent(3),
		string = "1,1,2+1,3,2"
	},
	{
		update = { 2, 3, 4, 2 },
		entropy = ent(1),
		string = "1,1,2,2,c"
	},
	{
		update = { 2, 3, 4, 1 },
		entropy = ent(4),
		string = "1,4,2"
	},
	{
		update = { 1, 3, 4, 5, 6, 7, 8, 5 },
		entropy = ent(1),
		string = "1,1,2+1,1,2,2,2,f0"
	},
	{
		update = { 6, 3, 2, 5, 6, 5, 3, 2 },
		entropy = ent(8),
		string = "2,2,c"
	},
	{
		update = {
			9, 16, 4, 3, 10, 12, 14, 2, 1, 13, 12, 13, 6, 16, 14, 7
		},
		entropy = ent(16),
		string = "2,1,2,c,38+2,2,2"
	},
	{
		update = { 14, 2, 9, 16, 6, 14, 11, 1, 2, 13, 2, 11, 6, 8, 16, 6 },
		entropy = ent(4),
		string = "1,1,2,2,f4c8+1,1,e98"
	},
	{
		update = {
			1, 2, 5, 10, 17, 26, 5, 18, 1, 18, 5, 26, 17, 10, 5, 2, 1, 2, 5,
			10, 17, 26, 5, 18, 1, 18, 5, 26, 17, 10, 5, 2
		},
		entropy = ent(2229534720),
		string = "1,1,f54d5294+1,1,f5554a94"
	}
}

for i, v in ipairs(test_cases) do
	local f = FD.new(v.update)
	local number = string.format("%03d", i)

	f:compute()

	if math.abs(f.entropy-v.entropy) > math.abs(epsilon*v.entropy) then
		error(msg(
			"Entropy test ",
			number,
			" failed.\n",
			"update = ",
			tostring(f.update),
			"\n",
			"expected = ",
			v.entropy,
			"\n",
			"output   = ",
			f.entropy
		))
	end

	if tostring(f) ~= v.string then
		error(msg(
			"String test ",
			number,
			" failed.\n",
			"update = ",
			tostring(f.update),
			"\n",
			"expected = ",
			v.string,
			"\n",
			"output   = ",
			tostring(f)
		))
	end

	print(msg("Test ", number, " passed."))
end
