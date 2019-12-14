import sys

suffix = ".par"
prefix = sys.argv[1]

input_name = prefix + suffix
output_name = "sort_" + prefix + suffix

res = []

with open(input_name, "r") as ifs:
	for line in ifs:
		i_line = int(line)
		res.append(i_line)

res.sort()

with open(output_name, "w") as ofs:
	for r in res:
		ofs.write(str(r))
		ofs.write("\n")
