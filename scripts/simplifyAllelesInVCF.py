import sys, gzip, random

if __name__ == "__main__":
	
	if sys.argv[1] == "-":
		fin = sys.stdin
	else:
		fin = open(sys.argv[1])

	for line in fin:
		if line.startswith("#"):
			print(line.strip())
		else:
			newINFO = None
			flagBreak = False
			newline = line.strip().split()

			if len(newline[3]) > 1:
				if "," in newline[3]:
					newline[3] = "INDEL"
					print("\t".join(newline))
					flagBreak = True
				else:
					newINFO = "origREF_%s" % newline[3] 
					newline[3] = newline[3][0]

			if len(newline[4]) > 1:
				if "," in newline[4]:
					newline[4] = "INDEL"
					print("\t".join(newline))
					flagBreak = True
				else:
					if newline[4][0] == newline[3][0]:
						letters = ["A","T","C","G"]
						letters.remove(newline[3][0])
						rand_nu = random.sample(letters, 1)[0]
						if newINFO is None:
							newINFO = "origALT_%s" % newline[4]
						else:
							newINFO = newINFO + "_origALT_%s" % newline[4]
						newline[4] = rand_nu
			else:
				if newINFO is not None:
					if newline[4][0] == newline[3][0]:
						letters = ["A","T","C","G"]
						letters.remove(newline[3][0])
						rand_nu = random.sample(letters, 1)[0]
						if newINFO is None:
							newINFO = "origALT_%s" % newline[4]
						else:
							newINFO = newINFO + "_origALT_%s" % newline[4]
						newline[4] = rand_nu

			if not flagBreak:
				if newINFO is not None:
					newline[7] = newline[7] + ";" + newINFO
					print("\t".join(newline))
				else:
					print(line.strip())


