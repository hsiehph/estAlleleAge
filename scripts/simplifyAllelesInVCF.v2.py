import sys, gzip, random

if __name__ == "__main__":
	
	if sys.argv[1] == "-":
		fin = sys.stdin
	else:
		fin = open(sys.argv[1])

	with open(sys.argv[2]) as ff:
		for line in ff:
			if line.startswith(">"):
				continue
			ancestralSeq = line.strip()


	for line in fin:
		if line.startswith("#"):
			print(line.strip())
		else:
			newINFO = None
			flagBreak = False
			newline = line.strip().split()
			ancestral_state = ancestralSeq[int(newline[1])]

			if len(newline[3]) > 1 or len(newline[4]) > 1:
				# skip sites with complex, multi-allelic INDELs
				if "," in newline[3] or "," in newline[4]:
					newINFO = "origREF_%s;origALT_%s;AA:%s" % (newline[3], newline[4], ancestral_state)
					newline[7] = newline[7] + ";" + newINFO
					newline[3] = "mINDEL"
					newline[4] = "mINDEL"
					print ("\t".join(newline))
					continue

				# work on simple biallelic INDELs
				else:
					if ancestral_state == newline[3][0] or ancestral_state == newline[4][0]:
						newINFO = "origREF_%s;origALT_%s;AA:%s" % (newline[3], newline[4], ancestral_state)
						newline[7] = newline[7] + ";" + newINFO
						if newline[3][0] == newline[4][0] or newline[3][0] == ancestral_state:
							newline[3] = ancestral_state
							letters = ["A","T","C","G"]
							letters.remove(ancestral_state)
							rand_nu = random.sample(letters, 1)[0]
							newline[4] = rand_nu
							print ("\t".join(newline))
						else:
							raise TypeError(newline)

			# simply print the non-INDELs sites
			else:
				newline[7] = newline[7] + ";" + "AA:%s" % ancestral_state
				print ("\t".join(newline))



