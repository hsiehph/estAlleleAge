import sys


if __name__ == "__main__":

	with open(sys.argv[1]) as fin:
		prev = None
		for line in fin:
			newline = line.strip().split()
			if prev is None:
				prev = newline
				continue

			diff = float(newline[2]) - float(prev[2])
			if diff < 1e-5:
				continue
			else:
				print ("\t".join(prev))
				prev = newline



