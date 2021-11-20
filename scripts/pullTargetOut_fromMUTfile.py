import sys, os, glob

if __name__ == "__main__":

	path2RelateResult = sys.argv[1]

	list_subDirs = glob.glob(path2RelateResult + "/*/*")

	for p in list_subDirs:
		region = p.split("/")[-1].split(".")[0]
		chrom, start, end = region.split("_")
		pos = (int(start) + int(end)) / float(2)
		mutFile = "%s/%s.relateOut.mut" % (p, region)
		with open(mutFile) as fin:
			for line in fin:
				if line.startswith("snp"):
					continue
				newline = line.strip().split(";")
				if int(newline[1]) == int(pos):
					print (line.strip())


