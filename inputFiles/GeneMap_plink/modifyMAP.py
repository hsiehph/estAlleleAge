import sys


if __name__ == "__main__":

	with open(sys.argv[1]) as fin:
		l_2markers = []
		for line in fin:
			newline = line.strip().split()
			if len(l_2markers) < 2:
				l_2markers.append(newline)
			else:
				diff = float(newline[2]) - float(l_2markers[1][2])
				if diff < 1e-6:
					newGenPos_midMarker = float(newline[2]) - 1e-6
					if float(l_2markers[0][2]) < newGenPos_midMarker:
						l_2markers[1][2] = str(newGenPos_midMarker)
						print("\t".join(l_2markers[0]))
						l_2markers[0] = l_2markers[1]
						l_2markers[1] = newline
					else:
						l_2markers[1][2] = str(newGenPos_midMarker + 0.5e-6)
						print("\t".join(l_2markers[0]))
						l_2markers[0] = l_2markers[1]
						l_2markers[1] = newline

				else:
					print("\t".join(l_2markers[0]))
					l_2markers[0] = l_2markers[1]
					l_2markers[1] = newline

		for item in l_2markers:
			print("\t".join(item))



