import numpy as np

pi = np.pi

def calculate_ff(d, a1, a2, a3, a4):

	ff = pi*(a1**2 + a2**2 + a3**2 + a4**2)/(d)**2
	# print "ff = ", ff
	# ff_round = round(ff,2)
	# print "ff = ", ff_round
	return ff#_round


if __name__ == "__main__":
    import sys
    calculate_ff(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
