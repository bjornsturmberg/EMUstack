import numpy as np

pi = np.pi

def calculate_ff(a1, a2, d):

	# print "a1 = ", a1
	# print "a2 = ", a2
	# print "d  = ", d

	ff = 2*pi*(a1**2 + a2**2)/(d)**2
	# print "ff = ", ff
	ff_round = round(ff,2)
	# print "ff = ", ff_round
	return ff_round


if __name__ == "__main__":
    import sys
    calculate_ff(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
