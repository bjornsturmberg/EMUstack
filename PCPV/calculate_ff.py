import numpy as np

pi = np.pi

def calculate_ff(d, a1, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0, a8=0, a9=0):

	ff = pi*(a1**2 + a2**2 + a3**2 + a4**2 + a5**2 + a6**2 + a7**2 + a8**2 + a9**2)/(d)**2
	# print "ff = ", ff
	# ff_round = round(ff,2)
	# print "ff = ", ff_round
	return ff#_round


if __name__ == "__main__":
    import sys
    calculate_ff(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
