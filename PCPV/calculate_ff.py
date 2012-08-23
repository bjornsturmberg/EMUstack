import numpy as np

pi = np.pi
sr = np.sqrt

def calculate_ff(d, a1, a2=0, a3=0, a4=0, a5=0, a6=0, a7=0, a8=0, a9=0, a10=0,
	a11=0, a12=0, a13=0, a14=0, a15=0, a16=0, el1 = 0):

	ff = pi*(a1**2*sr(1-el1) + a2**2 + a3**2 + a4**2 + a5**2 + a6**2 + a7**2 + a8**2 + a9**2
		+ a10**2 + a11**2 + a12**2 + a13**2 + a14**2 + a15**2 + a16**2)/(d)**2
	# print "ff = ", ff
	# ff_round = round(ff,2)
	# print "ff = ", ff_round
	return ff#_round


if __name__ == "__main__":
    import sys
    calculate_ff(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
