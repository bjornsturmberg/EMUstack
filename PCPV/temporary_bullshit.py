import numpy as np


def save_scat_mat(matrix, name, st, p, num_pw):
    # reshape matrices to be consistent with pcpv.exe output
    format_label_nu = '%04d' % st
    format_p        = '%04d' % p

    file_name = "st%(st)s_wl%(wl)s_%(mat_name)s.txt" % {
        'st' : format_label_nu, 'wl' : format_p, 'mat_name' : name }
    with open(file_name, 'w') as outfile:
        for k in range(num_pw):
            for i in range(num_pw):
                data = [i+1,  k+1, np.real(matrix[i,k]), np.imag(matrix[i,k]),
                    np.abs(matrix[i,k])**2]
                data = np.reshape(data, (1,5))
                np.savetxt(outfile, data, fmt=['%4i','%4i','%25.17G','%25.17G','%25.17G'], delimiter='')

def save_k_perps(anallo_list, num_pw):
    data_out = np.zeros((len(anallo_list), 2 + 2*len(anallo_list[0].beta)))

    #TODO: MAKE THIS SHIT REDUNDANT
    # Failing that,
    #TODO: check that beta is the same length for everything
    #TODO: check that everything is the same label_nu and num_pw

    for i, an in enumerate(anallo_list):
        data_out[i,:2] = (num_pw, an.light.Lambda)
        re, im = an.beta.real, an.beta.imag
        # beta = [beta[0].real, beta[0].imag, beta[1].real, ...]
        beta = np.vstack((re, im)).T.reshape(-1)
        data_out[i,2:] = beta

    filename = "beta_st%04d.txt" % an.thin_film.label_nu
    np.savetxt(filename, data_out, fmt='%25.17G', delimiter='')

def save_omegas(simmo_list):
    pass