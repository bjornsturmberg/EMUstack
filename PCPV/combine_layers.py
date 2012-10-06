

import numpy as np

def load_scat_mat(name, st, p):
    # reshape matrices to be consistent with pcpv.exe output
    format_title = '%04d' % st
    format_p     = '%04d' % p

    file_name = "st%(st)s_wl%(wl)s_%(mat_name)s.txt" % {
        'st' : format_title, 'wl' : format_p, 'mat_name' : name }
    data   = np.loadtxt(file_name)
    num_1  = max(data[:,0])
    num_2  = max(data[:,1])
    matrix = np.mat(data[:,2] + data[:,3]*(0+1j))
    matrix = np.reshape(matrix, (num_2, num_1))
    return matrix


def rmat_through_air(top_layer, bot_layer, wl):
    # P -> I for infinitesimal air layer
    r10 = load_scat_mat('R21', top_layer.label_nu, wl).T
    t01 = load_scat_mat('T12', top_layer.label_nu, wl).T
    r01 = load_scat_mat('R12', top_layer.label_nu, wl).T
    r02 = load_scat_mat('R12', bot_layer.label_nu, wl).T
    t10 = load_scat_mat('T21', top_layer.label_nu, wl).T
    I_mat    = np.matrix(np.eye(len(r01)),dtype='D')
    inv_term = (I_mat - r01*r02).I
    rmat = r10 + t01*r02*inv_term*t10
    return rmat

def tmat_through_air(top_layer, bot_layer, wl):
    # P -> I for infinitesimal air layer
    t02 = load_scat_mat('T12', bot_layer.label_nu, wl).T
    r01 = load_scat_mat('R12', top_layer.label_nu, wl).T
    r02 = load_scat_mat('R12', bot_layer.label_nu, wl).T
    t10 = load_scat_mat('T21', top_layer.label_nu, wl).T
    I_mat    = np.matrix(np.eye(len(r01)),dtype='D')
    inv_term = (I_mat - r01*r02).I
    tmat = t02*inv_term*t10
    return tmat


def save_TR_mats(matrix, name, p):
            # reshape matrices to be consistent with pcpv.exe output
            format_p     = '%04d' % p

            file_name = "wl%(wl)s_%(mat_name)s.txt" % {
                'wl' : format_p, 'mat_name' : name }
            num_pw = len(matrix)
            with file(file_name, 'w') as outfile:
                for k in range(num_pw):
                    for i in range(num_pw):
                        data = [i+1,  k+1, np.real(matrix[i,k]), np.imag(matrix[i,k]),
                            np.abs(matrix[i,k])**2]
                        data = np.reshape(data, (1,5))
                        np.savetxt(outfile, data, fmt=['%4i','%4i','%25.17f','%25.17f','%25.17f'], delimiter='')



def net_scat_mats(solar_cell, nu_wls):
    for p in range(nu_wls):
        p += 1

        # calculate list of scattering matrices of interfaces without air slices\
        # notice counting upwards from substrate!
        # ie r12s[0] = refelction from layer 0 off layer 1
        r12s = []
        r21s = []
        t12s = []
        t21s = []
        for st1, st2 in zip(solar_cell[:-1], solar_cell[1:]):
            r12s.append(rmat_through_air(st1, st2, p))
            r21s.append(rmat_through_air(st2, st1, p))
            t12s.append(tmat_through_air(st1, st2, p))
            t21s.append(tmat_through_air(st2, st1, p))#t12s[-1].T)

        # iterate through solar cell to find total Tnet, Rnet matrices
        rnet = r21s[0]
        tnet = t21s[0]
        for i in range(1,len(solar_cell)-1):
            P = load_scat_mat('P', solar_cell[i].label_nu, p).T
            I_mat   = np.matrix(np.eye(len(P)),dtype='D')
            inverse_term = (I_mat - r12s[i]*P*rnet*P).I
            repeated_term = P*inverse_term*t21s[i]
            tnet = tnet*repeated_term
            rnet = r21s[i] + t12s[i]*P*rnet*repeated_term

        save_TR_mats(tnet, 'Tnet', p)
        save_TR_mats(rnet, 'Rnet', p)


    # NEED to implement more involved selections from Tnet Rnet matrices, 
    #      eg for Left, Right pol and CD
    # atm inc in and all prop out
    # also want absorption in each layer etc

        # select elements of Tnet, Rnet matrices to calculate absorption etc.
        neq_PW   = solar_cell[0].nu_tot_ords
        select_ord_in  = solar_cell[-1].set_ord_in
        select_ord_out = solar_cell[0].set_ord_out

        Lambda_t = 0
        inc      = solar_cell[0].zero_ord
        # tnet = load_scat_mat('T12', solar_cell[0].label_nu, p)#.T
        for i in range(solar_cell[0].nu_prop_ords):
            #TM
            # Lambda_t = Lambda_t + abs(tnet[i,neq_PW+inc])**2 + abs(tnet[neq_PW+i,neq_PW+inc])**2
            #TE
            Lambda_t = Lambda_t + abs(tnet[i,inc])**2 + abs(tnet[neq_PW+i,inc])**2

        Lambda_r = 0
        inc      = solar_cell[-1].zero_ord
        for i in range(solar_cell[-1].nu_prop_ords):
            #TM
            # Lambda_r = Lambda_r + abs(rnet[i,neq_PW+inc])**2 + abs(rnet[neq_PW+i,neq_PW+inc])**2
            #TE
            Lambda_r = Lambda_r + abs(rnet[i,inc])**2 + abs(rnet[neq_PW+i,inc])**2

        
        absorption = 1 - Lambda_r - Lambda_t
        print 'r', Lambda_r
        print 't', Lambda_t
        print 'a', absorption

        # absorption = 1 - Lambda_r - Lambda_t*(solar_cell[0]/solar_cell[-1])
        # print 'a', absorption

