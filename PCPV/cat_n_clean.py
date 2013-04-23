# concatination of results and clean up

import glob
import os


def c_c_tra(pol = 0):
    """ concatenate transmission, reflection and absorption results 
    in one file each and remove individual lambda results """
    if pol == 0:
        filename1 = "wl*_A_Lambda.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Lambda.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "wl*_T_Lambda.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
        if len(Afiles)+len(Rfiles)+len(Tfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> Absorptance.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> Reflectance.txt")
            listoffiles3 = ""
            for fle in Tfiles:
                listoffiles3 += " " + fle
            os.system("cat" + listoffiles3 \
                + "> Transmittance.txt")
            for fle in Afiles + Rfiles + Tfiles:
                os.remove(fle)  
        # # T, R matricies
        # filename1 = "wl*_T_Lambda_MAT_sp.txt"
        # Afiles = glob.glob(filename1)
        # Afiles.sort()
        # filename2 = "wl*_R_Lambda_MAT_sp.txt"
        # Rfiles = glob.glob(filename2)
        # Rfiles.sort()
        # if len(Afiles)+len(Rfiles)>0:
        #     listoffiles1 = ""
        #     for fle in Afiles:
        #         listoffiles1 += " " + fle
        #     os.system("cat" + listoffiles1 \
        #         + "> T_Lambda_MAT_sp.txt")
        #     listoffiles2 = ""
        #     for fle in Rfiles:
        #         listoffiles2 += " " + fle
        #     os.system("cat" + listoffiles2 \
        #         + "> R_Lambda_MAT_sp.txt")
        #     for fle in Afiles + Rfiles:
        #         os.remove(fle)
    elif pol == 5:
        # Right hand circular polarisation
        filename1 = "wl*_A_Lambda_R.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Lambda_R.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "wl*_T_Lambda_R.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
        if len(Afiles)+len(Rfiles)+len(Tfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> Absorptance_R.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> Reflectance_R.txt")
            listoffiles3 = ""
            for fle in Tfiles:
                listoffiles3 += " " + fle
            os.system("cat" + listoffiles3 \
                + "> Transmittance_R.txt")
            for fle in Afiles + Rfiles + Tfiles:
                os.remove(fle)
        # Left hand circular polarisation
        filename1 = "wl*_A_Lambda_L.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Lambda_L.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "wl*_T_Lambda_L.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
        if len(Afiles)+len(Rfiles)+len(Tfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> Absorptance_L.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> Reflectance_L.txt")
            listoffiles3 = ""
            for fle in Tfiles:
                listoffiles3 += " " + fle
            os.system("cat" + listoffiles3 \
                + "> Transmittance_L.txt")
            for fle in Afiles + Rfiles + Tfiles:
                os.remove(fle)
        # Circular Dichroism
        filename1 = "wl*_A_Lambda_CD.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Lambda_CD.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "wl*_T_Lambda_CD.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
        if len(Afiles)+len(Rfiles)+len(Tfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> Absorptance_CD.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> Reflectance_CD.txt")
            listoffiles3 = ""
            for fle in Tfiles:
                listoffiles3 += " " + fle
            os.system("cat" + listoffiles3 \
                + "> Transmittance_CD.txt")
            for fle in Afiles + Rfiles + Tfiles:
                os.remove(fle)
        # T, R matricies
        filename1 = "wl*_T_Lambda_MAT_lr.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Lambda_MAT_lr.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        if len(Afiles)+len(Rfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> T_Lambda_MAT_lr.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> R_Lambda_MAT_lr.txt")
            for fle in Afiles + Rfiles:
                os.remove(fle)
    else:
        filename1 = "wl*_A_Lambda.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Lambda.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "wl*_T_Lambda.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
        if len(Afiles)+len(Rfiles)+len(Tfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> Absorptance.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> Reflectance.txt")
            listoffiles3 = ""
            for fle in Tfiles:
                listoffiles3 += " " + fle
            os.system("cat" + listoffiles3 \
                + "> Transmittance.txt")
            for fle in Afiles + Rfiles + Tfiles:
                os.remove(fle)   
        filename1 = "wl*_A_Layers.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "wl*_R_Layers.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "wl*_T_Layers.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
        if len(Afiles)+len(Rfiles)+len(Tfiles)>0:
            listoffiles1 = ""
            for fle in Afiles:
                listoffiles1 += " " + fle
            os.system("cat" + listoffiles1 \
                + "> Abs_Layers.txt")
            listoffiles2 = ""
            for fle in Rfiles:
                listoffiles2 += " " + fle
            os.system("cat" + listoffiles2 \
                + "> Ref_Layers.txt")
            listoffiles3 = ""
            for fle in Tfiles:
                listoffiles3 += " " + fle
            os.system("cat" + listoffiles3 \
                + "> Tran_Layers.txt")
            for fle in Afiles + Rfiles + Tfiles:
                os.remove(fle)



def c_c_omega(st):
    """         for nanostructured layers
        concatenate dispersion results in one file 
        and remove individual lambda results """

    format_st     = '%04d' % st
    filename1 = "st%s*_omega.txt" % format_st
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> omega_st%s.txt" % format_st)
        for fle in Afiles:
            os.remove(fle)
        
    filename1 = "st%s*_omega_pol.txt" % format_st
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> omega_pol_st%s.txt" % format_st)
        for fle in Afiles:
            os.remove(fle)
        
    filename1 = "st%s*_omega_Fz.txt" % format_st
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> omega_Fz_st%s.txt" % format_st)
        for fle in Afiles:
            os.remove(fle)
        
    filename1 = "st%s*_omega_Ft.txt" % format_st
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> omega_Ft_st%s.txt" % format_st)
        for fle in Afiles:
            os.remove(fle)



def c_c_beta(st):
    """             for thin film layers
        concatenate dispersion results in one file 
        and remove individual lambda results """

    format_st     = '%04d' % st
    filename1 = "st%s*_beta.txt" % format_st
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> beta_st%s.txt" % format_st)
        for fle in Afiles:
            os.remove(fle)


   
def c_c_detA():
    """ concatenate Fabry-Perot resonance determinant results in one file 
        and remove individual lambda results """

    filename1 = "detAe*.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> det_I-RPRP_E.txt")
        for fle in Afiles:
            os.remove(fle)

    filename2 = "detAo*.txt"
    Afiles = glob.glob(filename2)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles2 = ""
        for fle in Afiles:
            listoffiles2 += " " + fle
        os.system("cat" + listoffiles2 + "> det_I-RPRP_O.txt")
        for fle in Afiles:
            os.remove(fle)



def c_c_prop_modes():
    """ concatenate number of propagating modes (lossless) """

    filename1 = "SuperModes*.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    if len(Afiles)>0:
        listoffiles1 = ""
        for fle in Afiles:
            listoffiles1 += " " + fle
        os.system("cat" + listoffiles1 + "> SuperModes.txt")
        for fle in Afiles:
            os.remove(fle)
