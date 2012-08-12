# concatination of results and clean up

import glob
import os


def c_c_tra(pol = 'TM'):
    """ concatinate transmission, reflection and absorption results 
    in one file each and remove individual lambda results """

    if pol == 5:
        # Right hand circular polarisation
        filename1 = "p*_A_Lambda_R.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "p*_R_Lambda_R.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "p*_T_Lambda_R.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
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
        filename1 = "p*_A_Lambda_L.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "p*_R_Lambda_L.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "p*_T_Lambda_L.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
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
        filename1 = "p*_A_Lambda_CD.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "p*_R_Lambda_CD.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "p*_T_Lambda_CD.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
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
    else:
        filename1 = "p*_A_Lambda.txt"
        Afiles = glob.glob(filename1)
        Afiles.sort()
        filename2 = "p*_R_Lambda.txt"
        Rfiles = glob.glob(filename2)
        Rfiles.sort()
        filename3 = "p*_T_Lambda.txt"
        Tfiles = glob.glob(filename3)
        Tfiles.sort()
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



def c_c_omega():
    """ concatinate dispersion results in one file 
        and remove individual lambda results """

    filename1 = "p*_omega.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    listoffiles1 = ""
    for fle in Afiles:
        listoffiles1 += " " + fle
    os.system("cat" + listoffiles1 + "> omega.txt")
    for fle in Afiles:
        os.remove(fle)
        
    filename1 = "p*_omega_pol.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    listoffiles1 = ""
    for fle in Afiles:
        listoffiles1 += " " + fle
    os.system("cat" + listoffiles1 + "> omega_pol.txt")
    for fle in Afiles:
        os.remove(fle)
        
    filename1 = "p*_omega_Fz.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    listoffiles1 = ""
    for fle in Afiles:
        listoffiles1 += " " + fle
    os.system("cat" + listoffiles1 + "> omega_Fz.txt")
    for fle in Afiles:
        os.remove(fle)
        
    filename1 = "p*_omega_Ft.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    listoffiles1 = ""
    for fle in Afiles:
        listoffiles1 += " " + fle
    os.system("cat" + listoffiles1 + "> omega_Ft.txt")
    for fle in Afiles:
        os.remove(fle)


   
def c_c_detA():
    """ concatinate Fabry-Perot resonance determinant results in one file 
        and remove individual lambda results """

    filename1 = "detAe*.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    listoffiles1 = ""
    for fle in Afiles:
        listoffiles1 += " " + fle
    os.system("cat" + listoffiles1 + "> det_I-RPRP_E.txt")
    for fle in Afiles:
        os.remove(fle)

    filename2 = "detAo*.txt"
    Afiles = glob.glob(filename2)
    Afiles.sort()
    listoffiles2 = ""
    for fle in Afiles:
        listoffiles2 += " " + fle
    os.system("cat" + listoffiles2 + "> det_I-RPRP_O.txt")
    for fle in Afiles:
        os.remove(fle)



def c_c_prop_modes():
    """ concatinate number of propagating modes (lossless) """

    filename1 = "SuperModes*.txt"
    Afiles = glob.glob(filename1)
    Afiles.sort()
    listoffiles1 = ""
    for fle in Afiles:
        listoffiles1 += " " + fle
    os.system("cat" + listoffiles1 + "> SuperModes.txt")
    for fle in Afiles:
        os.remove(fle)
