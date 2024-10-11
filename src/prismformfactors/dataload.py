def load_data(filepath:str):
    """Loads all data from a .dat file as three lists : q, I(q), and dI(q)"""
    cursor = open(filepath,"r")
    q = []
    Iq = []
    Errorq = []
    for ligne in cursor:
        if not(ligne[0].isdigit()) or ("nan" in ligne.split()) or float(ligne.split()[1])==0. :
            pass
        else :
            data = ligne.split()
            q.append(float(data[0]))
            Iq.append(float(data[1]))
            Errorq.append(float(data[2]))
    cursor.close()
    return q, Iq, Errorq