

def getMuGrid(mu_size):
    muVals, weightVals = legpoly.legendre.leggauss(mu_size)
    return muVals

def alpha(mu, mu_p):
    return 3*(mu*mu_p)**2-mu**2-mu_p**2+3

def beta(mu, mu_p):
    return 5*(mu*mu_p+(mu*mu_p)**3)-3*(mu**3*mu_p-mu*mu_p**3)

def epsilon(mu, mu_p):
    return -3*alpha(mu, mu_p)+8+2*beta(mu, mu_p)-8*mu*mu_p

def n_e(z):
    return (1-z**2)**3

def Temp(z):
    return (1-z**2)

def eta(z, x):
    # T = Temp(z)
    # res = n_e(z)**2/np.sqrt(T)*np.exp(-np.power(10, x)/T)
    res = (1-z**2)**(5.5)*np.exp(-np.power(10, x)/(1-z**2))
    return np.nan_to_num(res)

def weight(mu):
    muVals, weightVals = legpoly.legendre.leggauss(muSize)
    return weightVals

