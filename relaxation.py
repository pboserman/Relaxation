



def alpha(mu, mu_p):
    return 3*(mu*mu_p)**2-mu**2-mu_p**2+3

def beta(mu, mu_p):
    return 5*(mu*mu_p+(mu*mu_p)**3)-3*(mu**3*mu_p-mu*mu_p**3)

def epsilon(mu, mu_p):
    return -3*alpha(mu, mu_p)+8+2*beta(mu, mu_p)-8*mu*mu_p

