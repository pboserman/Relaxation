

def G(z, mu):
    return tau*n_e(z)/mu

def A(z, x, mu, mu_p):
    if len(z.shape) != 4:
        z = z.reshape(-1, 1, 1, 1)
    if len(x.shape) != 4:
        x = x.reshape(1, -1, 1, 1)
    if len(mu.shape) != 4:
        mu = mu.reshape(1, 1, -1, 1)
    if len(mu_p.shape) != 4:
        mu_p = mu_p.reshape(1, 1, 1, -1)

    # print(x.shape, z.shape, mu.shape)
    return tau*n_e(z)/mu*(3./16.)*(alpha(mu, mu_p)*weight(mu_p))



def F(z, x, mu):
    return eta(z, x)/mu

zGrid_3 = zGrid.reshape(-1, 1, 1)
xGrid_3 = xGrid.reshape(1, -1, 1)
muGrid_3 = muGrid.reshape(1, 1, -1)

zGrid_4 = zGrid.reshape(-1, 1, 1, 1)
xGrid_4 = xGrid.reshape(1, -1, 1, 1)
muGrid_4 = muGrid.reshape(1, 1, -1, 1)
muGridp_4 = muGrid.reshape(1, 1, 1, -1)

def E_k(z, x, mu):
    I_i = I[1:-1]
    I_ip1 = I[2:]
    I_im1 = I[:-2]
    
    E_1 = I_ip1 - I_i + delZ*G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*(I_ip1+I_i) - delZ*F(z+0.5*delZ, x, mu) - delZ*np.sum(A(z, x, mu, muGridp_4)/2*np.expand_dims(I_i+I_ip1, axis = -2), axis = -1)
    E_2 = I_i - I_im1 + delZ*G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*(I_im1+I_i) - delZ*F(z-0.5*delZ, x, mu) - delZ*np.sum(A(z, x, mu, muGridp_4)/2*np.expand_dims(I_i+I_im1, axis = -2), axis = -1)
    
    return E_1*(muGrid < 0) + E_2*(muGrid > 0)

def S_k(z, x, mu, mu_p, which):
    if which == 'p':
        return (G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ - 1)*(mu==mu_p) - A(z+0.5*delZ, x, mu, mu_p)/2*delZ
    if which == 'n':
        return (G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ + 1)*(mu==mu_p) - A(z-0.5*delZ, x, mu, mu_p)/2*delZ
    else:
        print("Option not available")
        exit(1)

def S_km1(z, x, mu, mu_p):
    return (G(z-0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ - 1)*(mu==mu_p) - A(z-0.5*delZ, x, mu, mu_p)/2*delZ

def S_kp1(z, x, mu, mu_p):
    return (G(z+0.5*delZ, mu).repeat(xSize, axis = 1)/2*delZ + 1)*(mu==mu_p) - A(z+0.5*delZ, x, mu, mu_p)/2*delZ

