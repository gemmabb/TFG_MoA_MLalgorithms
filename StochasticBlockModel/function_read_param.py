#!/usr/bin/env python
# coding: utf-8

# In[ ]:


def read_parameters(file, K, L, T, n_cond, n_moa, n_gen):
    """Function that read the .dat file containing the optimized parameters (after having run the code with pypy).
       :param file: string with the name of the .dat file
       :param K: number of groups of conditions
       :param L: number of groups of moas
       :param T: number of groups of genes
       :param n_cond: total number of conditions
       :param n_moa: total number of moas
       :param n_gen: total number of genes
       :return tupple with the parameters, following this order --> theta(0), eta(1), phi(2), pr(3), qr(4)
    """
    with open(file) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip() #remove all the \n
    lines.pop(0) #remove this first line (iteration + likelihood)
    ### Theta ###
    theta = []
    for i in range(0, n_cond):
        theta.append(lines[i].split()) 
    #convert strings into floats:
    for c in range(n_cond):
        for alpha in range(K):
            theta[c][alpha] = float(theta[c][alpha])
    ### Eta ###
    eta = []
    for i in range(n_cond, n_cond+n_moa):
        eta.append(lines[i].split())
    #convert strings into floats:
    for m in range(n_moa):
        for beta in range(L):
            eta[m][beta] = float(eta[m][beta])
    ### Phi ###
    phi = []
    for i in range(n_cond+n_moa, n_cond+n_moa+n_gen):
        phi.append(lines[i].split())
    #convert strings into floats:
    for g in range(n_gen):
        for gamma in range(T):
            phi[g][gamma] = float(phi[g][gamma])
    ### P mat ###
    pr = []
    for i in range(n_cond+n_moa+n_gen, n_cond+n_moa+n_gen+K):
        pr.append(lines[i].split())

    #convert strings into floats:
    for alpha in range(K):
        for beta in range(L):
            pr[alpha][beta] = float(pr[alpha][beta])
    ### Q mat ###
    qr_lines = lines[n_cond+n_gen+n_moa+K:]
    qr_lines = [x for x in qr_lines if x != ''] #remove ''
    qr = []
    i = 0
    for e in range(3): #links
        sub_qr = []
        for l in qr_lines[i:i+K]:
            sub_qr.append(l.split())
        qr.append(sub_qr)
        i+=K
    #convert strings into floats:
    for e in range(3):
        for alpha in range(K):
            for gamma in range(T):
                qr[e][alpha][gamma] = float(qr[e][alpha][gamma])
                
    return(theta, eta, phi, pr, qr)

def read_parameters_withCells(file, K, L, T, F, n_cond, n_moa, n_gen, n_cell):
    """Function that read the .dat file containing the optimized parameters (after having run the code with pypy).
       :param file: string, name of the .dat file
       :param K: integer, number of groups of conditions
       :param L: integer, number of groups of moas
       :param T: integer, number of groups of genes
       :param F: integer, number of groups of cells
       :param n_cond: integer, total number of conditions
       :param n_moa: integer, total number of moas
       :param n_gen: integer, total number of genes
       :param n_cell: integer, total number of cells
       :return tupple with the parameters, following this order --> theta(0), eta(1), phi(2), pr(3), qr(4)
    """
    with open(file) as f:
        lines = f.readlines()
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip() #remove all the \n
    lines.pop(0) #remove this first line (iteration + likelihood)
    ### Theta ###
    theta = []
    for i in range(0, n_cond):
        theta.append(lines[i].split()) 
    #convert strings into floats:
    for c in range(n_cond):
        for alpha in range(K):
            theta[c][alpha] = float(theta[c][alpha])
    ### Eta ###
    eta = []
    for i in range(n_cond, n_cond+n_moa):
        eta.append(lines[i].split())
    #convert strings into floats:
    for m in range(n_moa):
        for beta in range(L):
            eta[m][beta] = float(eta[m][beta])
    ### Phi ###
    phi = []
    for i in range(n_cond+n_moa, n_cond+n_moa+n_gen):
        phi.append(lines[i].split())
    #convert strings into floats:
    for g in range(n_gen):
        for gamma in range(T):
            phi[g][gamma] = float(phi[g][gamma])
    ### Xi ###
    xi = []
    for i in range(n_cond+n_moa+n_gen, n_cond+n_moa+n_gen+n_cell):
        xi.append(lines[i].split())
    #convert strings into floats:
    for c in range(n_cell):
        for delta in range(F):
            xi[c][delta] = float(xi[c][delta])
    ### P mat ###
    pr = []
    for i in range(n_cond+n_moa+n_gen+n_cell, n_cond+n_moa+n_gen+n_cell+K):
        pr.append(lines[i].split())
    #convert strings into floats:
    for alpha in range(K):
        for beta in range(L):
            pr[alpha][beta] = float(pr[alpha][beta])
    ### Q mat ###
    qr_lines = lines[n_cond+n_moa+n_gen+n_cell+K:n_cond+n_moa+n_gen+n_cell+K+3*(K+1)] #3-> #possible links cond-gene (-1,0,1)
    qr_lines = [x for x in qr_lines if x != ''] #remove ''
    qr = []
    i = 0
    for e in range(3): #links
        sub_qr = []
        for l in qr_lines[i:i+K]:
            sub_qr.append(l.split())
        qr.append(sub_qr)
        i+=K
    #convert strings into floats:
    for e in range(3):
        for alpha in range(K):
            for gamma in range(T):
                qr[e][alpha][gamma] = float(qr[e][alpha][gamma])
    ### S mat ###
    sr = []
    for i in range(n_cond+n_moa+n_gen+n_cell+K+3*(K+1), n_cond+n_moa+n_gen+n_cell+K+3*(K+1)+K):
        sr.append(lines[i].split())
    #convert strings into floats:
    for alpha in range(K):
        for delta in range(F):
            sr[alpha][delta] = float(sr[alpha][delta])
                
    return(theta, eta, phi, xi, pr, qr, sr)

