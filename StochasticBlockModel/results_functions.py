#!/usr/bin/env python
# coding: utf-8

# In[1]:


def probabilities(theta, eta, pr, K, L, n_cond, n_moa, cond_test):
    '''Function that creates a dataframe with the probabilities computed for each condition and moa pair of 
    the test set. The dataframe will be as large as the number of pairs to be predicted.
    :param theta: list, parameter theta (condition membership vectors) 
    :param eta: list, parameter theta (moa membership vectors) 
    :param pr: list, parameter p matrix (probability matrix p)
    :param K: integer, number of condition groups
    :param L: integer, number of moa groups
    :param n_cond: integer, number of conditions in the whole dataset
    :param n_moa: integer, number of moas in the whole dataset
    :param cond_test: list of integers, conditions used as the test set
    :return df_probas: dataframe with the probabilities(1) for each condition-moa pair
    '''
    ### creating data structures for saving the predictions ###
    c_list = []
    m_list = []
    proba_list = []
    for c in cond_test: #fold 1 --> cond_test_index[0], so... fold-1
        for m in range(n_moa):
            proba_1 = 0
            for alpha in range(K):
                for beta in range(L):
                    proba_1 += theta[c][alpha]*pr[alpha][beta]*eta[m][beta]
            c_list.append(c)
            m_list.append(m)
            proba_list.append(proba_1)
    df_probas = pd.DataFrame({'cond': c_list, 'moa': m_list, 'proba1': proba_list})
    return(df_probas)

def predictions_and_scores(fold, file_cond_test, file_parameters, K, L, T, n_cond, n_moa, n_gen, runs):
    '''Function that computes the average probabilities given all the parameters of each run. Besides 
    computing the probabilities(1) it also calculates the scores: precision, recall and logistic loss.
    :param fold: integer, fold that have been used as the training set
    :param file_cond_test: txt file with integers, conditions used as the test set
    :param file_parameters: dat file, last parameters (theta, eta, phi, pr and qr)
    :param K: integer, number of condition groups
    :param L: integer, number of moa groups
    :param T: integer, number of gene groups
    :param n_cond: integer, number of conditions in the whole dataset
    :param n_moa: integer, number of moas in the whole dataset
    :param n_gen: integer, number of genes in the whole dataset
    :param runs: integer, number of runs
    :return (probas, log_loss_score)
    '''
    import pandas as pd
    #reading the list with the conditions for the test:
    fh=open(file_cond_test,'r')
    igot=fh.readlines()
    cond_test = []
    for line in igot:
        cond_test.append(int(line))
    fh.close()

    #reading the parameters + computing the probabilities + creating the dataframe of probabilities:
    import function_read_param as reading
    for i in range(1, runs+1):
        (theta,eta,phi,pr,qr) = reading.read_parameters(file_parameters+str(i)+'.dat',K,L,T,n_cond,n_moa,n_gen)
        if i==1: #create new df
            probas = probabilities(theta, eta, pr, K, L, n_cond, n_moa, cond_test)
            probas['proba1_r1'] = probas['proba1']
            probas = probas.drop(['proba1'], axis=1) 
        else: #update df
            probas['proba1_r'+str(i)] = probabilities(theta, eta, pr, K, L, n_cond, n_moa, cond_test)['proba1']
    cols = ['proba1_r'+str(i) for i in range(1,runs+1)] #columns to be used when computing the avg probability
    probas['avg_proba1'] = probas[cols].mean(axis=1) #adding a new column in the df with the avg probability
    
    #if a condition is a control sample, their predictions must be 0.0
    features = pd.read_csv("/Users/belbordes/Desktop/MoA Project/lish-moa/train_features.csv", index_col=0)
    ctrl_conds = features[features['cp_type'] == 'ctl_vehicle'].index.tolist()
    ### dictionaries (condition id <----> condition integer)
    fh=open('/Users/belbordes/Desktop/MoA Project/Code_revision_Toñi/Toñi/IDData_agodoy/conditionsID.csv','r')
    igot=fh.readlines()
    ids = []
    ints = []
    for line in igot:
        about = line.strip().split(',')
        ids.append(about[0])
        ints.append(int(about[1]))
    fh.close()
    dictionary = dict(zip(ids, ints))
    dictionary_r = dict(zip(ints, ids))
    #get the integers of the control conditions:
    ctrl_conds_int = []
    for c in ctrl_conds:
        ctrl_conds_int.append(dictionary[c])
    #probabilities associated with these conditions will be written by us with 0.0
    probas['treated'] = [1 for i in range(probas.shape[0])] #adding a column
    probas.loc[(probas.cond).isin(ctrl_conds_int),'treated']=0 #ctrl will be 0
    probas[probas.treated==0]
    probas['avg_proba1_ctrl0'] = probas['avg_proba1']*probas['treated'] #ctrl samples will be x0 --> avg proba = 0.0
    
    #reading the actual links:
    df_test = pd.read_csv('/Users/belbordes/Desktop/MoA Project/Code_revision_Toñi/Gemma_1abril/5foldCV_data/F'+
                          str(fold)+'_test_dataframe.csv', index_col=0)
    
    #scores: logistic loss
    #we must create a dataframe with the probabilities following the same structure as the test dataframe
    print('\nCalculating logistic loss...')
    df_probas = pd.DataFrame(columns=[i for i in range(206)], index=cond_test)
    for i in range(probas.shape[0]):
        cond = probas.cond[i]
        moa = probas.moa[i]
        pr = probas.avg_proba1_ctrl0[i]
        df_probas.loc[cond, moa] = pr
    from sklearn.metrics import log_loss #scikit function
    log_loss_list = []
    for i in range(206):
        log_loss_list.append(log_loss(df_test.iloc[:, i], df_probas.iloc[:, i], labels = [0, 1]))
    import numpy as np
    log_loss_arr = np.array(log_loss_list)
    log_loss_score = np.mean(log_loss_arr)
    print('log loss:', log_loss_score)
    
    return(probas, log_loss_score)

def predictions_and_scores_withCells(fold, file_cond_test, file_parameters, K, L, T, F, n_cond, 
                                     n_moa, n_gen, n_cell, runs):
    '''Function that computes the average probabilities given all the parameters of each run. Besides 
    computing the probabilities(1) it also calculates the scores: precision, recall and logistic loss.
    :param fold: integer, fold that have been used as the training set
    :param file_cond_test: txt file with integers, conditions used as the test set
    :param file_parameters: dat file, last parameters (theta, eta, phi, pr and qr)
    :param K: integer, number of condition groups
    :param L: integer, number of moa groups
    :param T: integer, number of gene groups
    :param F: integer, number of cell groups
    :param n_cond: integer, number of conditions in the whole dataset
    :param n_moa: integer, number of moas in the whole dataset
    :param n_gen: integer, number of genes in the whole dataset
    :param n_cell: integer, number of cells in the whole dataset
    :param runs: integer, number of runs
    :return (probas, log_loss_score)
    '''
    import pandas as pd
    #reading the list with the conditions for the test:
    fh=open(file_cond_test,'r')
    igot=fh.readlines()
    cond_test = []
    for line in igot:
        cond_test.append(int(line))
    fh.close()

    #reading the parameters + computing the probabilities + creating the dataframe of probabilities:
    import function_read_param as reading
    for i in range(1, runs+1):
        (theta,eta,phi,xi,pr,qr,sr) = reading.read_parameters_withCells(file_parameters+str(i)+'.dat',K,L,T,F,
                                                                  n_cond,n_moa,n_gen)
        if i==1: #create new df
            probas = probabilities(theta, eta, pr, K, L, n_cond, n_moa, cond_test)
            probas['proba1_r1'] = probas['proba1']
            probas = probas.drop(['proba1'], axis=1) 
        else: #update df
            probas['proba1_r'+str(i)] = probabilities(theta, eta, pr, K, L, n_cond, n_moa, cond_test)['proba1']
    cols = ['proba1_r'+str(i) for i in range(1,runs+1)] #columns to be used when computing the avg probability
    probas['avg_proba1'] = probas[cols].mean(axis=1) #adding a new column in the df with the avg probability
    
    #if a condition is a control sample, their predictions must be 0.0
    features = pd.read_csv("/Users/belbordes/Desktop/MoA Project/lish-moa/train_features.csv", index_col=0)
    ctrl_conds = features[features['cp_type'] == 'ctl_vehicle'].index.tolist()
    ### dictionaries (condition id <----> condition integer)
    fh=open('/Users/belbordes/Desktop/MoA Project/Code_revision_Toñi/Toñi/IDData_agodoy/conditionsID.csv','r')
    igot=fh.readlines()
    ids = []
    ints = []
    for line in igot:
        about = line.strip().split(',')
        ids.append(about[0])
        ints.append(int(about[1]))
    fh.close()
    dictionary = dict(zip(ids, ints))
    dictionary_r = dict(zip(ints, ids))
    #get the integers of the control conditions:
    ctrl_conds_int = []
    for c in ctrl_conds:
        ctrl_conds_int.append(dictionary[c])
    #probabilities associated with these conditions will be written by us with 0.0
    probas['treated'] = [1 for i in range(probas.shape[0])] #adding a column
    probas.loc[(probas.cond).isin(ctrl_conds_int),'treated']=0 #ctrl will be 0
    probas[probas.treated==0]
    probas['avg_proba1_ctrl0'] = probas['avg_proba1']*probas['treated'] #ctrl samples will be x0 --> avg proba = 0.0
    
    #reading the actual links:
    df_test = pd.read_csv('/Users/belbordes/Desktop/MoA Project/Code_revision_Toñi/Gemma_1abril/5foldCV_data/F'+
                          str(fold)+'_test_dataframe.csv', index_col=0)
    
    #scores: logistic loss
    #we must create a dataframe with the probabilities following the same structure as the test dataframe
    print('\nCalculating logistic loss...')
    df_probas = pd.DataFrame(columns=[i for i in range(206)], index=cond_test)
    for i in range(probas.shape[0]):
        cond = probas.cond[i]
        moa = probas.moa[i]
        pr = probas.avg_proba1_ctrl0[i]
        df_probas.loc[cond, moa] = pr
    from sklearn.metrics import log_loss #scikit function
    log_loss_list = []
    for i in range(206):
        log_loss_list.append(log_loss(df_test.iloc[:, i], df_probas.iloc[:, i], labels = [0, 1]))
    import numpy as np
    log_loss_arr = np.array(log_loss_list)
    log_loss_score = np.mean(log_loss_arr)
    print('log loss:', log_loss_score)
    
    return(probas, log_loss_score)

