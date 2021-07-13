# calling TAD using deTOKI
# migrate to notebook compatible 
import numpy as np
import pandas as pd
from sklearn import decomposition
import math
import multiprocessing

## Run NMF in several times with random initialisation and output consensus matrix
def corate(A,n,time):
    # Input:
    #   A : matrix, window-sized submatrix
    #   n : number of components, optimized by 
    #       {average_length_of_TAD, length_of_window}
    #   time : times of NMF running
    # Output:
    #   consensus NMF matrix
    S=np.zeros([np.shape(A)[0],np.shape(A)[1]])
    for i in range(time):
        K=np.zeros([np.shape(A)[0],np.shape(A)[1]])
        estimator=decomposition.NMF(n_components=n, init='random',random_state=i, max_iter=10_000)
        estimator.fit(A)
        B=estimator.transform(A)
        index=B.argmax(axis=1)
        for j in range(n):
            for s in np.where(index==j)[0]:
                for t in np.where(index==j)[0]:
                    K[s,t]=1
        S=S+K
    return S.astype(np.float64)/time
## Define silhouette coefficient of consensus matrix and detected TAD boundaries
def silhou(R, pos):
    n=np.shape(R)[0]
    silhou=0
    for i in range(len(pos)-1):
        for j in range(pos[i],pos[i+1]):
            a=np.sum((1-R)[j,pos[i]:pos[i+1]])
            b=np.sum((1-R)[i,:])-a
            silhou+=(-a/(pos[i+1]-pos[i])+b/(n+pos[i]-pos[i+1]))/max((a/(pos[i+1]-pos[i]),b/(n+pos[i]-pos[i+1])))
    return silhou/n

## Calculate clustering rate of each bin
def IS(R, length):
    # Calculate clustering rate(number of contacts 
    # locatingin  inter_domain triangle field of the bin)

    # Input:
    #   R: contact matrix, typically consensus_NMF_matrix
    #   length : length of bin's flanking regions considered
    bias=np.zeros([np.shape(R)[0]])
    for i in range(1,np.shape(R)[0]-1):
        bias[i]=np.mean(
            R[max(0,i-length):(i+1), i:min(i+length+1,np.shape(R)[0])]
            )
    return bias

## Find the best n_components by comparing silhouette coefficient
def bestco(F, reso, length, delta, min_val, max_val):
    # argument picking
    # find best n(number of components)
    # judged by silhouette coefficient

    # n ~ number of TADs in 1 window

    # rely on:
    #   corate, silhou, zero
    # Input:
    #   F, reso, min_val, max_val: for n determining
    #   length, delta: for boundary calling
    # Output:
    #   final consensus NMF matrix, and its n
    x=-1
    R0=0
    n1=0
    for n in range(max(int(np.shape(F)[0]*reso/max_val),1),int(np.shape(F)[0]*reso/min_val)+1):
        # get NMF_consensus_matrix
        R1=corate(F,n,10)
        # calc tad under above n/NMF_consensus_matrix
        x1=silhou(R1, zero(R1, n-1, length, delta))
        if x1>=x:
            R0=R1
            x=x1
            n1=n
    # final R0, best consensus NMF matrix
    # fianl n1, best n_component used in NMF
    return R0,n1

## Find bins with local minimal clustering rate and global comparative low clustering rate （these bins are detected TAD boundaries)
def zero(R, t, length, delta):
    # pick at most t TAD boundary from contact matrix
    # Input:
    #   R : small piece of matrix
    #   t : number of output boundaries

    # find local valley of IS result

    #bias: number of contacts in "anti-diag" field
    bias=IS(R, length)
    #delta_list: right_bias - left_bias
    delta_list=[]
    for i in range(delta,len(bias)-delta):
        delta_list.append(-sum(bias[(i-delta):i])+sum(bias[i:(i+delta)]))
    # zero: position of local valleys
    zero=[0]
    for i in range(len(delta_list)-1):
        if delta_list[i]<0 and delta_list[i+1]>=0:
            #zero.append(i+5-np.argmin([bias[i+x] for x in range(5,0,-1)]))
            zero.append(i+delta) # why "+delta"?
    zero.append(len(bias))
    zero=sorted(np.unique(zero))

    # pick boundaries near strong TAD

    # enrich: subset of zero, picked by index
    enrich=[0]
    # strength: max summit value upstream and downstream the valley
    #   aka inner tad strength of tad's alongside of the boundary(zero pos)
    strength=[]
    # dicard left and right end
    for j in range(1,len(zero)-1):
        strength.append(
            max(
                # upstream summit
                max( bias[zero[j-1]:(zero[j]+1)] ) - bias[zero[j]],
                # downstream summit
                max( bias[(zero[j]-1):zero[j+1]] ) - bias[zero[j]]
                )
            )
    strength=np.array(strength)
    # select globally from zeros, using strength
    # top t, and value>0.3
    index = sorted(list(set(np.argsort(-strength)[:t])|set(np.where(strength>0.3)[0])))
    # store result in enrich, adding left and right end
    for k in index:
        enrich.append(zero[k+1])
    enrich.append(len(bias))
    return np.array(enrich)
## Split huge contact matrix to windows and detected TAD boundaries in each window
def part_zero(F:np.ndarray, core:int, window:float, reso:int, delta, length, min_val, max_val):
    # depend on : bestco, zero
    
    #Input:
    #   F: full contact_matrix
    #   core: number of threads
    #   window, reso, length, delta -> bestco
    #   length, delta -> zero
    
    pos=[]
    # length of chromosome/input_matrix, in num of reso_bins
    n=np.shape(F)[0]

    # deTOKI's core function

    global task
    def task(i):
        # i : work id
        # (implicit)F :  full mat
        # (implicit)window : window size (uint: reso_bins)
        # (implicit)n : length of full chromosome 

        # P : a window-size-tile of F, different near chromosome ends
        P=F[max(0,window//2*i-window//4):min(n,window//2*i+window-window//4),
            max(0,window//2*i-window//4):min(n,window//2*i+window-window//4)]
        if np.sum(P)<100:
            # poor signal
            # windows with less than 100 contacts yield 0 tad boundaries
            return []
        # find best NMF_consensus_matrix
        R, t = bestco(P, reso, length, delta, min_val, max_val)
        if t==0:
            # all n candidates are not good
            return []
        # call TAD on consensus matrix
        p=zero(R, t-1, length, delta)
        pos=[]
        # transform TAD boundary from local(in window) index to global index
        for j in p:
            if j in range(window//2*i-max(0,window//2*i-window//4),window//2*(i+1)-max(0,window//2*i-window//4)):
                pos.append(j+max(0,window//2*i-window//4))
        return pos

    # multi-thread calc
    
    pool = multiprocessing.Pool(processes=core)
    res_list=[]
    for i in range(math.ceil(2*n/window)):
        # iterate by window's mid-points
        # window size is 8mb,
        # each window overlaps neighbor windows by 4mb,
        # yield 2*n/window mid-points
        result = pool.apply_async(task, args=(i,))
        res_list.append(result)
    pool.close()
    pool.join()

    # async get and join every window's tad boundary

    for i in range(math.ceil(2*n/window)):
        result=res_list[i]
        pos=np.append(pos,result.get())
    return pos.astype('int32')
## entry function
def tad(F, core:int=1, reso:int=40, min_val:int=600, max_val:int=1000, window_size_raw:int=8000, delta_raw:int=100):
#Input:
#   F: input contact matrix : np.ndarray
#   reso: resolution of the matrix (kbp)
#   min, max: tad mean size (kbp)
#   split: window size (kbp) 
#   core: treads used

    # length of flanking regions considered
    # by default 400kbp
    length=400//reso
    # small intervals for calculating curve's summints and valleys
    delta=int(math.ceil(delta_raw/reso))
    # tiling windows' size(number of res_bins), for parallel calculation 
    window_size = window_size_raw//reso

    return part_zero(F, core, window, reso, delta, length, min_val, max_val)