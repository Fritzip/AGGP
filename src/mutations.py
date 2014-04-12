import numpy as np
#import random as rand

#====================== Internal functions =============================

#Compress the matrix, based on the fact it's symetrical, diagonal is null and matrix is binary 
def compress(m):
    dim=m.shape[0]
    mLin=[]
    for i in range(dim-1):
        mLin+=list(m[i,i+1:])
    return mLin

#Recognize the original dimension of a compress squared matrix
def recognize_dim_from_compr(mLin):
    dim=1
    l=len(mLin)
    while dim<=l:
        l-=dim
        dim+=1
    return dim

#build a list of indices of the beginnings of each raws from the original matrix in the compressed one
def raws_begin_index(dim):
    return [0]+[np.sum(np.arange(i,dim)) for i in range(dim-1,1,-1)]

def uncompress(mLin):
    m=[]
    dim=recognize_dim_from_compr(mLin)
    rbi=raws_begin_index(dim)
    for i in range(dim-1):
        m+=[0]*(i+1)+mLin[rbi[i]:rbi[i]+dim-i-1]
    m+=[0]*dim
    m=np.array(m).reshape((dim,dim))
    return symmetrize(m)


#Symmetrize, ie copy superior triangle in inferior one.
def symmetrize(m):
    tri=np.triu(m,k=1)
    m = tri+np.transpose(tri)
    return(m)

#======================= Usable functions ========================

#insertion in the compressed matrix
def ins_in_compr(m,ins_index,ins_bit):
    #index < sum(range(dim))
    mc = compress(m)
    return uncompress(mc[:ins_index]+[ins_bit]+mc[ins_index:-1])

#deletion in the compressed matrix
def del_in_compr(m,del_index):
    #index < sum(range(dim))
    mc = compress(m)
    return uncompress(mc[:del_index]+mc[del_index+1:]+[0])

def substitution(m,i,j):
    #i<dim-1 and j>=i
    m[i,j]=(m[i,j]+1)%2
    return symmetrize(m)

#insertion in the complete adjacency matrix
def insertion(m,i,j,ins_bit):        
    #i<dim-1 and dim-1>=j>=i
    for x in range(m.shape[0]-1,j,-1):
        m[i,x]=m[i,x-1]
    for y in range(0,i):
        m[y,j]=m[y+1,j]
    m[i,j]=ins_bit
    return symmetrize(m)

#insertion in the complete adjacency matrix
def deletion(m,i,j):        
    #i<dim-1 and dim-1>=j>=i
    for x in range(j,m.shape[0]-1):
        m[i,x] = m[i,x+1]
    m[i,m.shape[0]-1] = 0
    for y in range(i,0,-1):
        m[y,j] = m[y-1,j]
    m[0,j] = 0
    
    return symmetrize(m)

#============================== Tests ================

if __name__=='__main__':
    m = np.random.randint(2,size=(6,6))
    print "random matrix\n",m
    m = symmetrize(m)
    print "symmetrized one\n",m
    mc = compress(m)
    print "compressed one\n",mc
    m = uncompress(mc)
    print "uncompressed one\n",m
    m = ins_in_compr(m,1,1)
    print "ins in comp\n",m
    m = del_in_compr(m,1)
    print "del in compr\n",m
    m = substitution(m,i=1,j=2)
    print "substitution\n",m
    m = insertion(m,i=2,j=3,ins_bit=1)
    print "insertion\n",m
    m = deletion(m,i=2,j=3)
    print "deletion\n",m
