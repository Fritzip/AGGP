import numpy as np

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

#insertion in the compressed matrix
def ins_in_compr(m,ins_index,ins_bit):
    mc = compress(m)
    return uncompress(mc[:ins_index]+[ins_bit]+mc[ins_index:-1])

#deletion in the compressed matrix
def del_in_compr(m,del_index):
    mc = compress(m)
    return uncompress(mc[:del_index]+mc[del_index+1:]+[0])

def substitution(m,i,j):
    #i<dim-1 and j>=i
    m[i,j]=(m[i,j]+1)%2

#insertion in the complete adjacency matrix
def insertion(m,i,j):        
    #i<dim-1 and j>=i
    for x in range(dim-1,i,-1):
       pass 


m = np.random.randint(2,size=(4,4))
print "random matrix\n",m
m = symmetrize(m)
print "symmetrized one\n",m
mc = compress(m)
print "compressed one\n",mc
m = uncompress(mc)
print "uncompressed one\n",m
m = ins_in_compr(m,1,1)
print m
m = del_in_compr(m,1)
print m
