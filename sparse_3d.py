import numpy as np
from scipy.sparse import coo_matrix

# domain size but excluding boundary
J = I = K = 4

def n2ijk(n):
    # n=IJk+Ij+i
    k = int(( n - ( n%(I*J) ) ) / (I*J))
    i = int((n - k*I*J)%I)
    j = int((n - I*J*k - i)/I)
    return i, j, k

# n = 0,...,I*J*K-1
# i = 0,...,I-1
# j = 0,...,J-1
# k = 0,...,K-1


sparse = np.zeros((K*J*I, K*J*I), dtype=np.short)
rows   = np.array([],dtype=np.short)
cols   = np.array([],dtype=np.short)
data   = np.array([],dtype=np.short)


for m in range(I*J*K):
  # fix mask
  mask_i, mask_j, mask_k = n2ijk(m)
  for n in range(I*J*K):
  # swipe entire domain but excluding boundary
      i, j, k = n2ijk(n)
    
      # volume
      if 1<=i<=I-2 and 1<=j<=J-2 and 1<=k<=K-2:
         if   abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      # 4 surfaces
      elif  i==0   and 1<=j<=J-2 and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif  i==I-1 and 1<=j<=J-2 and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif  1<=i<=I-2 and    j==0   and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and    j==J-1 and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and 1<=j<=J-2 and    k==0  :
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and 1<=j<=J-2 and    k==K-1:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0



      # 12 lines
      elif    i==0   and    j==0   and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==0   and 1<=j<=J-2 and    k==0  :
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and    j==0   and    k==0  :
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==I-1 and    j==J-1 and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==I-1 and 1<=j<=J-2 and    k==K-1:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and    j==J-1 and    k==K-1:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==0   and    j==J-1 and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==0   and 1<=j<=J-2 and    k==K-1:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and    j==0   and    k==K-1:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==I-1 and    j==0   and 1<=k<=K-2:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif    i==I-1 and 1<=j<=J-2 and    k==0  :
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0

      elif 1<=i<=I-2 and    j==J-1 and    k==0  :
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 1
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0


      # 8 points
      else:
         if   abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
            sparse[m][n] = 2
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
            sparse[m][n] = -6
            rows = np.append(rows, np.array([m]))
            cols = np.append(cols, np.array([n]))
            data = np.append(data, np.array([sparse[m][n]]))
         else:
            sparse[m][n] = 0




sparse = coo_matrix((data, (rows, cols)), shape=(I*J*K, I*J*K)).toarray()



# plot 2D sparse matrx of 2D domain
for m in range(I*J*K):
   for n in range(I*J*K):
        print("%+d" % sparse[m][n], end ='')
   print("")
