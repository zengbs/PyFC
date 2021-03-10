import numpy as np

# domain size but excluding boundary
J = I = K = 5

def n2ijk(n):
    # n=IJk+Ij+i
    k = int(( n - ( n%(I*J) ) ) / (I*J))
    i = int((n - k*I*J)%I)
    j = int((n - I*J*k - i)/I)
    return i, j, k

# n = 0,...,I*J-1
# i = 0,...,I-1
# j = 0,...,J-1


sparse = np.zeros((K*J*I, K*J*I), dtype=np.short)

for m in range(I*J*K):
  # fix mask
  mask_i, mask_j, mask_k = n2ijk(m)
  for n in range(I*J*K):
  # swipe entire domain but excluding boundary
      i, j, k = n2ijk(n)
      if   abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
         sparse[m][n] = -6
      elif abs(i-mask_i)==1 and abs(j-mask_j)==0 and abs(k-mask_k)==0:
         sparse[m][n] = 1
      elif abs(i-mask_i)==0 and abs(j-mask_j)==1 and abs(k-mask_k)==0:
         sparse[m][n] = 1
      elif abs(i-mask_i)==0 and abs(j-mask_j)==0 and abs(k-mask_k)==1:
         sparse[m][n] = 1
      else:
         sparse[m][n] = 0



# plot 2D sparse matrx of 2D domain
for m in range(I*J*K):
   for n in range(I*J*K):
        print("%+d" % sparse[m][n], end ='')
   print("")
