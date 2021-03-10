import numpy as np

# domain size but excluding boundary
J = I = 5

def ij2n(i,j):
    return i+j*I

def n2ij(n):
    i = n%I
    j = int((n-i)/J)
    return i, j



sparse = np.zeros((J*I, J*I), dtype=np.short)

for m in range(I*J):
  # fix mask
  mask_i, mask_j = n2ij(m)
  for n in range(I*J):
  # swipe entire domain but excluding boundary
      i, j = n2ij(n)
      if abs(i-mask_i)==0 and abs(j-mask_j)==0:
         sparse[m][n] = -4
      elif abs(i-mask_i)==0 and abs(j-mask_j)==1:
         sparse[m][n] = 1
      elif abs(i-mask_i)==1 and abs(j-mask_j)==0:
         sparse[m][n] = 1
      else:
         sparse[m][n] = 0



# plot 2D sparse matrx of 2D domain
for m in range(I*J):
   for n in range(I*J):
        print("%+d" % sparse[m][n], end ='')
   print("")
