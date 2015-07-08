"""
Extracts a pair of features from MNIST SVD training set
"""
f_rd = open ("/scratch/cluster/aditya/DBN_research/rp_deep/original_code/F.txt")
X = []
Y = []

for line in f_rd: 
    words = line.strip().split()
    X += [float(words[0])] #.append(words[0]+'\n')
    Y += [float(words[1])] #.append(words[1]+'\n')

min_x = min(X)
max_x = max(X)
X_norm = [(x-min_x)/(max_x-min_x) for x in X]
min_y = min(Y)
max_y = max(Y)
Y_norm = [(y-min_y)/(max_y-min_y) for y in Y]

f_wrx = open("del_X.txt", "w")
for x in X_norm:
  f_wrx.write("%s\n" % x)
f_wrx.close()

f_wry = open("del_Y.txt", "w")
for y in Y_norm:
  f_wry.write("%s\n" % y)
f_wry.close()
