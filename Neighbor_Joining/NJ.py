#! /usr/bin/env python3
##################
#####Neighbor Joining Algorithm for dissimilarity matrix
#####Yifan Ethan Chen
#####Math 127 Groupwork 1
#####09/15/17 Cal
import numpy as np
######input info
n = 5
a = np.array([[0,0.31,1.01,0.75,1.03],
             [0.31,0,1.00,0.69,0.90],
             [1.01,1.00,0,0.61,0.42],
             [0.75,0.69,0.61,0,0.37],
             [1.03,0.90,0.42,0.37,0]])
##########
rep=[]
order = [str(i) for i in list(np.arange(1,6))]
######step1
def compMij(array):
######compute R and Rn+Rm matrix
    Rl=[]
    for row in array:
        r=np.sum(row)
        Rl.append(r)
    R=[]
    for i in range(len(order)):
        for j in range(len(order)):
            if i == j:
                R.append(0)
            else:
                R.append(Rl[i]+Rl[j])
    R=np.reshape(np.asarray(R),(len(order),len(order)))
    #print (Rl)
    #print (R)
######compute Mij matrix
    Mij=(n-2)*array-R
    print ('Mij table:')
    print (order)
    print (Mij)
######find the min location 
    mini, minj = np.where(Mij == np.min(Mij))[0]
    print ('We\'re gonna join ' + order[mini] + ' and ' + order[minj] + ' in this step and the node is (' + order[mini] + ',' + order[minj] +')')
    order.append('('+order[mini]+','+order[minj]+')')
    return Rl, mini, minj, order
#######step2
def calGS(array, Rl, mini, minj):
######calculate the group distance
    GSi=1/(n-2)*(Rl[mini]-array[mini][minj])
    GSj=1/(n-2)*(Rl[minj]-array[mini][minj])
######3 point formula
    vSi=np.around(1/2*(array[mini][minj]+GSi-GSj), decimals=4)
    vSj=array[mini][minj]-vSi
    print ('distance to the parent: ' + order[mini] +': ' + str(vSi) + ' ;' + order[minj] +': ' + str(vSj))
    rep.append((order[mini],order[minj],str(vSi),str(vSj)))
#######step3
def newdisa(array, mini, minj):
######calculate the new distance
    addeddis=[]
    for i in range(len(order)-1):
        if i == mini or i == minj:
            pass
        else:
            dis=1/2*(array[i][mini]+array[i][minj]-array[mini][minj])
            addeddis.append(dis)
######create new dis matrix
    array=np.delete(array, mini, 0)
    array=np.delete(array, minj-1, 0)
    array=np.delete(array, mini, 1)
    newarr=np.delete(array, minj-1, 1)
    del order[mini], order[minj-1]
    newarr=np.c_[newarr,addeddis]
    addeddis.append(0)
    newarr=np.r_['0,2',newarr,addeddis]
    print ('updated dis matrix')
    print (order)
    print (newarr)
    return order, newarr
######endgame
def endgame(array):
    dis0=1/2*(array[0][1]+array[0][2]-array[1][2])
    dis1=1/2*(array[0][1]+array[1][2]-array[0][2])
    dis2=1/2*(array[1][2]+array[0][2]-array[0][1])
    rep.append((order[0],order[1],order[2],str(dis0),str(dis1),str(dis2)))
    print ('final results:')
    print (rep)
#######
for times in range(n):
    Rl, mini, minj, order=compMij(a)
    calGS(a,Rl,mini, minj)
    order, a=newdisa(a,mini, minj)
    if len(order)==3:
        break
endgame(a)