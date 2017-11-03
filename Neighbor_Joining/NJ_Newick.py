#! /usr/bin/env python3
##################
#####Neighbor Joining Algorithm with Newick output for dissimilarity matrix
#####Yifan Ethan Chen
#####Math 127 Groupwork 2
#####10/11/17 Cal
import numpy as np
import copy
#####input info
n = 5     #change to the matrix n
input_order = 'A,B,C,D,E'
a = np.array([[0,0.31,1.01,0.75,1.03],
             [0.31,0,1.00,0.69,0.90],
             [1.01,1.00,0,0.61,0.42],
             [0.75,0.69,0.61,0,0.37],
             [1.03,0.90,0.42,0.37,0]])
'''a = np.array([[0,0.31,1.01,0.75,1.03,0.44],
             [0.31,0,1.00,0.69,0.90,0.51],
             [1.01,1.00,0,0.61,0.42,0.62],
             [0.75,0.69,0.61,0,0.37,0.73],
             [1.03,0.90,0.42,0.37,0,0.25],
             [0.44,0.51,0.62,0.73,0.25,0]])'''
##########
#rep=[]
allresult=[]
order_for_output = input_order.split(',')
order = ['_('+i+')_' for i in order_for_output]


#####step1
def compMij(array):
#####compute R and Rn+Rm matrix
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
#####compute Mij matrix
    Mij=(n-2)*array-R
    print ('Mij table:')
    print (order_for_output)
    print (Mij)
#####find the min location 
    mini, minj = np.where(Mij == np.min(Mij))[0]
    print ('We\'re gonna join ' + order_for_output[mini] + ' and ' + order_for_output[minj] + ' in this step and the node is (' + order_for_output[mini] + ',' + order_for_output[minj] +')')
    order.append('_('+order[mini]+','+order[minj]+')_')
    order_for_output.append('('+order_for_output[mini]+','+order_for_output[minj]+')')
    return Rl, mini, minj, order


#####step2
def calGS(array, Rl, mini, minj):
#####calculate the group distance
    GSi=1/(n-2)*(Rl[mini]-array[mini][minj])
    GSj=1/(n-2)*(Rl[minj]-array[mini][minj])
#####3 point formula
    vSi=np.around(1/2*(array[mini][minj]+GSi-GSj), decimals=4)
    vSj=array[mini][minj]-vSi
    print ('distance to the parent: ' + order_for_output[mini] +': ' + str(vSi) + ' ;' + order_for_output[minj] +': ' + str(vSj))
    #rep.append((order[mini]+':',str(vSi),',', order[minj]+':',str(vSj),','))
    allresult.append((order[mini],order[mini][1:-1]+':'+str(vSi)))
    allresult.append((order[minj],order[minj][1:-1]+':'+str(vSj)))


#####step3
def newdisa(array, mini, minj):
#####calculate the new distance
    addeddis=[]
    for i in range(len(order)-1):
        if i == mini or i == minj:
            pass
        else:
            dis=1/2*(array[i][mini]+array[i][minj]-array[mini][minj])
            addeddis.append(dis)
#####create new dis matrix
    array=np.delete(array, mini, 0)
    array=np.delete(array, minj-1, 0)
    array=np.delete(array, mini, 1)
    newarr=np.delete(array, minj-1, 1)
    del order[mini], order[minj-1], order_for_output[mini], order_for_output[minj-1]
    newarr=np.c_[newarr,addeddis]
    addeddis.append(0)
    newarr=np.r_['0,2',newarr,addeddis]
    print ('updated dis matrix')
    print (order_for_output)
    print (newarr)
    return order, newarr


#####endgame
def endgame(array):
    dis0=1/2*(array[0][1]+array[0][2]-array[1][2])
    dis1=1/2*(array[0][1]+array[1][2]-array[0][2])
    dis2=1/2*(array[1][2]+array[0][2]-array[0][1])
    #rep.append((order[0]+':',str(dis0),',',order[1]+':',str(dis1),',',order[2]+':',str(dis2),','))
    allresult.append((order[0],order[0][1:-1]+':'+str(dis0)))
    allresult.append((order[1],order[1][1:-1]+':'+str(dis1)))
    allresult.append((order[2],order[2][1:-1]+':'+str(dis2)))
    final_result = order[0]+','+order[1]+','+order[2]
    #print ('final results:')
    #print (rep)
    return final_result


#####convert to newick format
def do_newick(newick_seq):
    strn = 0
    bracket_num = 0
    for str_string in newick_seq:
        try:
            if str_string == '_' and newick_seq[strn+1] == '(':
                if bracket_num == 0:
                    bracket_num_left = copy.deepcopy(strn)
                bracket_num += 1
            elif str_string == '_' and newick_seq[strn-1] == ')':
                bracket_num -= 1
                if bracket_num == 0:
                    bracket_num_right = copy.deepcopy(strn)
                    break
        except:
            if str_string == '_' and newick_seq[strn-1] == ')':
                bracket_num -= 1
                if bracket_num == 0:
                    bracket_num_right = copy.deepcopy(strn)
                    break
        strn += 1
#####find the replacement
    for group in allresult:
        if group[0] == newick_seq[bracket_num_left:bracket_num_right+1]:
            group_1 = group[1]
            break
    newick_seq = newick_seq.replace(newick_seq[bracket_num_left:bracket_num_right+1], group_1)
    return newick_seq


#####
for times in range(n):
    Rl, mini, minj, order=compMij(a)
    calGS(a,Rl,mini, minj)
    order, a=newdisa(a,mini, minj)
    if len(order)==3:
        break
final_result = endgame(a)
while '_' in final_result:
    final_result = do_newick(final_result)
for i in input_order.split(','):
    final_result = final_result.replace('('+i+')', i)
print ('The final result in Newick Fomat:')
print ('('+final_result+');')