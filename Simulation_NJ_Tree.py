#! /usr/bin/env python3
##################
#####Simulate the evolution data, using Neighbor Joining Algorithm and draw the tree
#####Yifan Ethan Chen
#####Math 127 Groupwork 2
#####10/20/17 Cal
from scipy.io import loadmat
from Bio import Phylo
from io import StringIO
import numpy as np
import copy, sys, random, matplotlib
import matplotlib.pyplot as plt


#####set up parameters
mutation_rate = 0.01   #set up your mutation rate
time_step = 1   #set up your time step
sequence_want = 10   #the total sequences you want 
report_dict = {}   #all new sequences will be here {'sequence': sequence number}

#the tree should be like this:
#             root
#             /
#            /-----start sequence
#           /
#          /----- sequence 1
#         /
#        /----- sequence 2
#       /
#      /----- sequence 3
#     /
#   sequence 4


#####some functions
####simulate with Jukes-Cantor matrix
def JC_simulate(sequence, mutation_rate, time_step):
    for i in range(time_step):
        result_seq = []
        for base in sequence:
            if random.random() < mutation_rate:
                base_list = ['A', 'T', 'C', 'G']
                base_list.remove(base)
                random.shuffle(base_list)
                new_base = base_list[random.randint(0, 2)]
                result_seq.append(new_base)
            else:
                result_seq.append(base)
        sequence = ''.join(result_seq)
    return sequence


####for sorting new sequences by their evolution time
def evolution_num(s): 
    return int(s[0])


####for drawing tree
def plot_tree(treedata, output_file):
    handle = StringIO(treedata)
    tree = Phylo.read(handle, 'newick')
    matplotlib.rc('font', size=6)
    fig = plt.figure(figsize=(10, 20), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    plt.savefig(output_file+'.png')
    #plt.savefig(output_file+'.pdf', format='PDF')


####for neighbor joining algorithm
###step1
def compMij(array):
###compute R and Rn+Rm matrix
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
###compute Mij matrix
    Mij=(n-2)*array-R
    print ('Mij table:')
    print (order_for_output)
    print (Mij)
###find the min location 
    mini, minj = np.where(Mij == np.min(Mij))[0]
    print ('We\'re gonna join ' + order_for_output[mini] + ' and ' + order_for_output[minj] + ' in this step and the node is (' + order_for_output[mini] + ',' + order_for_output[minj] +')')
    order.append('_('+order[mini]+','+order[minj]+')_')
    order_for_output.append('('+order_for_output[mini]+','+order_for_output[minj]+')')
    return Rl, mini, minj, order

###step2
def calGS(array, Rl, mini, minj):
###calculate the group distance
    GSi=1/(n-2)*(Rl[mini]-array[mini][minj])
    GSj=1/(n-2)*(Rl[minj]-array[mini][minj])
###3 point formula
    vSi=np.around(1/2*(array[mini][minj]+GSi-GSj), decimals=4)
    vSj=array[mini][minj]-vSi
    print ('distance to the parent: ' + order_for_output[mini] +': ' + str(vSi) + ' ;' + order_for_output[minj] +': ' + str(vSj))
    allresult.append((order[mini],order[mini][1:-1]+':'+str(vSi)))
    allresult.append((order[minj],order[minj][1:-1]+':'+str(vSj)))

###step3
def newdisa(array, mini, minj):
###calculate the new distance
    addeddis=[]
    for i in range(len(order)-1):
        if i == mini or i == minj:
            pass
        else:
            dis=1/2*(array[i][mini]+array[i][minj]-array[mini][minj])
            addeddis.append(dis)
###create new dis matrix
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

###endgame
def endgame(array):
    dis0=1/2*(array[0][1]+array[0][2]-array[1][2])
    dis1=1/2*(array[0][1]+array[1][2]-array[0][2])
    dis2=1/2*(array[1][2]+array[0][2]-array[0][1])
    allresult.append((order[0],order[0][1:-1]+':'+str(dis0)))
    allresult.append((order[1],order[1][1:-1]+':'+str(dis1)))
    allresult.append((order[2],order[2][1:-1]+':'+str(dis2)))
    final_result = order[0]+','+order[1]+','+order[2]
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
####find the replacement
    for group in allresult:
        if group[0] == newick_seq[bracket_num_left:bracket_num_right+1]:
            group_1 = group[1]
            break
    newick_seq = newick_seq.replace(newick_seq[bracket_num_left:bracket_num_right+1], group_1)
    return newick_seq


#####here we use dnt sequence in file flhivdata.mat
data = loadmat('flhivdata.mat')
target_sequence = str(data['dnt'])[3:-2]   #here we use dnt sequence in file flhivdata.mat
target_sequence_len = len(target_sequence)
report_dict.setdefault(target_sequence, 0)
target_sequence = JC_simulate(target_sequence, mutation_rate, time_step)

#####generate new sequences into report_dict
for i in range(1, sequence_want):
    report_sequence = JC_simulate(target_sequence, mutation_rate, time_step)
    while report_sequence in report_dict:
        report_sequence = JC_simulate(target_sequence, mutation_rate, time_step)
    report_dict.setdefault(report_sequence, i)
    target_sequence = JC_simulate(target_sequence, mutation_rate, time_step)
    
#####write new sequences into a file named new_sequence.fa
fout = open('new_sequence.fa', 'w')
for sequence in report_dict:
    fout.write('>'+str(report_dict[sequence])+'\n')
    seqlist = [sequence[i:i+80]+'\n' for i in range(0,len(sequence),80)]
    for subseq in seqlist:
        fout.write(subseq)
fout.close()

#####true tree for simulation data
simulate_newick = '0'
for i in range(1,sequence_want-2):
    simulate_newick = '('+simulate_newick+','+str(i)+')'
simulate_newick = '('+simulate_newick+','+str(sequence_want-2)+','+str(sequence_want-1)+');'
print ('Newick format for the true tree without distance:')
print (simulate_newick)
plot_tree(simulate_newick, 'true_tree')

#####calculate Jukes-Cantor distance (dissimilarity) matrix
sequence_list = []
for sequence in report_dict:
    sequence_list.append((str(report_dict[sequence]), sequence))
sequence_list = sorted(sequence_list, key = evolution_num)
input_order = ','.join([str(i) for i in range(sequence_want)])
JC_dis = []
for i in sequence_list:
    sequence_i = i[1]
    for j in sequence_list:
        if i == j:
            JC_dis.append(0)
        else:
            sequence_j = j[1]
            mutation_site_num = 0
            for base_num in range(target_sequence_len):
                if sequence_j[base_num] != sequence_i[base_num]:
                    mutation_site_num += 1
            distance = -3/4*np.log(1-4/3*mutation_site_num/target_sequence_len)
            JC_dis.append(distance)
JC_dis_matrix = np.reshape(np.asarray(JC_dis),(sequence_want,sequence_want))

#####doing Neighbor Joining Algorithm
n =  sequence_want   #change to the matrix n
a = copy.deepcopy(JC_dis_matrix)

##########
#rep=[]
allresult=[]
order_for_output = input_order.split(',')
order = ['_('+i+')_' for i in order_for_output]
print ('The Jukes-Cantor distance (dissimilarity) matrix for this simulation as follows:')
print (order_for_output)
print (JC_dis_matrix)

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
final_result = '('+final_result+');'
print ('The final result in Newick Fomat:')
print (final_result)
plot_tree(final_result, 'NJ_tree')

sys.exit()