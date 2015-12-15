import os
from pylab import plot,show
from numpy import vstack,array
from numpy.random import rand
from scipy.cluster.vq import kmeans,vq
from scipy import stats
from scipy.stats import linregress

def read_in_dotplot_files(filein):
    fin=open(filein)
    off_dis=int(filein.split('.')[-1])
    out=[[],[]]
    for line in fin:
        pin=line.strip().split()
        out[0].append(int(pin[0])-off_dis)
        out[1].append(int(pin[1]))
    fin.close()
    #out_range=min([max(out[0]),max(out[1])])+1
    #out1=[[],[]]
    #for x in range(len(out[0])):
    #    if out[0][x]<out_range and out[1][x]<out_range:
    #        out1[0].append(out[0][x])
    #        out1[1].append(out[1][x])            
    return out

def tranform_diagnal_to_horizonal(out):
    #transformation: x2=x+y, y2=x-y
    #out=read_in_dotplot_files(filein)
    out2=[[],[]]
    for x in range(len(out[0])):
        out2[0].append(out[0][x]+out[1][x])
        out2[1].append(out[0][x]-out[1][x])
    return out2

def tranform_horizonal_to_diagnal(out):
    out2=[[],[]]
    for x in range(len(out[0])):
        out2[0].append(int(float(out[0][x]+out[1][x])/2.0))
        out2[1].append(int(float(out[0][x]-out[1][x])/2.0))
    return out2

def cluster_numbers(list,dis_cff):
    out=[[]]
    for k1 in sorted(list):
        if out[-1]==[]:
            out[-1].append(k1)
        else:
            if k1-out[-1][-1]<dis_cff:
                out[-1].append(k1)
            else:
                out.append([k1])
    return out

lenght_cff=100
dots_num_cff=20
clu_dis_cff=5
def cluter_to_diagnal(out):
    #search for clusters according to diagnal. based on dis of a dot to diagnal:
    out2=tranform_diagnal_to_horizonal(out)
    out2_hash1={}
    for k1 in range(len(out2[1])):
        if not out2[1][k1] in out2_hash1.keys():
            out2_hash1[out2[1][k1]]=[]
        out2_hash1[out2[1][k1]].append(out2[0][k1])
    cluster_a=cluster_numbers(out2_hash1.keys(),clu_dis_cff)
    cluster_b=[]
    for x in cluster_a:
        cluster_b.append([])
        for y in x:
            cluster_b[-1]+=out2_hash1[y]
    cluster_2_a=[]
    cluster_2_b=[]
    cluster_2_rest=[[],[]]
    for x in range(len(cluster_b)):
        if max(cluster_b[x])-min(cluster_b[x])>lenght_cff and len(cluster_b[x])>dots_num_cff:
            cluster_2_a.append([])
            for y in cluster_a[x]:
                cluster_2_a[-1]+=[y for i in out2_hash1[y]]
            cluster_2_b.append(cluster_b[x])
        else:
            cluster_2_rest[0]+=cluster_a[x]
            cluster_2_rest[1]+=cluster_b[x]
    diagnal_segs=[]
    for x in range(len(cluster_2_a)):
        diagnal_segs.append(tranform_horizonal_to_diagnal([cluster_2_b[x],cluster_2_a[x]]))
    #then check for possible anti-diagnal patterns (inversions)
    out2_hash2={}
    for x in cluster_2_rest[0]:
        for y in out2_hash1[x]:
            if not y in out2_hash2.keys():
                out2_hash2[y]=[]
            out2_hash2[y].append(x)
    cluster_c=cluster_numbers(cluster_2_rest[1],clu_dis_cff)
    cluster_d=[]
    for x in cluster_c:
        cluster_d.append([])
        for y in x:
            cluster_d[-1]+=out2_hash2[y]
    cluster_2_c=[]
    cluster_2_d=[]
    for x in range(len(cluster_d)):
        if max(cluster_d[x])-min(cluster_d[x])>lenght_cff and len(cluster_d[x])>dots_num_cff:
            cluster_2_c.append([])
            for y in cluster_c[x]:
                cluster_2_c[-1]+=[y for i in out2_hash2[y]]
            cluster_2_d.append(cluster_d[x])
    anti_diagnal_segs=[]
    for x in range(len(cluster_2_d)):
        anti_diagnal_segs.append(tranform_horizonal_to_diagnal([cluster_2_c[x],cluster_2_d[x]]))
    return [diagnal_segs,anti_diagnal_segs]

def K_2_means_cluster(data_list):
    #eg of data_list: [[],[]]
    out=[]
    if not data_list==[[],[]]:
        linear_info_ini=linregress(data_list[1],data_list[0]) #apply linear regression on original data
        array_diagnal=array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
        std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
        whitened = whiten(array_diagnal)
        centroids, variance=kmeans(whitened,2)
        idx,_= vq(whitened,centroids)
        #group1=[[int(i*std_rec[0]) for i in whitened[idx==0,0]],[int(i*std_rec[1]) for i in whitened[idx==0,1]]]
        #group2=[[int(i*std_rec[0]) for i in whitened[idx==1,0]],[int(i*std_rec[1]) for i in whitened[idx==1,1]]]
        group1=[[int(i) for i in array_diagnal[idx==0,0]],[int(i) for i in array_diagnal[idx==0,1]]]
        linear_info_1=linregress(group1[1],group1[0]) #apply linear regression on group1
        if linear_info_1[4]<linear_info_ini[4]/10 and abs(linear_info_1[0]-1)<abs(linear_info_ini[0]-1)/10:
            out.append(group1)
        group2=[[int(i) for i in array_diagnal[idx==1,0]],[int(i) for i in array_diagnal[idx==1,1]]]
        linear_info_2=linregress(group2[1],group2[0]) #apply linear regression on group2
        if linear_info_2[4]<linear_info_ini[4]/10 and abs(linear_info_2[0]-1)<abs(linear_info_ini[0]-1)/10:
            out.append(group2)
        if out==[]:
            out.append(data_list)
    else:
        out.append([[],[]])
    return out

def K_2_anti_diag_cluster(data_list):
    #eg of data_list: [[],[]]
    out=[]
    if not data_list==[[],[]]:
        linear_info_ini=linregress(data_list[1],data_list[0]) #apply linear regression on original data
        array_diagnal=array([[data_list[0][x],data_list[1][x]] for x in range(len(data_list[0]))])
        std_rec=[scipy.std(data_list[0]),scipy.std(data_list[1])]
        whitened = whiten(array_diagnal)
        centroids, variance=kmeans(whitened,2)
        idx,_= vq(whitened,centroids)
        #group1=[[int(i*std_rec[0]) for i in whitened[idx==0,0]],[int(i*std_rec[1]) for i in whitened[idx==0,1]]]
        #group2=[[int(i*std_rec[0]) for i in whitened[idx==1,0]],[int(i*std_rec[1]) for i in whitened[idx==1,1]]]
        group1=[[int(i) for i in array_diagnal[idx==0,0]],[int(i) for i in array_diagnal[idx==0,1]]]
        linear_info_1=linregress(group1[1],group1[0]) #apply linear regression on group1
        if linear_info_1[4]<linear_info_ini[4]/10 and abs(linear_info_1[0]+1)<abs(linear_info_ini[0]+1)/10:
            out.append(group1)
        group2=[[int(i) for i in array_diagnal[idx==1,0]],[int(i) for i in array_diagnal[idx==1,1]]]
        linear_info_2=linregress(group2[1],group2[0]) #apply linear regression on group2
        if linear_info_2[4]<linear_info_ini[4]/10 and abs(linear_info_2[0]+1)<abs(linear_info_ini[0]+1)/10:
            out.append(group2)
        if out==[]:
            out.append(data_list)
    else:
        out.append([[],[]])
    return out

def X_means_cluster(data_list):
    temp_result=K_2_means_cluster(data_list)
    if len(temp_result)==1:
        return temp_result[0]
    else:
        return [X_means_cluster(i) for i in temp_result]

def X_anti_diag_cluster(data_list):
    temp_result=K_2_anti_diag_cluster(data_list)
    if len(temp_result)==1:
        return temp_result[0]
    else:
        return [X_means_cluster(i) for i in temp_result]

def cluster_qc(out):
    diagnal_info=cluter_to_diagnal(out)
    out2=[]
    for x in diagnal_info:
        out2.append([[],[]])
        for y in x:
            out2[-1][0]+=y[0]
            out2[-1][1]+=y[1]
    if not out2[0]==[[],[]]:
        out_diagonal=X_means_cluster(out2[0])
        if len(out_diagonal)==2 and len(out_diagonal[0])>2:
            out_diagonal=[out_diagonal]
    else:
        out_diagonal=[[],[]]
    if not out2[1]==[[],[]]:
        out_anti_diagonal=X_anti_diag_cluster(out2[1])
        if len(out_anti_diagonal)==2 and len(out_anti_diagonal[0])>2:
            out_anti_diagonal=[out_anti_diagonal]
    else:
        out_anti_diagonal=[[],[]]
    return [out_diagonal,out_anti_diagonal]

filein='Triple.chr1_144894229_144906851.dotplot.ref.left.935c09dd_131293_18613.start.0'
filein='Triple.chr1_17008871_17013610.dotplot.ref.935c09dd_40088_1325.start.9'
filein='Triple.chr1_72766298_72811876.dotplot.ref.left.45a41f96_105491_0.start.0'
out=read_in_dotplot_files(filein)
out_final=cluster_qc(out)

