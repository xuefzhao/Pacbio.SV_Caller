#!/usr/bin/env python

#!python
#sys.argv=command.split()
import os
import re
import sys
if len(sys.argv)<2:
	print_read_me()
else:
	function_name=sys.argv[1]
	window_size_1=10
	def read_vcf_in(filein):
		data_hash={}
		fin=open(filein)
		for line in fin:
			pin=line.strip().split()
			if not pin[0][0]=='#':
				if not pin[2] in data_hash.keys():
					data_hash[pin[2]]=[]
				if pin[2][-1]=='C':
					if '[' in pin[4] or ']' in pin[4]:
						data_hash[pin[2]].append([pin[0],pin[1],pin[4]])
				else:
					data_hash[pin[2]].append(pin[4:])
		fin.close()
		return data_hash
	def produce_structure_to_vali(data_hash,k1):
		out=[]
		for x in range(len(data_hash[k1])/2):
			if data_hash[k1][x][2][0] in ['a','t','g','c','n','A','T','G','C','N'] and data_hash[k1][x][2][-1] in ['[',']']:
				if '[' in data_hash[k1][x][2]:
					block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1])-flank_length,int(data_hash[k1][x][1])]
					block_b_temp=data_hash[k1][x][2].split('[')[1].split(':')
					block_b=[block_b_temp[0],int(block_b_temp[1]),int(block_b_temp[1])+flank_length]
					out.append(['a','b',block_a,block_b])
				elif ']' in data_hash[k1][x][2]:
					block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1])-flank_length,int(data_hash[k1][x][1])]
					block_b_temp=data_hash[k1][x][2].split(']')[1].split(':')
					block_b=[block_b_temp[0],int(block_b_temp[1])-flank_length,int(block_b_temp[1])]
					out.append( ['a','b^',block_a,block_b])
			if data_hash[k1][x][2][-1] in ['a','t','g','c','n','A','T','G','C','N'] and data_hash[k1][x][2][0] in ['[',']']:
				if '[' in data_hash[k1][x][2]:
					block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1]),int(data_hash[k1][x][1])+flank_length]
					block_b_temp=data_hash[k1][x][2].split('[')[1].split(':')
					block_b=[block_b_temp[0],int(block_b_temp[1]),int(block_b_temp[1])+flank_length]
					out.append( ['a','b^',block_a,block_b])
				elif ']' in data_hash[k1][x][2]:
					block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1]),int(data_hash[k1][x][1])+flank_length]
					block_b_temp=data_hash[k1][x][2].split(']')[1].split(':')
					block_b=[block_b_temp[0],int(block_b_temp[1])-flank_length,int(block_b_temp[1])]
					out.append( ['a','b',block_a,block_b])
		return out
	def complementary(seq):
		seq2=[]
		for i in seq:
				if i in 'ATGCN':
						seq2.append('ATGCN'['TACGN'.index(i)])
				elif i in 'atgcn':
						seq2.append('atgcn'['tacgn'.index(i)])
		return ''.join(seq2)
	def reverse(seq):
		return seq[::-1]
	def chop_pacbio_read_long(bps):
		flank_length=0
		block_length={}
		len_cff=bps[-1]-bps[1]
		fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
		out=[]
		out2=[]
		out3=[]
		for line in fbam:
				pbam=line.strip().split()
				#out3.append(int(pbam[3])-int(bps[1])+flank_length)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
						align_start=align_info[0]
						miss_bp=align_info[1]
						align_pos=int(pbam[3])
						target_read=pbam[9][align_start:]
						if len(target_read)>flank_length+len_cff:
							out.append(target_read[:len_cff+2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
							#out3.append(pbam[:9])
		fbam.close()
		return [out,out2,out3]
	def chop_x_for_fasta(x):
		if len(x)>60:
				lines=len(x)/60
				out=[]
				for k in range(lines):
						out.append(x[(k*60):((k+1)*60)])
				out.append(x[((k+1)*60):])
				return out
		else:
				return [x]
	def alt_sv_to_list(alt_sv):
		out=[[],[]]
		for x in alt_sv.split('/')[0]:
				if not x=='^':
						out[0].append(x)
				else:
						out[0][-1]+=x
		for x in alt_sv.split('/')[1]:
				if not x=='^':
						out[1].append(x)
				else:
						out[1][-1]+=x
		return out
	def ref_hash_modi(ref_hash):
		out={}
		for x in ref_hash.keys():
				out[x]=ref_hash[x]
				out[x+'^']=reverse(complementary(ref_hash[x]))
		return out
	def Pacbio_prodce_ref_alt_short(ref,flank_length,txt_file):
		fin=open(txt_file)
		pin=fin.readline().strip().split()
		fin.close()
		ref_sv=pin[0]
		ref_sv='/'.join(['1'+ref_sv.split('/')[0]+'2','1'+ref_sv.split('/')[0]+'2'])
		alt_sv=pin[1]
		alt_sv='/'.join(['1'+alt_sv.split('/')[0]+'2','1'+alt_sv.split('/')[1]+'2'])
		chrom=pin[2]
		bps=pin[2:]
		bps=[bps[0]]+[str(int(bps[1])-flank_length)]+bps[1:]+[str(int(bps[-1])+flank_length)]
		ref_hash={}
		rec=0
		for x in ref_sv.split('/')[0]:
				rec+=1
				fref=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,chrom,int(bps[rec]), int(bps[rec+1])))
				fref.readline().strip().split()
				seq=''
				while True:
						pref=fref.readline().strip().split()
						if not pref: break
						seq+=pref[0]
				fref.close()
				ref_hash[x]=seq
		ref_hash=ref_hash_modi(ref_hash)
		alt_sv_list=alt_sv_to_list(alt_sv)
		ref_seq=''.join([ref_hash[x] for x in ref_sv.split('/')[0]])
		alt_1_seq=''.join([ref_hash[x] for x in alt_sv_list[0]])
		alt_2_seq=''.join([ref_hash[x] for x in alt_sv_list[1]])
		return [ref_seq,alt_1_seq,alt_2_seq]
	def Pacbio_prodce_ref_alt_long(ref,new_structure_single,txt_file):
		#eg: new_structure_single=new_structure[0]
		x=new_structure_single
		struc_hash={}
		struc_hash[x[0][0]]=''
		struc_hash[x[0][0]+'^']=''
		fread1=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,x[2][0],x[2][1],x[2][2]))
		pin=fread1.readline().strip().split()
		for line in fread1:
			pin=line.strip().split()
			struc_hash[x[0][0]]+=pin[0]
		fread1.close()
		struc_hash[x[0][0]+'^']=reverse(complementary(struc_hash[x[0][0]]))
		struc_hash[x[1][0]]=''
		struc_hash[x[1][0]+'^']=''
		fread1=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,x[3][0],x[3][1],x[3][2]))
		pin=fread1.readline().strip().split()
		for line in fread1:
			pin=line.strip().split()
			struc_hash[x[1][0]]+=pin[0]
		fread1.close()
		struc_hash[x[1][0]+'^']=reverse(complementary(struc_hash[x[1][0]]))
		fo1=open('.'.join(txt_file.split('.')[:-1])+'.alt.fa','w')
		print >>fo1, '_'.join([str(i) for i in [x[0],x[1]]+x[2]+x[3]])
		new_seq=''
		for y in x[:2]:
			new_seq+=struc_hash[y]
		for y in chop_x_for_fasta(new_seq):
			print >>fo1,y
		fo1.close()
		os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,x[2][0],x[2][1],x[2][1]+2*flank_length,'.'.join(txt_file.split('.')[:-1])+'.ref1.fa'))
		os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,x[3][0],x[3][1],x[3][1]+2*flank_length,'.'.join(txt_file.split('.')[:-1])+'.ref2.fa'))
		return ['.'.join(txt_file.split('.')[:-1])+'.alt.fa','.'.join(txt_file.split('.')[:-1])+'.ref1.fa','.'.join(txt_file.split('.')[:-1])+'.ref2.fa']
	def key_name_define(vcf_data_hash,k3):
		if '_'.join(k3+['C']) in vcf_data_hash.keys():
			key_name='_'.join(k3+['C'])
		elif '_'.join(k3+['S']) in vcf_data_hash.keys():
			key_name='_'.join(k3+['S'])
		else:
			key_name='Null'
		return key_name
	def txt_rec_file_writing(out_path,SVelter_Name,info):
		txt_file=out_path+SVelter_Name+'.txt'
		fo=open(txt_file,'w')
		print >>fo, '\t'.join(info)
		fo.close()
	def eu_dis_calcu_1(fo_ref,rsquare_ref,align_off,delta):
		#this function calculates total distance of dots to diagnal, with direction considered; smaller = better ;  especially designed for duplcations;
		temp_data1=[[],[]]
		for x in fo_ref:
			temp_data1[0].append(int(x[0])-align_off)
			temp_data1[1].append(int(x[1]))
		if not temp_data1[0]==[]:
			out=sum([temp_data1[1][x]-temp_data1[0][x] for x in range(len(temp_data1[0]))])
			rsquare_ref.append(abs(out))
	def eu_dis_calcu_2(fo_ref,rsquare_ref,align_off,delta):
		#this function calculates total distance of dots to diagnal; smaller = better ; not suitable for duplications
		temp_data1=[[],[]]
		for x in fo_ref:
			temp_data1[0].append(int(x[0])-align_off)
			temp_data1[1].append(int(x[1]))
		if not temp_data1[0]==[]:
			out=sum([abs(temp_data1[1][x]-temp_data1[0][x]) for x in range(len(temp_data1[0]))])
			rsquare_ref.append(abs(out))
	def eu_dis_calcu_3a(fo_ref,rsquare_ref,align_off,delta):
		#this function calculates total distance of dots to diagnal; smaller = better
		fo1=open(fo_ref+'.txt')
		temp_data1=[[],[]]
		for line in fo1:
			po1=line.strip().split()
			temp_data1[0].append(int(po1[0])-align_off)
			temp_data1[1].append(int(po1[1]))
		fo1.close()
		rec1=1
		for x in range(len(temp_data1[0])):
			if abs(temp_data1[1][x]-temp_data1[0][x])<float(temp_data1[0][x])/10.0:
				rec1+=1
		rsquare_ref.append(float(len(temp_data1[0]))/float(rec1))
	def eu_dis_calcu_3(fo_ref,rsquare_ref,align_off,delta):
		#count total #dots/#dots locates close to diagnal; smaller = better
		#either diagnal or reverse diganal
		fo1=open(fo_ref+'.txt')
		temp_data1=[[],[]]
		for line in fo1:
			po1=line.strip().split()
			temp_data1[0].append(int(po1[0])-align_off)
			temp_data1[1].append(int(po1[1]))
		fo1.close()
		temp_data1.append(temp_data1[1][::-1])
		rec1=1
		for x in range(len(temp_data1[0])):
			if abs(temp_data1[1][x]-temp_data1[0][x])<float(temp_data1[0][x])/10.0:
				rec1+=1
		rec2=1
		for x in range(len(temp_data1[0])):
			if abs(temp_data1[2][x]-temp_data1[0][x])<float(temp_data1[0][x])/10.0:
				rec2+=1
		rsquare_ref.append(float(len(temp_data1[0]))/float(rec1))
	def dup_decide(structure):
		flag=0
		for x in structure:
			if not x=='^':
				if structure.count(x)>1:
					flag+=1
		return flag
	def chromos_readin(ref):
		fin=open(ref+'.fai')
		chromos=[]
		for line in fin:
				pin=line.strip().split()
				chromos.append(pin[0])
		fin.close()
		return chromos
	def remove_files_short(txt_file):
		for x in os.listdir('/'.join(txt_file.split('/')[:-1])):
			if txt_file.split('/')[-1].replace('.txt','') in x:
				if not x.split('.')[-1]=='png':
					if not x.split('.')[-1]=='rsquare':
						if not x==txt_file.split('/')[-1]:
							if not x.split('.')[-2]=='start':
								os.system(r'''rm %s'''%('/'.join(txt_file.split('/')[:-1])+'/'+x))
	def chop_x_for_fasta(x):
		if len(x)>60:
			lines=len(x)/60
			out=[]
			for k in range(lines):
				out.append(x[(k*60):((k+1)*60)])
			out.append(x[((k+1)*60):])
			return out
		else:
			return [x]
	def cigar_modify(cigar):
		temp=[[]]
		for x in cigar:
			if ord(x)>64 and ord(x)<91:
				temp[-1].append(x)
				temp.append([])
			else:
				if temp[-1]==[]:
					temp[-1].append(x)
				else:
					temp[-1][-1]+=x
		temp.remove([])
		for x in temp:
			if x[1]=='S' and not int(x[0])>min_length:
				x[1]='I'
		return ''.join(''.join(y) for y in temp)
	def cigar_integrate(pin):
		pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
		cigar=cigar_modify(pin[5])
		pos_start=int(pin[3])
		cigars=[]
		for m in pcigar.finditer(cigar):
			cigars.append((m.groups()[0],m.groups()[1]))
		read_start=0
		read_end=0
		Seq_pos=[]
		Read_Pos=[]
		Map_Pos=[]
		rec_seq=0
		rec_map=int(pin[3])
		for x in cigars:
			if x[1]=='I':
				rec_seq1=rec_seq
				rec_seq=rec_seq1+int(x[0])
				Seq_pos.append(pin[9][rec_seq1:rec_seq])
				Read_Pos.append([rec_seq1,rec_seq])
				Map_Pos.append([rec_map,rec_map])
			if x[1]=='D' or x[1]=='N':
				Seq_pos.append('')
				Read_Pos.append([rec_seq,rec_seq])
				rec_map1=rec_map
				rec_map=rec_map1+int(x[0])
				Map_Pos.append([rec_map1,rec_map])
			if x[1]=='M':
				rec_seq1=rec_seq
				rec_seq=rec_seq1+int(x[0])
				rec_map1=rec_map
				rec_map=rec_map1+int(x[0])
				Seq_pos.append(pin[9][rec_seq1:rec_seq])
				Read_Pos.append([rec_seq1,rec_seq])
				Map_Pos.append([rec_map1,rec_map])
			if x[1]=='S':
				rec_seq1=rec_seq
				rec_seq=rec_seq1+int(x[0])
				rec_map1=rec_map
				rec_map=rec_map1
				Seq_pos.append(pin[9][rec_seq1:rec_seq])
				Read_Pos.append([rec_seq1,rec_seq])
				Map_Pos.append([rec_map1,rec_map])
		Cigar_group=[]
		for x in cigars:
			if int(x[0])>min_length and not x[1]=='M':
				Cigar_group.append([x])
				Cigar_group.append([])
			else:
				if Cigar_group==[]:
					Cigar_group.append([x])
				else:
					Cigar_group[-1]+=[x]
		Seq2_pos=[]
		Read2_Pos=[]
		Map2_Pos=[]
		rec=-1
		for x in Cigar_group:
			Seq2_pos.append([])
			Read2_Pos.append([])
			Map2_Pos.append([])
			for y in x:
				rec+=1
				Seq2_pos[-1].append(Seq_pos[rec])
				Read2_Pos[-1].append(Read_Pos[rec])
				Map2_Pos[-1].append(Map_Pos[rec])
		chopped_reads=[''.join(x) for x in Seq2_pos if not x==[]]
		chopped_Read_Pos=[[x[0][0],x[-1][1]] for x in Read2_Pos if not x==[]]
		chopped_Map_Pos=[[x[0][0],x[-1][1]] for x in Map2_Pos if not x==[]]
		return [chopped_reads,chopped_Read_Pos,chopped_Map_Pos,[x for x in Cigar_group if not x==[]]]
	def path_mkdir(path):
			if not os.path.isdir(path):
					os.system(r'''mkdir %s'''%(path))
	def path_modify(path):
		if not path[-1]=='/':
			path+='/'
		return path
	def readSeq(filename):
		"""reads in a FASTA sequence"""
		stream = open(filename)
		seq = []
		for line in stream:
			if line.startswith(">"):
				continue
			seq.append(line.rstrip())
		return "".join(seq)
	def quality(hits):
		"""determines the quality of a list of hits"""
		slope1 = 1.0e6 / (825000 - 48000)
		slope2 = 1.0e6 / (914000 - 141000)
		offset1 = 0 - slope1*48000
		offset2 = 0 - slope2*141000
		goodhits = []
		for hit in hits:
			upper = slope1 * hit[0] + offset1
			lower = slope2 * hit[0] + offset2
			if lower < hit[1] < upper:
				goodhits.append(hit)
		return goodhits
	def makeDotplot(filename, hits, lenSeq1, lenSeq2):
		"""generate a dotplot from a list of hits
		   filename may end in the following file extensions:
			 *.ps, *.png, *.jpg
		"""
		x, y = zip(* hits)
		slope1 = 1.0e6 / (825000 - 48000)
		slope2 = 1.0e6 / (914000 - 141000)
		offset1 = 0 - slope1*48000
		offset2 = 0 - slope2*141000
		hits2 = quality(hits)
		print "%.5f%% hits on diagonal" % (100 * len(hits2) / float(len(hits)))
		# create plot
		p = plotting.Gnuplot()
		p.enableOutput(False)
		p.plot(x, y, xlab="sequence 2", ylab="sequence 1")
		p.plotfunc(lambda x: slope1 * x + offset1, 0, 1e6, 1e5)
		p.plotfunc(lambda x: slope2 * x + offset2, 0, 1e6, 1e5)
		# set plot labels
		p.set(xmin=0, xmax=lenSeq2, ymin=0, ymax=lenSeq1)
		p.set(main="dotplot (%d hits, %.5f%% hits on diagonal)" %
			  (len(hits), 100 * len(hits2) / float(len(hits))))
		p.enableOutput(True)
		# output plot
		p.save(filename)
		return p
	invert_base = { 'A' : 'T', 'T' : 'A', 'C' : 'G', 'G' : 'C','N' : 'N'}
	def subkeys(key, nth_base, inversions):
		subkeys = []
		keylen = len(key)
		# speed tip from:
		# http://wiki.python.org/moin/PythonSpeed/PerformanceTips#String_Concatenation
		if nth_base == 1:
			subkeys = [key]
		elif nth_base != 0:
			for k in range(nth_base):
				substr_list = [key[j] for j in range(keylen) if (j % nth_base == k)]
				subkeys.append("".join(substr_list))
		else:
			# nth_base = 0 is a special case for third base mismatches
			# for every codon, only include the first 2 bases in the hash
			subkeys = ["".join([key[i] for i in range(len(key)) if i % 3 != 2])]
		if inversions:
			for i in range(len(subkeys)):
				subkeys.append("".join([invert_base[c] for c in reversed(subkeys[i])]))
		return subkeys
	def kmerhits(seq1, seq2, kmerlen, nth_base=1, inversions=False):
		# hash table for finding hits
		lookup = {}
		# store sequence hashes in hash table
		#print "hashing seq1..."
		seq1len = len(seq1)
		for i in xrange(seq1len - kmerlen + 1):
			key = seq1[i:i+kmerlen]
			for subkey in subkeys(key, nth_base, inversions):
				lookup.setdefault(subkey, []).append(i)
		# match every nth base by 
		# look up hashes in hash table
		#print "hashing seq2..."
		hits = []
		for i in xrange(len(seq2) - kmerlen + 1):
			key = seq2[i:i+kmerlen]
			# only need to specify inversions for one seq
			for subkey in subkeys(key, nth_base, False):
				subhits = lookup.get(subkey, [])
				if subhits != []:
					# store hits to hits list
					for hit in subhits:
						hits.append((i, hit))
					# break out of loop to avoid doubly counting
					# exact matches
					break
		return hits
	def dotdata(kmerlen,seq1, seq2):
		# parse command-line arguments
		if len(sys.argv) < 4:
			print "you must call program as:  "
			print "   python ps1-dotplot.py <KMER LEN> <FASTA 1> <FASTA 2> <PLOT FILE>"
			print "   PLOT FILE may be *.ps, *.png, *.jpg"
			sys.exit(1)
		# match every nth base. Or set to 0 to allow third base mismatches
		nth_base = 1
		inversions = True
		hits = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
		return hits
	def dotplot(kmerlen,seq1,seq2,plotfile):
		nth_base = 1
		inversions = True
		hits = kmerhits(seq1, seq2, kmerlen, nth_base, inversions)
		p = makeDotplot(plotfile, hits, len(seq1), len(seq2))
	def write_dotfile(dotdata_for_record,txt_file,rec_start):
		fout='.'.join(txt_file.split('.')[:-1])+'.ref.'+rec_name+'.start.'+str(rec_start)
		fo=open(fout,'w')
		for x in dotdata_for_record[0]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
		fout='.'.join(txt_file.split('.')[:-1])+'.alt1.'+rec_name+'.start.'+str(rec_start)
		fo=open(fout,'w')
		for x in dotdata_for_record[1]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
		fout='.'.join(txt_file.split('.')[:-1])+'.alt2.'+rec_name+'.start.'+str(rec_start)
		fo=open(fout,'w')
		for x in dotdata_for_record[2]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
	def write_dotfile_left(dotdata_for_record,txt_file):
		fout='.'.join(txt_file.split('.')[:-1])+'.ref.left.'+rec_name+'.start.'+str(0)
		fo=open(fout,'w')
		for x in dotdata_for_record[0]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
		fout='.'.join(txt_file.split('.')[:-1])+'.alt1.left.'+rec_name+'.start.'+str(0)
		fo=open(fout,'w')
		for x in dotdata_for_record[1]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
		fout='.'.join(txt_file.split('.')[:-1])+'.alt2.left.'+rec_name+'.start.'+str(0)
		fo=open(fout,'w')
		for x in dotdata_for_record[2]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
	def write_dotfile_right(dotdata_for_record,txt_file):
		fout='.'.join(txt_file.split('.')[:-1])+'.ref.right.'+rec_name+'.start.'+str(0)
		fo=open(fout,'w')
		for x in dotdata_for_record[0]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
		fout='.'.join(txt_file.split('.')[:-1])+'.alt1.right.'+rec_name+'.start.'+str(0)
		fo=open(fout,'w')
		for x in dotdata_for_record[1]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
		fout='.'.join(txt_file.split('.')[:-1])+'.alt2.right.'+rec_name+'.start.'+str(0)
		fo=open(fout,'w')
		for x in dotdata_for_record[2]:
			print >>fo, ' '.join([str(i) for i in x])
		fo.close()
	def write_pacval_individual_stat_simple_short(txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2):
		fo=open('.'.join(txt_file.split('.')[:-1])+'.pacval','w')
		print >>fo, ' '.join(['alt','ref'])
		if len(rsquare_ref)>0:
			"""pick min Dis score from two alleles"""
			rsquare_alt_min=[min([rsquare_alt1[x],rsquare_alt2[x]]) for x in range(len(rsquare_ref))]
			"""if denominator=0, assign a small number"""
			for x in range(len(rsquare_ref)):
				if rsquare_alt_min[x]==0:
					rsquare_alt_min[x]=0.0001*rsquare_ref[x]
			""" modify stat from Dis to Dis(PacBio Read| Orignal Ref)/Dis(PacBio Read| Altered Ref)-1"""
			rsquare_out_ref=[0 for i in rsquare_ref]
			rsquare_out_alt_min=[float(rsquare_ref[x])/float(rsquare_alt_min[x])-1 for x in range(len(rsquare_alt_min))]
			for x in range(len(rsquare_out_ref)):
				print >>fo, ' '.join([str(i) for i in [rsquare_out_alt_min,rsquare_out_ref]])
		fo.close()
	def write_pacval_individual_stat_simple_long(txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2):
		fo=open('.'.join(txt_file.split('.')[:-1])+'.pacval','w')
		print >>fo, ' '.join(['alt','ref'])
		if len(rsquare_ref)>0:
			fo=open('.'.join(txt_file.split('.')[:-1])+'.pacval','w')
			"""pick Max Counts score from two counts"""
			rsquare_alt_max=[max([rsquare_alt1[x],rsquare_alt2[x]]) for x in range(len(rsquare_ref))]
			"""if denominator=0, assign a small number"""
			for x in range(len(rsquare_ref)):
				if rsquare_ref[x]==0:
					rsquare_ref[x]=0.0001*rsquare_alt_max[x]
			""" modify stat from Counts to 1-Counts(PacBio Read| Orignal Ref)/Dis(PacBio Read| Altered Ref)"""
			rsquare_out_ref=[0 for i in rsquare_ref]
			rsquare_out_alt_max=[1-float(rsquare_ref[x])/float(rsquare_alt_max[x]) for x in range(len(rsquare_alt_max))]
			for x in range(len(rsquare_out_ref)):
				print >>fo, ' '.join([str(i) for i in [rsquare_out_alt_max,rsquare_out_ref]])
		fo.close()
	def write_null_individual_stat(txt_file):
		fout='.'.join(txt_file.split('.')[:-1])+'.pacval'
		fo=open(fout,'w')
		print >>fo, ' '.join(['alt','ref'])
		fo.close()
	def name_hash_modify(name_hash):
		out=[]
		for x in name_hash:
			if '/' in x:
				y=x.replace('/','-')
			else:
				y=x
			out.append(y)
		return out
	if function_name=='simple':
		import os
		import sys
		import getopt
		import plotting
		import re
		import pickle
		import time
		import datetime
		import random
		import numpy
		import glob
		import numpy as np
		import scipy
		from scipy import stats
		import plotting
		def modify_read(pin):
			start=int(pin[3])
			real_start=100
			real_end=int(bps[-1])-int(bps[1])+real_start
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigar=pin[5]
			map_start=int(pin[3])
			pos_start=0
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			pos_rec=[]
			for x in cigars:
				if x[1]=='M':
					map_start+=int(x[0])
					pos_start+=int(x[0])
				if x[1]=='I':
					pos_start+=int(x[0])
				if x[1]=='D':
					map_start+=int(x[0])
				if x[1]=='S':
					pos_start+=int(x[0])
				if map_start>real_start:
					if pos_rec==[]:
						pos_rec.append(pos_start)
				if map_start>real_end:
					if len(pos_rec)==1:
						pos_rec.append(pos_start)
			return pin[9][pos_rec[0]:]
		def Draw_dotplot(sam_in):
			fin=open(sam_in)
			test=0
			tread=''
			start=300
			pin_rec=''
			for line in fin:
				pin=line.strip().split()
				if not pin[0][0]=='@':
					if len(pin[9])>test:
						if int(pin[3])<start:
							tread=pin[9]
							test=len(pin[9])
							pin_rec=pin
			fin.close()
			fref2=open(sam_in.replace('.sam','.dotplot.2.fa'),'w')
			print >>fref2, '>'+sam_in.split('/')[-1].replace('.sam','')
			print >>fref2,modify_read(pin_rec)
			fref2.close()
			#os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,bps[0],int(bps[1]),int(bps[-1]),sam_in.replace('.sam','.dotplot.1.fa')))
			dotmatcher = "/home/remills/data/local/bin/dotmatcher"
			os.system(r'''%s -asequence %s -bsequence %s -graph svg -threshold %s'''%(dotmatcher,sam_in.replace('.sam','.dotplot.2.fa'),sam_in.replace('.sam','.fa'),threshold))
			os.system(r'''rsvg-convert -f pdf -o dotmatcher.pdf dotmatcher.svg''')
			os.system(r'''mv %s %s'''%('./dotmatcher.svg',sam_in.replace('.sam','.svg')))
			os.system(r'''mv %s %s'''%('./dotmatcher.pdf',sam_in.replace('.sam','.pdf')))
			os.system(r'''%s -asequence %s -bsequence %s -graph svg -threshold %s'''%(dotmatcher,sam_in.replace('.sam','.dotplot.2.fa'),txt_in.replace('.txt','.ref.fa'),threshold))
			os.system(r'''rsvg-convert -f pdf -o dotmatcher.pdf dotmatcher.svg''')
			os.system(r'''mv %s %s'''%('./dotmatcher.svg',sam_in.replace('.sam','.2.svg')))
			os.system(r'''mv %s %s'''%('./dotmatcher.pdf',sam_in.replace('.sam','.2.pdf')))
		def cigar2reaadlength(cigar):
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			MapLen=0
			for n in cigars:
				if n[1]=='M' or n[1]=='D' or n[1]=='N':
					MapLen+=int(n[0])
			return MapLen
		def eu_dis_calcu_check_dup(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info):
			structure1=info[0].split('/')[0]
			structure2=info[1].split('/')[0]
			structure3=info[1].split('/')[1]
			if sum([int(dup_decide(structure1)),int(dup_decide(structure2)),int(dup_decide(structure3))])>0:
					eu_dis_calcu_1(fo_ref,rsquare_ref,y2,delta)
					eu_dis_calcu_1(fo1_alt,rsquare_alt1,y2,delta)
					eu_dis_calcu_1(fo2_alt,rsquare_alt2,y2,delta)
			else:
					eu_dis_calcu_2(fo_ref,rsquare_ref,y2,delta)
					eu_dis_calcu_2(fo1_alt,rsquare_alt1,y2,delta)
					eu_dis_calcu_2(fo2_alt,rsquare_alt2,y2,delta)
		def bps_check(bps):
			flag=0
			for x in bps[1:]:
				if x in chromosomes:
					return 1
			for x in range(len(bps)-2):
				if int(bps[x+2])-int(bps[x+1])>10**6:
					flag+=1
			return flag
		def cigar2alignstart(cigar,start,bps,flank_length):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[1])-flank_length: break
			return [read_rec,int(align_rec)-int(bps[1])+flank_length]
		def cigar2alignstart_2(cigar,start,bps):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[1]): break
			return [read_rec,int(align_rec)-int(bps[1])]
		def Capitalize_ref(fa2):
			fin=open(fa2)
			data=[]
			for line in fin:
				pin=line.strip().split()
				data+=pin
			fin.close()
			for x in range(1,len(data)):
				data[x]=data[x].upper()
			fin=open(fa2,'w')
			for x in data:
				print >>fin,x
			fin.close()
		def Capitalize_ref_forward_reverse(fa2):
			return [fa2.upper(),]
			fin=open(fa2)
			data=[]
			for line in fin:
				pin=line.strip().split()
				data+=pin
			fin.close()
			for x in range(1,len(data)):
				data[x]=data[x].upper()
			fin=open(fa2,'w')
			for x in data:
				print >>fin,x
			fin.close()
			data2=''.join(data[1:])[::-1]
			fin=open('.'.join(fa2.split('.')[:-1]+['reverse.fa']),'w')
			print >>fin, data[0]
			for x in chop_x_for_fasta(data2):
				print >>fin,x
			fin.close()
		def Capitalize_all(list):
			for x in list:
				Capitalize_ref(x)
		def calcu_eu_dis_long(SVelter_Name,info,new_structure):
			flank_length=500
			txt_rec_file_writing(out_path,SVelter_Name,info[2:])
			txt_file=out_path+SVelter_Name+'.txt'
			out=[[],[],[]]
			out_al_a=[[],[],[]]
			out_al_b=[[],[],[]]
			for x_structure in new_structure:
				fa_file_names=Pacbio_prodce_ref_alt_long(ref,x_structure,txt_file)
				ref_info1=[x_structure[2][0],x_structure[2][1],x_structure[2][1]+2*flank_length]
				ref_info2=[x_structure[3][0],x_structure[3][1],x_structure[3][1]+2*flank_length]
				delta=50
				rsquare_alt=[]
				rsquare_ref1=[]
				rsquare_ref2=[]
				if bps_check(ref_info1)==0:
					bps=ref_info1
					end_point=bps[-1]
					all_reads1=chop_pacbio_read_long(bps)
				if bps_check(ref_info2)==0:
					bps=ref_info2
					end_point=bps[-1]
					all_reads2=chop_pacbio_read_long(bps)
				read_hash=all_reads1[0]+all_reads2[0]
				miss_hash=all_reads1[1]+all_reads2[1]
				name_hash=name_hash_modify(all_reads1[2]+all_reads2[2])
				rec_len=0
				fa1=out_path+'temp.fa'
				Capitalize_ref(fa_file_names[0])
				Capitalize_ref(fa_file_names[1])
				Capitalize_ref(fa_file_names[2])
				miss_rec=-1
				if not read_hash==[]:
					if len(read_hash)>20:
						new_read_hash=random.sample(range(len(read_hash)),20)
						read_hash2=[read_hash[i] for i in new_read_hash]
						miss_hash2=[miss_hash[i] for i in new_read_hash]
						name_hash2=[name_hash[i] for i in new_read_hash]
						read_hash=read_hash2
						miss_hash=miss_hash2
						name_hash=name_hash2
					name_rec=['0','0']
					fo_alt=txt_file.replace('.txt','.dotplot.alt')
					fo1_ref=txt_file.replace('.txt','.dotplot.ref1')
					fo2_ref=txt_file.replace('.txt','.dotplot.ref2')
					fa2=txt_file.replace('.txt','.alt.fa')
					fa3=txt_file.replace('.txt','.ref1.fa')
					fa4=txt_file.replace('.txt','.ref2.fa')
					fa5=txt_file.replace('.txt','.ref.fa')
					for x in read_hash:
						if len(x)>1000:
							x=x[:1000]
						miss_rec+=1
						y=x
						y2=miss_hash[miss_rec]
						y3=name_hash[miss_rec]
						fo=open(out_path+'temp.fa','w')
						print >>fo, '>temp'
						for z in chop_x_for_fasta(y):
								print >>fo, z
						fo.close()
						os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,fa2,fo_alt))
						os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,fa3,fo1_ref))
						os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,fa4,fo2_ref))
						eu_dis_calcu_bothend(fo_alt,fo1_ref,fo2_ref,rsquare_alt,rsquare_ref1,rsquare_ref2,y2,delta,info)
						if not len(rsquare_alt)==len(rsquare_ref1)==len(rsquare_ref2):
								min_len=min([len(rsquare_alt),len(rsquare_ref1),len(rsquare_ref2)])
								rsquare_alt=rsquare_alt[:min_len]
								rsquare_ref1=rsquare_ref1[:min_len]
								rsquare_ref2=rsquare_ref2[:min_len]
						if not rsquare_alt==[]:
							if rsquare_alt[-1]-max([abs(rsquare_ref1[-1]),abs(rsquare_ref2[-1])])<rec_len:
								rec_len=rsquare_alt[-1]-max([abs(rsquare_ref1[-1]),abs(rsquare_ref2[-1])])
								rec_start=y2
								rec_name=y3
								os.system(r'''cp %s %s'''%(fo_alt+'.txt',fo_alt+'.'+y3))
								os.system(r'''cp %s %s'''%(fo1_ref+'.txt',fo1_ref+'.'+y3))
								os.system(r'''cp %s %s'''%(fo2_ref+'.txt',fo2_ref+'.'+y3))
								os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.fa')))	 
								name_rec[0]=y3
					if	rec_len==0:
						os.system(r'''cp %s %s'''%(fo_alt+'.txt',fo_alt+'.'+y3))
						os.system(r'''cp %s %s'''%(fo1_ref+'.txt',fo1_ref+'.'+y3))
						os.system(r'''cp %s %s'''%(fo2_ref+'.txt',fo2_ref+'.'+y3))
						rec_start=y2
						rec_name=y3
						os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.fa')))
					test_fin=os.popen(r'''wc -l %s'''%(fo_alt+'.'+rec_name))
					test_pin=test_fin.readline().strip().split()
					test_fin.close()
					if not test_pin[0]=='0':
						os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),fa2,fo_alt+rec_name+'.png'))
						os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),fa3,fo1_ref+rec_name+'.png'))
						os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),fa4,fo2_ref+rec_name+'.png'))
						os.system(r'''mv %s %s'''%(fo_alt+'.'+rec_name,fo_alt+'.'+rec_name+'.start.'+str(rec_start)))
						os.system(r'''mv %s %s'''%(fo1_ref+'.'+rec_name,fo1_ref+'.'+rec_name+'.start.'+str(rec_start)))
						os.system(r'''mv %s %s'''%(fo2_ref+'.'+rec_name,fo2_ref+'.'+rec_name+'.start.'+str(rec_start)))
				remove_files_short(txt_file)
				out1=[[1 for i in rsquare_alt],[abs(float(rsquare_ref1[i])/float(rsquare_alt[i])) for i in range(len(rsquare_alt)) if not rsquare_alt[i]==0],[abs(float(rsquare_ref2[i])/float(rsquare_alt[i])) for i in range(len(rsquare_alt)) if not rsquare_alt[i]==0]]
				out_al1=[[],[],[]]
				out_al2=[[],[],[]]
				for x in range(len(out1[2])):
					if out1[1][x]>out1[2][x]:
						out_al1[0].append(out1[0][x])
						out_al1[1].append(out1[1][x])
						out_al1[2].append(out1[2][x])
					else:
						out_al2[0].append(out1[0][x])
						out_al2[1].append(out1[1][x])
						out_al2[2].append(out1[2][x])
				out[0] +=out1[0]
				out_al_a[2]+=out_al1[2]
				out_al_b[1]+=out_al2[1]
			return [out_al_a[2]+out_al_b[1],out[0]]
		def delta_calculate(bps):
			if (int(bps[2])-int(bps[1]))>1000:
				delta=(int(bps[2])-int(bps[1]))/10
			else:
				delta=50
			return delta
		def flank_length_calculate(bps):
			if int(bps[-1])-int(bps[1])<100:
				flank_length=2*(int(bps[-1])-int(bps[1]))
			else:
				if int(bps[-1])-int(bps[1])<500:
					flank_length=int(bps[-1])-int(bps[1])
				else:
					flank_length=500
			return flank_length
		def bl_len_hash_calculate(bps,ref_sv):
			bl_len_hash={}
			for x in ref_sv.split('/')[0]:
				bl_len_hash[x]=int(bps[ord(x)-97+2])-int(bps[ord(x)-97+1])
			return bl_len_hash
		def end_point_calculate(alt_sv,bl_len_hash):
			end_point=0
			for x in alt_sv.split('/')[0]:
				if not x=='^':
					end_point+=bl_len_hash[x]
			return end_point
		def read_hash_minimize(all_reads):
			read_hash=all_reads[0]
			miss_hash=all_reads[1]
			name_hash=name_hash_modify(all_reads[2])							
			if len(read_hash)>20:
				new_read_hash=random.sample(range(len(read_hash)),20)
				read_hash2=[read_hash[i] for i in new_read_hash]
				miss_hash2=[miss_hash[i] for i in new_read_hash]
				name_hash2=[name_hash[i] for i in new_read_hash]
				read_hash=read_hash2
				miss_hash=miss_hash2
				name_hash=name_hash2
			return [read_hash,miss_hash,name_hash]
		def chop_pacbio_read_simple_short(flank_length,info):
			bps=info[2:]
			block_length={}
			for x in range(len(info[2:])-2):
				block_length[chr(97+x)]=int(info[x+4])-int(info[x+3])
			alA_len=numpy.sum([block_length[x] for x in info[1].split('/')[0] if not x=='^'])
			alB_len=numpy.sum([block_length[x] for x in info[1].split('/')[1] if not x=='^'])
			alRef_len=int(info[-1])-int(info[3])
			len_cff=max([alA_len,alB_len,alRef_len])
			if len_cff>10**4:
				len_cff=min([alA_len,alB_len,alRef_len])
			fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
			out=[]
			out2=[]
			out3=[]
			for line in fbam:
				pbam=line.strip().split()
				#out3.append(int(pbam[3])-int(bps[1])+flank_length)
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
						align_start=align_info[0]
						miss_bp=align_info[1]
						align_pos=int(pbam[3])
						target_read=pbam[9][align_start:]
						if len(target_read)>flank_length+len_cff:
							out.append(target_read[:max([alA_len,alB_len,alRef_len])+2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
							#out3.append(pbam[:9])
					elif int(pbam[3])<int(bps[1])+1:
						align_info=cigar2alignstart_2(pbam[5],int(pbam[3]),bps)
						align_start=align_info[0]
						miss_bp=align_info[1]+flank_length
						target_read=pbam[9][align_start:]
						if len(target_read)>len_cff:
							out.append(target_read[:max([alA_len,alB_len,alRef_len])+flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
							#out3.append(pbam[:9])
			fbam.close()
			return [out,out2,out3]
		def calcu_eu_dis_simple_short(out_path,sample_name,info):
			#eg of info:['a/a', 'a^/a^', '8', '1412203', '1412359']
			k1=info[0]
			k2=info[1]
			k3=info[2:]
			case_name='_'.join(k3)
			txt_file=out_path+sample_name+'.'+case_name+'.txt'
			fo=open(txt_file,'w')
			print >>fo, ' '.join([str(i) for i in [k1,k2]+k3])
			fo.close()
			temp_bam=txt_file.replace('.txt','.bam')
			temp_sam=txt_file.replace('.txt','.sam')
			ref_sv=info[0]
			alt_sv=info[1]
			chrom=info[2]
			bps=info[2:]
			delta=delta_calculate(bps)
			flank_length=flank_length_calculate(bps)
			if bps_check(bps)==0:
				bl_len_hash=bl_len_hash_calculate(bps,ref_sv)
				end_point=end_point_calculate(alt_sv,bl_len_hash)
				all_reads=chop_pacbio_read_simple_short(flank_length,info)
				if not all_reads[0]==[]:
					all_reads_new=read_hash_minimize(all_reads)
					read_hash=all_reads_new[0]
					miss_hash=all_reads_new[1]
					name_hash=name_hash_modify(all_reads_new[2])
					rsquare_ref=[]
					rsquare_alt1=[]
					rsquare_alt2=[]
					rec_len=0
					rec_start=0
					rec_name='0'
					if not read_hash==[]:
						seqs=Pacbio_prodce_ref_alt_short(ref,flank_length,txt_file)
						seq2=seqs[0]
						seq3=seqs[1]
						seq4=seqs[2]
						miss_rec=-1
						dotdata_for_record=[]
						seq1=''
						for x in read_hash:
							miss_rec+=1
							y=x
							y2=miss_hash[miss_rec]
							y3=name_hash[miss_rec]
							dotdata_ref=dotdata(window_size,y,seq2)
							dotdata_alt1=dotdata(window_size,y,seq3)
							dotdata_alt2=dotdata(window_size,y,seq4)							
							eu_dis_calcu_check_dup(dotdata_ref,dotdata_alt1,dotdata_alt2,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
							if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
								min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
								rsquare_ref=rsquare_ref[:min_len]
								rsquare_alt1=rsquare_alt1[:min_len]
								rsquare_alt2=rsquare_alt2[:min_len]
							if not rsquare_alt1==[]:
								if min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]<rec_len:
									rec_len=min([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
									rec_start=y2
									rec_name=y3
									seq1=y
									dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
						if	rec_len==0:
							dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
							rec_name=y3
							rec_start=y2
							seq1=y
						if len(dotdata_for_record[0])>0:
							dotplot(window_size,seq1,seq2,out_path+sample_name+'.'+case_name+'.ref.'+rec_name+'.png')
							dotplot(window_size,seq1,seq3,out_path+sample_name+'.'+case_name+'.alt1.'+rec_name+'.png')
							dotplot(window_size,seq1,seq4,out_path+sample_name+'.'+case_name+'.alt2.'+rec_name+'.png')
							write_dotfile(dotdata_for_record,txt_file,rec_start)
							write_pacval_individual_stat_simple_short(txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2)
				else:
					write_null_individual_stat(txt_file)
					remove_files_short(txt_file)
		def calcu_eu_dis_simple_long(out_path,sample_name,info):
			k1=info[0]
			k2=info[1]
			k3=info[2:]
			case_name='_'.join(k3)
			fo=open(out_path+sample_name+'.'+case_name+'.txt','w')
			print >>fo, ' '.join([str(i) for i in [k1,k2]+k3])
			fo.close()
			txt_file=out_path+sample_name+'.'+case_name+'.txt'
			temp_bam=txt_file.replace('.txt','.bam')
			temp_sam=txt_file.replace('.txt','.sam')
			ref_sv=info[0]
			alt_sv=info[1]
			chrom=info[2]
			bps=info[2:]
			delta=delta_calculate(bps)
			if bps_check(bps)==0:
				bl_len_hash={}
				for x in ref_sv.split('/')[0]:
					bl_len_hash[x]=int(bps[ord(x)-97+2])-int(bps[ord(x)-97+1])
				end_point=0
				for x in alt_sv.split('/')[0]:
					if not x=='^':
						end_point+=bl_len_hash[x]
				flank_length=500
				all_reads_left=chop_pacbio_read_left(bps,flank_length)
				all_reads_right=chop_pacbio_read_right(bps,flank_length)
				read_hash_left=all_reads_left[0]
				miss_hash_left=all_reads_left[1]
				name_hash_left=name_hash_modify(all_reads_left[2])
				read_hash_right=all_reads_right[0]
				miss_hash_right=all_reads_right[1]
				name_hash_right=name_hash_modify(all_reads_right[2])
				rsquare_ref=[]
				rsquare_alt1=[]
				rsquare_alt2=[]
				seqs=Pacbio_prodce_ref_alt_short(ref,flank_length,txt_file)
				seq2=seqs[0]
				seq3=seqs[1]
				seq4=seqs[2]
				if not read_hash_left==[]:
					miss_rec=-1
					rec_len=0
					if len(read_hash_left)>20:
						read_num_sample=random.sample(range(len(read_hash_left)),20)
						new_read_hash_left=[read_hash_left[i] for i in read_num_sample]
						new_miss_hash_left=[miss_hash_left[i] for i in read_num_sample]
						new_name_hash_left=[name_hash_left[i] for i in read_num_sample]
						read_hash_left=new_read_hash_left
						miss_hash_left=new_miss_hash_left
						name_hash_left=new_name_hash_left
					for x in read_hash_left:
						miss_rec+=1
						y=x
						y3=name_hash_left[miss_rec]
						dotdata_ref=dotdata(window_size,y,seq2)
						dotdata_alt1=dotdata(window_size,y,seq3)
						dotdata_alt2=dotdata(window_size,y,seq4)							
						temp_rsquare=eu_dis_calcu_simple_long(dotdata_ref,delta)
						if not temp_rsquare==0:
							rsquare_ref.append(temp_rsquare)
						temp_rsquare=eu_dis_calcu_simple_long(dotdata_alt1,delta)
						if not temp_rsquare==0:
							rsquare_alt1.append(temp_rsquare)
						temp_rsquare=eu_dis_calcu_simple_long(dotdata_alt2,delta)
						if not temp_rsquare==0:
							rsquare_alt2.append(temp_rsquare)
						if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
							min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
							rsquare_ref=rsquare_ref[:min_len]
							rsquare_alt1=rsquare_alt1[:min_len]
							rsquare_alt2=rsquare_alt2[:min_len]
						if not rsquare_alt1==[]:
							if max([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]>rec_len:
								rec_len=max([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
								rec_name=y3
								dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
					if rec_len>0:
						if len(dotdata_for_record[0])>0:
							dotplot(window_size,seq1,seq2,out_path+sample_name+'.'+case_name+'.ref.left.'+rec_name+'.png')
							dotplot(window_size,seq1,seq3,out_path+sample_name+'.'+case_name+'.alt1.left.'+rec_name+'.png')
							dotplot(window_size,seq1,seq4,out_path+sample_name+'.'+case_name+'.alt2.left.'+rec_name+'.png')
							write_dotfile_left(dotdata_for_record,txt_file)
				if not read_hash_right==[]:
					miss_rec=-1
					rec_len=0
					if len(read_hash_right)>20:
						read_num_sample=random.sample(range(len(read_hash_right)),20)
						new_read_hash_right=[read_hash_right[i] for i in read_num_sample]
						new_miss_hash_right=[miss_hash_right[i] for i in read_num_sample]
						new_name_hash_right=[name_hash_right[i] for i in read_num_sample]
						read_hash_right=new_read_hash_right
						miss_hash_right=new_miss_hash_right
						name_hash_right=new_name_hash_right
					for x in read_hash_right:
						miss_rec+=1
						y=x
						y3=name_hash_right[miss_rec]
						dotdata_ref=dotdata(window_size,y,seq2[::-1])
						dotdata_alt1=dotdata(window_size,y,seq3[::-1])
						dotdata_alt2=dotdata(window_size,y,seq4[::-1])							
						temp_rsquare=eu_dis_calcu_simple_long(dotdata_ref,delta)
						if not temp_rsquare==0:
							rsquare_ref.append(temp_rsquare)
						temp_rsquare=eu_dis_calcu_simple_long(dotdata_alt1,delta)
						if not temp_rsquare==0:
							rsquare_alt1.append(temp_rsquare)
						temp_rsquare=eu_dis_calcu_simple_long(dotdata_alt2,delta)
						if not temp_rsquare==0:
							rsquare_alt2.append(temp_rsquare)
						if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
							min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
							rsquare_ref=rsquare_ref[:min_len]
							rsquare_alt1=rsquare_alt1[:min_len]
							rsquare_alt2=rsquare_alt2[:min_len]
						if not rsquare_alt1==[]:
							if max([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]>rec_len:
								rec_len=max([rsquare_alt1[-1],rsquare_alt2[-1]])-rsquare_ref[-1]
								rec_name=y3
								dotdata_for_record=[dotdata_ref,dotdata_alt1,dotdata_alt2]
					if rec_len>0:
						if not len(dotdata_for_record[0])>0:
							dotplot(window_size,seq1,seq2[::-1],out_path+sample_name+'.'+case_name+'.ref.right.'+rec_name+'.png')
							dotplot(window_size,seq1,seq3[::-1],out_path+sample_name+'.'+case_name+'.alt1.right.'+rec_name+'.png')
							dotplot(window_size,seq1,seq4[::-1],out_path+sample_name+'.'+case_name+'.alt2.right.'+rec_name+'.png')
							write_dotfile_right(dotdata_for_record,txt_file)
				write_pacval_individual_stat_simple_long(txt_file,rsquare_ref,rsquare_alt1,rsquare_alt2)
				remove_files_short(txt_file)
		def eu_dis_calcu_simple_long(fo_ref,delta):
			temp_data1=[[],[]]
			for x in fo_ref:
				temp_data1[0].append(int(x[0]))
				temp_data1[1].append(int(x[1]))
			temp_data2=[[],[]]
			for x in range(len(temp_data1[0])):
				if numpy.abs(temp_data1[0][x]-temp_data1[1][x])<temp_data1[0][x]/10:
					temp_data2[0].append(temp_data1[0][x])
					temp_data2[1].append(temp_data1[1][x])			
			return len(temp_data2[0])
		def cigar2alignstart_left(cigar,start,bps):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[1])-flank_length: break
			return [read_rec,int(align_rec)-int(bps[1])+flank_length]
		def cigar2alignstart_right(cigar,start,bps):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[2])+flank_length: break
			return [read_rec,-int(align_rec)+(int(bps[2])+flank_length)]
		def cigar2alignstart_2(cigar,start,bps):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[1]): break
			return [read_rec,int(align_rec)-int(bps[1])]
		def chop_pacbio_read_left(bps,flank_length):
			fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[1])+flank_length))
			out=[]
			out2=[]
			out3=[]
			for line in fbam:
				pbam=line.strip().split()
				if not pbam[0]=='@': 
					if int(pbam[3])<int(bps[1])-flank_length+1:
						align_info=cigar2alignstart_left(pbam[5],int(pbam[3]),bps)
						align_start=align_info[0]
						miss_bp=align_info[1]
						target_read=pbam[9][align_start:]
						if len(target_read)>2*flank_length:
							out.append(target_read[:2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
			fbam.close()
			return [out,out2,out3]
		def chop_pacbio_read_right(bps,flank_length):
			fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[2])-flank_length,int(bps[2])+flank_length))
			out=[]
			out2=[]
			out3=[]
			for line in fbam:
				pbam=line.strip().split()
				if not pbam[0]=='@': 
					if int(pbam[3])+int(cigar2reaadlength(pbam[5]))>int(bps[2])+flank_length-1:
						align_info=cigar2alignstart_right(pbam[5],int(pbam[3]),bps)
						align_end=align_info[0]
						miss_bp=align_info[1]
						target_read=pbam[9][:align_end]
						if len(target_read)>2*flank_length:
							out.append(target_read[::-1][:2*flank_length])
							out2.append(miss_bp)
							out3.append(pbam[0])
			fbam.close()
			return [out,out2,out3]
		def rsquare_file_analysis(rsquare_file):
			fin=open(rsquare_file)
			pin=fin.readline().strip().split()
			if not pin[0]=='alt' or not pin[1]=='ref':
				return 'Error'
			else:
				total_num=0
				pass_qc_mum=0
				for line in fin:
					pin=line.strip().split()
					total_num+=1
					if float(pin[0])>float(pin[1]):
						pass_qc_mum+=1
			fin.close()
			if not total_num==0:
				return float(pass_qc_mum)/float(total_num)
			else:
				return 'Nan'
		def txt_file_analysis(txt_file):
			fin=open(txt_file)
			pin=fin.readline().strip().split()
			fin.close()
			return pin
		def Pacbio_Vali_Simple_Prep():
			global out_prefix
			out_prefix=dict_opts['--output-prefix']
			global out_path
			if not '--output-path' in dict_opts.keys():
				out_path='/'.join(out_prefix.split('/')[:-1])+'/'
			else:
				out_path=path_modify(dict_opts['--output-path'])
			path_mkdir(out_path)
			filein=dict_opts['--vali-file']
			fin=open(filein)
			global sample_name
			sample_name=filein.split('/')[-1].split('.')[0]
			global case_hash
			case_hash={}
			for line in fin:
				pin=line.strip().split()
				if not pin==[]:
					if not pin[0] in case_hash.keys():
						case_hash[pin[0]]={}
					if not pin[1] in case_hash[pin[0]].keys():
						case_hash[pin[0]][pin[1]]=[]
					case_hash[pin[0]][pin[1]].append(pin[2:])
			fin.close()
			global window_size
			if '--window-size' in dict_opts.keys():
				window_size=int(dict_opts['--window-size'])
			else:
				window_size=10
			global start
			global delta
			start=0
			delta=50
			global bam_in
			global ref
			global chromosomes
			bam_in=dict_opts['--bam-file']
			ref=dict_opts['--reference']
			chromosomes=chromos_readin(ref)
		def Pacbio_Vali_Simple_parametre_readin():
			opts,args=getopt.getopt(sys.argv[2:],'i:',['vcf-file=','window-size=','bam-file=','reference=','vali-file=','output-prefix=','output-path='])
			global dict_opts
			dict_opts=dict(opts)
		Pacbio_Vali_Simple_parametre_readin()
		Pacbio_Vali_Simple_Prep()
		min_length=50
		min_read_compare=20
		for k1 in case_hash.keys():
			for k2 in case_hash[k1].keys():
				for k3 in case_hash[k1][k2]:
					current_info=[k1,k2]+k3
					if int(k3[2])-int(k3[1])<5000:
						calcu_eu_dis_simple_short(out_path,sample_name,[k1,k2]+k3)
					else:
						calcu_eu_dis_simple_long(out_path,sample_name,[k1,k2]+k3)
		fo=open(out_prefix+'.stat','w')
		print >>fo, ' '.join(['score','ref','alt','pos'])
		for k1 in os.listdir(out_path):
			if k1.split('.')[-1]=='rsquare':
				qc_score_structure=rsquare_file_analysis(out_path+k1)
				info=txt_file_analysis(out_path+'.'.join(k1.split('.')[:-1]+['txt']))
				print >>fo, ' '.join([str(i) for i in [qc_score_structure]+info[:2]+['_'.join(info[2:])]])
		fo.close()
	if function_name=='complex':
		import getopt
		import re
		import pickle
		import time
		import datetime
		import random
		import numpy
		import glob
		import numpy as np
		import scipy
		from scipy import stats
		def path_name_check(Path_Name):
			if not Path_Name[-1]=='/':
				Path_Name+='/'
			return Path_Name
		def csv_info_readin(SVelter_in):
			hash1={}
			fin=open(SVelter_in)
			for line in fin:
				pin=line.strip().split()
				if not pin==[]:
					if not pin[2] in hash1.keys():
						hash1[pin[2]]={}
					if not int(pin[3]) in hash1[pin[2]].keys():
						hash1[pin[2]][int(pin[3])]={}
					if not int(pin[-1]) in hash1[pin[2]][int(pin[3])].keys():
						hash1[pin[2]][int(pin[3])][int(pin[-1])]=[]
					if not pin in hash1[pin[2]][int(pin[3])][int(pin[-1])]:
						hash1[pin[2]][int(pin[3])][int(pin[-1])].append(pin)
			fin.close()
			hash2={}
			for k1 in hash1.keys():
				hash2[k1]=[]
				for k2 in sorted(hash1[k1].keys()):
					for k3 in sorted(hash1[k1][k2].keys()):
						for k4 in hash1[k1][k2][k3]:
							hash2[k1].append([k2,k3]+k4)
			return hash2
		def alt_SV_produce(k1):
			if 'DEL' in k1.split('.') or 'del' in k1.split('.') or 'Del' in k1.split('.'):
				return ''
			elif 'DUP' in k1.split('.') or 'dup' in k1.split('.') or 'Dup' in k1.split('.'):
				return 'aa'
			elif 'INV' in k1.split('.') or 'inv' in k1.split('.') or 'Inv' in k1.split('.'):
				return 'a^'
			else:
				return 'a'
		def Other_Algorithm_result_readin(Path_Delly):
			out_hash={}
			if not Path_Delly=='Null':
				ref='a/a'
				for k1 in os.listdir(Path_Delly):
					if k1.split('.')[-1]=='vcf':
						os.system(r'''SV.Simple.Output.Process.py vcf-to-bed --input %s'''%(Path_Delly+k1))
					if k1.split('.')[-1]=='bedpe':
						os.system(r'''SV.Simple.Output.Process.py bedpe-to-bed --input %s --reference %s'''%(Path_Delly+k1,dict_opts['--reference']))
				for k1 in os.listdir(Path_Delly):
					if k1.split('.')[-1]=='bed':
						fin=open(Path_Delly+k1)
						alt_allele=alt_SV_produce(k1)
						if not alt_allele=='a':
							for line in fin:
								pin=line.strip().split()
								if pin[-1]=='homo':
									alt='/'.join([alt_allele,alt_allele])
								else:
									alt='/'.join([alt_allele,'a'])
								if not pin[0] in out_hash.keys():
									out_hash[pin[0]]={}
								if not int(pin[1]) in out_hash[pin[0]].keys():
									out_hash[pin[0]][int(pin[1])]={}
								if not int(pin[2]) in out_hash[pin[0]][int(pin[1])].keys():
									out_hash[pin[0]][int(pin[1])][int(pin[2])]=[]
								out_hash[pin[0]][int(pin[1])][int(pin[2])].append([ref,alt])
						fin.close()
				out_hash2={}
				for k1 in out_hash.keys():
					out_hash2[k1]=[]
					for k2 in sorted(out_hash[k1].keys()):
						for k3 in sorted(out_hash[k1][k2].keys()):
							for k4 in out_hash[k1][k2][k3]:
								out_hash2[k1].append([k2,k3]+k4)
				return out_hash2
			else:
				return out_hash
		def compare_hash():
			out={}
			for k1 in SVelter_hash.keys():
				out[k1]=[]
				for k2 in SVelter_hash[k1]:
					temp=[[],[],[]]
					if k1 in Delly_hash.keys():
						for k3 in Delly_hash[k1]:
							if k3[1]<k2[0]: continue
							elif k3[0]>k2[1]: break
							else:
								if float(sorted(k2[:2]+k3[:2])[2]-sorted(k2[:2]+k3[:2])[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.1:
									if not k3 in temp[0]:
										temp[0].append(k3)
					if k1 in Lumpy_hash.keys():
						for k3 in Lumpy_hash[k1]:
							if k3[1]<k2[0]: continue
							elif k3[0]>k2[1]: break
							else:
								if float(sorted(k2[:2]+k3[:2])[2]-sorted(k2[:2]+k3[:2])[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.1:
									if not k3 in temp[1]:
										temp[1].append(k3)
					if k1 in Pindel_hash.keys():
						for k3 in Pindel_hash[k1]:
							if k3[1]<k2[0]: continue
							elif k3[0]>k2[1]: break
							else:
								if float(sorted(k2[:2]+k3[:2])[2]-sorted(k2[:2]+k3[:2])[1])/float(max([k2[1]-k2[0],k3[1]-k3[0]]))>0.1:
									if not k3 in temp[1]:
										temp[2].append(k3)
					#if not temp[0]==[] and not temp[1]==[]:
					out[k1].append([[k2]]+temp)
			return out
		def bps_check(bps):
			flag=0
			for x in range(len(bps)-2):
				if int(bps[x+2])-int(bps[x+1])>10**6:
					flag+=1
			return flag
		def flank_length_decide(bps):
			if int(bps[-1])-int(bps[1])<100:
				flank_length=2*(int(bps[-1])-int(bps[1]))
			else:
				if int(bps[-1])-int(bps[1])<500:
					flank_length=int(bps[-1])-int(bps[1])
				else:
					flank_length=500
			return flank_length
		def modify_read(pin):
			start=int(pin[3])
			real_start=100
			real_end=int(bps[-1])-int(bps[1])+real_start
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigar=pin[5]
			map_start=int(pin[3])
			pos_start=0
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			pos_rec=[]
			for x in cigars:
				if x[1]=='M':
					map_start+=int(x[0])
					pos_start+=int(x[0])
				if x[1]=='I':
					pos_start+=int(x[0])
				if x[1]=='D':
					map_start+=int(x[0])
				if x[1]=='S':
					pos_start+=int(x[0])
				if map_start>real_start:
					if pos_rec==[]:
						pos_rec.append(pos_start)
				if map_start>real_end:
					if len(pos_rec)==1:
						pos_rec.append(pos_start)
			return pin[9][pos_rec[0]:]
		def Draw_dotplot(sam_in):
			fin=open(sam_in)
			test=0
			tread=''
			start=300
			pin_rec=''
			for line in fin:
				pin=line.strip().split()
				if not pin[0][0]=='@':
					if len(pin[9])>test:
						if int(pin[3])<start:
							tread=pin[9]
							test=len(pin[9])
							pin_rec=pin
			fin.close()
			fref2=open(sam_in.replace('.sam','.dotplot.2.fa'),'w')
			print >>fref2, '>'+sam_in.split('/')[-1].replace('.sam','')
			print >>fref2,modify_read(pin_rec)
			fref2.close()
			#os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,bps[0],int(bps[1]),int(bps[-1]),sam_in.replace('.sam','.dotplot.1.fa')))
			dotmatcher = "/home/remills/data/local/bin/dotmatcher"
			os.system(r'''%s -asequence %s -bsequence %s -graph svg -threshold %s'''%(dotmatcher,sam_in.replace('.sam','.dotplot.2.fa'),sam_in.replace('.sam','.fa'),threshold))
			os.system(r'''rsvg-convert -f pdf -o dotmatcher.pdf dotmatcher.svg''')
			os.system(r'''mv %s %s'''%('./dotmatcher.svg',sam_in.replace('.sam','.svg')))
			os.system(r'''mv %s %s'''%('./dotmatcher.pdf',sam_in.replace('.sam','.pdf')))
			os.system(r'''%s -asequence %s -bsequence %s -graph svg -threshold %s'''%(dotmatcher,sam_in.replace('.sam','.dotplot.2.fa'),txt_in.replace('.txt','.ref.fa'),threshold))
			os.system(r'''rsvg-convert -f pdf -o dotmatcher.pdf dotmatcher.svg''')
			os.system(r'''mv %s %s'''%('./dotmatcher.svg',sam_in.replace('.sam','.2.svg')))
			os.system(r'''mv %s %s'''%('./dotmatcher.pdf',sam_in.replace('.sam','.2.pdf')))
		def cigar2reaadlength(cigar):
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			MapLen=0
			for n in cigars:
				if n[1]=='M' or n[1]=='D' or n[1]=='N':
					MapLen+=int(n[0])
			return MapLen
		def rsquare_calcu(fo_ref,rsquare_ref):
			fo1=open(fo_ref+'.txt')
			temp_data1=[[],[]]
			for line in fo1:
				po1=line.strip().split()
				temp_data1[0].append(int(po1[0]))
				temp_data1[1].append(int(po1[1]))
			fo1.close()
			slope, intercept, r_value, p_value, std_err = stats.linregress(temp_data1[0],temp_data1[1])
			rsquare_ref.append(r_value**2)
		def eu_dis_calcu(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info):
			structure1=info[2].split('/')[0]
			structure2=info[3].split('/')[0]
			structure3=info[3].split('/')[1]
			if sum([int(dup_decide(structure1)),int(dup_decide(structure2)),int(dup_decide(structure3))])>0:
					eu_dis_calcu_1(fo_ref,rsquare_ref,y2,delta)
					eu_dis_calcu_1(fo1_alt,rsquare_alt1,y2,delta)
					eu_dis_calcu_1(fo2_alt,rsquare_alt2,y2,delta)
			else:
					eu_dis_calcu_2(fo_ref,rsquare_ref,y2,delta)
					eu_dis_calcu_2(fo1_alt,rsquare_alt1,y2,delta)
					eu_dis_calcu_2(fo2_alt,rsquare_alt2,y2,delta)
		def eu_dis_calcu_complex_long(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info):
			structure1=info[2].split('/')[0]
			structure2=info[3].split('/')[0]
			structure3=info[3].split('/')[1]
			eu_dis_calcu_3(fo_ref,rsquare_ref,y2,delta)
			eu_dis_calcu_3(fo1_alt,rsquare_alt1,y2,delta)
			eu_dis_calcu_3(fo2_alt,rsquare_alt2,y2,delta)
		def remove_files(txt_file):
			os.system(r'''rm %s'''%(txt_file.replace('.txt','*.fa')))			
		def bps_check(bps):
			flag=0
			for x in range(len(bps)-2):
				if int(bps[x+2])-int(bps[x+1])>10**6:
					flag+=1
			return flag
		def cigar2alignstart(cigar,start,bps,flank_length):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[1])-flank_length: break
			return [read_rec,int(align_rec)-int(bps[1])+flank_length]
		def cigar2alignstart_2(cigar,start,bps):
			#eg cigar2alignstart(pbam[5],int(pbam[3]),bps)
			import re
			pcigar=re.compile(r'''(\d+)([MIDNSHP=X])''')
			cigars=[]
			for m in pcigar.finditer(cigar):
				cigars.append((m.groups()[0],m.groups()[1]))
			read_rec=0
			align_rec=start
			for x in cigars:
				if x[1]=='S':
					read_rec+=int(x[0])
				if x[1]=='M':
					read_rec+=int(x[0])
					align_rec+=int(x[0])
				if x[1]=='D':
					align_rec+=int(x[0])
				if x[1]=='I':
					read_rec+=int(x[0])
				if align_rec>int(bps[1]): break
			return [read_rec,int(align_rec)-int(bps[1])]
		def chop_pacbio_read(bps,flank_length,info):
			block_length={}
			for x in range(len(info[4:])-2):
				block_length[chr(97+x)]=int(info[x+6])-int(info[x+5])
			alA_len=numpy.sum([block_length[x] for x in info[3].split('/')[0] if not x=='^'])
			alB_len=numpy.sum([block_length[x] for x in info[3].split('/')[1] if not x=='^'])
			alRef_len=int(info[-1])-int(info[5])
			len_cff=max([alA_len,alB_len,alRef_len])
			if len_cff>10**4:
				len_cff=min([alA_len,alB_len,alRef_len])
			fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
			out=[]
			out2=[]
			out3=[]
			for line in fbam:
				pbam=line.strip().split()
				#out3.append(int(pbam[3])-int(bps[1])+flank_length)
				if not pbam[0]=='@': 
				  if int(pbam[3])<int(bps[1])-flank_length+1:
					align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
					align_start=align_info[0]
					miss_bp=align_info[1]
					align_pos=int(pbam[3])
					target_read=pbam[9][align_start:]
					if len(target_read)>flank_length+len_cff:
					  out.append(target_read[:max([alA_len,alB_len,alRef_len])+2*flank_length])
					  out2.append(miss_bp)
					  out3.append(pbam[0])
					  #out3.append(pbam[:9])
				  elif int(pbam[3])<int(bps[1])+100+1:
					align_info=cigar2alignstart_2(pbam[5],int(pbam[3]),bps)
					align_start=align_info[0]
					miss_bp=align_info[1]+flank_length
					target_read=pbam[9][align_start:]
					if len(target_read)>len_cff:
					  out.append(target_read[:max([alA_len,alB_len,alRef_len])+flank_length])
					  out2.append(miss_bp)
					  out3.append(pbam[0])
					  #out3.append(pbam[:9])
			fbam.close()
			return [out,out2,out3]
		def Capitalize_ref(fa2):
			fin=open(fa2)
			data=[]
			for line in fin:
				pin=line.strip().split()
				data+=pin
			fin.close()
			for x in range(1,len(data)):
				data[x]=data[x].upper()
			fin=open(fa2,'w')
			for x in data:
				print >>fin,x
			fin.close()
		def calcu_eu_dis_algorithm(SVelter_Name,info):
			#eg: info=k2[0][0]
			txt_rec_file_writing(out_path,SVelter_Name,info[2:])
			txt_file=out_path+SVelter_Name+'.txt'
			ref_sv=info[2]
			alt_sv=info[3]
			chrom=info[4]
			bps=info[4:]
			if (int(bps[-1])-int(bps[1]))>1000:
				delta=(int(bps[-1])-int(bps[1]))/10
			else:
				delta=50
			rsquare_ref=[]
			rsquare_alt1=[]
			rsquare_alt2=[]
			if bps_check(bps)==0:
					bl_len_hash={}
					for x in ref_sv.split('/')[0]:
						bl_len_hash[x]=int(bps[ord(x)-97+2])-int(bps[ord(x)-97+1])
					end_point=0
					for x in alt_sv.split('/')[0]:
						if not x=='^':
							end_point+=bl_len_hash[x]
					flank_length=flank_length_decide(bps)
					all_reads=chop_pacbio_read(bps,flank_length,info)
					read_hash=all_reads[0]
					miss_hash=all_reads[1]
					name_hash=name_hash_modify(all_reads[2])
					rec_len=0
					os.system(r'''Pacbio.produce.ref.alt.ref.py --fl %d --ref %s --sv %s'''%(flank_length,ref,txt_file))
					fa1=out_path+'temp.fa'
					Capitalize_ref(txt_file.replace('.txt','.ref.fa'))
					Capitalize_ref(txt_file.replace('.txt','.alt1.fa'))
					Capitalize_ref(txt_file.replace('.txt','.alt2.fa'))
					miss_rec=-1
					if not read_hash==[]:
							if len(read_hash)>20:
									new_read_hash=random.sample(range(len(read_hash)),20)
									read_hash2=[read_hash[i] for i in new_read_hash]
									miss_hash2=[miss_hash[i] for i in new_read_hash]
									name_hash2=[name_hash[i] for i in new_read_hash]
									read_hash=read_hash2
									miss_hash=miss_hash2
									name_hash=name_hash2
							name_rec=['0','0','0','0']
							for x in read_hash:
								miss_rec+=1
								y=x
								y2=miss_hash[miss_rec]
								y2_name=name_hash[miss_rec]
								fo=open(out_path+'temp.fa','w')
								print >>fo, '>temp'
								for z in chop_x_for_fasta(y):
										print >>fo, z
								fo.close()
								fo_ref=txt_file.replace('.txt','.dotplot.ref')
								fo1_alt=txt_file.replace('.txt','.dotplot.alt1')
								fo2_alt=txt_file.replace('.txt','.dotplot.alt2')
								os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,txt_file.replace('.txt','.ref.fa'),fo_ref))
								os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,txt_file.replace('.txt','.alt1.fa'),fo1_alt))
								os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,txt_file.replace('.txt','.alt2.fa'),fo2_alt))
								eu_dis_calcu(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
								#Smaller the dis is, better the prediction is
								if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
										min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
										rsquare_ref=rsquare_ref[:min_len]
										rsquare_alt1=rsquare_alt1[:min_len]
										rsquare_alt2=rsquare_alt2[:min_len]
								if min([abs(rsquare_alt1[-1]),abs(rsquare_alt2[-1])])-rsquare_ref[-1] <rec_len:
										rec_len=min([abs(rsquare_alt1[-1]),abs(rsquare_alt2[-1])])-abs(rsquare_ref[-1])
										rec_start=y2
										if abs(rsquare_alt1[-1])<abs(rsquare_alt2[-1]):
											os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest.allel1'))
											os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest.allel1'))
											os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest.allel1'))
											os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.allel1.fa')))	 
											name_rec[0]= y2_name
											name_rec[2]= y2
										else:
											os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest.allel2'))
											os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest.allel2'))
											os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest.allel2'))
											os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.allel2.fa')))	
											name_rec[1]= y2_name
											name_rec[3]= y2
							if not os.path.isfile(fo_ref+'longest.allel2'):
										os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest.allel2'))
										os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest.allel2'))
										os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest.allel2'))
										os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.allel2.fa')))	
							if not os.path.isfile(fo_ref+'longest.allel1'):						 
										os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest.allel1'))
										os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest.allel1'))
										os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest.allel1'))
										os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.allel1.fa')))								
							if os.path.isfile(txt_file.replace('.txt','.sample.allel1.fa')) and not name_rec[0]=='0':
								os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.allel1.fa'),txt_file.replace('.txt','.ref.fa'),fo_ref+'.allel1.'+name_rec[0]+'.png'))
								os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.allel1.fa'),txt_file.replace('.txt','.alt1.fa'),fo1_alt+'.allel1.'+name_rec[0]+'.png'))
								os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.allel1.fa'),txt_file.replace('.txt','.alt2.fa'),fo2_alt+'.allel1.'+name_rec[0]+'.png'))
								os.system(r'''cp %s %s '''%(fo_ref+'.allel1.'+name_rec[0]+'.png.txt',fo_ref+'.allel1.'+name_rec[0]+'.start.'+str(name_rec[2])))
								os.system(r'''cp %s %s '''%(fo1_alt+'.allel1.'+name_rec[0]+'.png.txt',fo1_alt+'.allel1.'+name_rec[0]+'.start.'+str(name_rec[2])))
								os.system(r'''cp %s %s '''%(fo2_alt+'.allel1.'+name_rec[0]+'.png.txt',fo2_alt+'.allel1.'+name_rec[0]+'.start.'+str(name_rec[2])))
							if os.path.isfile(txt_file.replace('.txt','.sample.allel2.fa')) and not name_rec[1]=='0':
								os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.allel2.fa'),txt_file.replace('.txt','.ref.fa'),fo_ref+'.allel2.'+name_rec[1]+'.png'))
								os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.allel2.fa'),txt_file.replace('.txt','.alt1.fa'),fo1_alt+'.allel2.'+name_rec[1]+'.png'))
								os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.allel2.fa'),txt_file.replace('.txt','.alt2.fa'),fo2_alt+'.allel2.'+name_rec[1]+'.png'))
								os.system(r'''cp %s %s '''%(fo_ref+'.allel2.'+name_rec[1]+'.png.txt',fo_ref+'.allel1.'+name_rec[1]+'.start.'+str(name_rec[3])))
								os.system(r'''cp %s %s '''%(fo1_alt+'.allel2.'+name_rec[1]+'.png.txt',fo1_alt+'.allel1.'+name_rec[1]+'.start.'+str(name_rec[3])))
								os.system(r'''cp %s %s '''%(fo2_alt+'.allel2.'+name_rec[1]+'.png.txt',fo2_alt+'.allel1.'+name_rec[1]+'.start.'+str(name_rec[3])))
			remove_files(txt_file)
			out1=[[1 for i in rsquare_ref],[abs(float(rsquare_alt1[i])/float(rsquare_ref[i])) for i in range(len(rsquare_ref)) if not rsquare_ref[i]==0],[abs(float(rsquare_alt2[i])/float(rsquare_ref[i])) for i in range(len(rsquare_ref)) if not rsquare_ref[i]==0]]
			out_al1=[[],[],[]]
			out_al2=[[],[],[]]
			for x in range(len(out1[2])):
					if out1[1][x]>out1[2][x]:
						out_al1[0].append(out1[0][x])
						out_al1[1].append(out1[1][x])
						out_al1[2].append(out1[2][x])
					else:
						out_al2[0].append(out1[0][x])
						out_al2[1].append(out1[1][x])
						out_al2[2].append(out1[2][x])
			return [out1[0],out_al1[2]+out_al2[1]]
		def multi_calcu_eu(Delly_Name,multi_info):
			#eg: multi_info=k2[1]
			Delly_Info=[[],[]]
			for k3 in multi_info:
				temp_info=k3+[k1]+[str(i) for i in k3[:2]]
				temp_result=calcu_eu_dis_algorithm(Delly_Name,temp_info)
				for x in range(len(Delly_Info)):
					Delly_Info[x]+=temp_result[x]
			return Delly_Info
		def heter_info_Filter(SVelter_Info):
			SVelter_Info2=[[],[]]
			for x in range(len(SVelter_Info[0])):
					if SVelter_Info[0][x]>SVelter_Info[1][x]:
							SVelter_Info2[0].append(SVelter_Info[0][x])
							SVelter_Info2[1].append(SVelter_Info[1][x])
			if not SVelter_Info2==[[],[]]:
					return SVelter_Info2
			else:
					return SVelter_Info
		def read_csv_in(filein):
			data_hash={}
			fin=open(filein)
			for line in fin:
				pin=line.strip().split()
				if not pin[0][0]=='#':
					if not pin[2] in data_hash.keys():
						data_hash[pin[2]]=[]
					if pin[2][-1]=='C':
						if '[' in pin[4] or ']' in pin[4]:
							data_hash[pin[2]].append([pin[0],pin[1],pin[4]])
					else:
						data_hash[pin[2]].append(pin[4:])
			fin.close()
			return data_hash
		def produce_structure_to_vali(data_hash,k1):
			out=[]
			for x in range(len(data_hash[k1])/2):
				if data_hash[k1][x][2][0] in ['a','t','g','c','n','A','T','G','C','N'] and data_hash[k1][x][2][-1] in ['[',']']:
					if '[' in data_hash[k1][x][2]:
						block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1])-flank_length,int(data_hash[k1][x][1])]
						block_b_temp=data_hash[k1][x][2].split('[')[1].split(':')
						block_b=[block_b_temp[0],int(block_b_temp[1]),int(block_b_temp[1])+flank_length]
						out.append(['a','b',block_a,block_b])
					elif ']' in data_hash[k1][x][2]:
						block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1])-flank_length,int(data_hash[k1][x][1])]
						block_b_temp=data_hash[k1][x][2].split(']')[1].split(':')
						block_b=[block_b_temp[0],int(block_b_temp[1])-flank_length,int(block_b_temp[1])]
						out.append( ['a','b^',block_a,block_b])
				if data_hash[k1][x][2][-1] in ['a','t','g','c','n','A','T','G','C','N'] and data_hash[k1][x][2][0] in ['[',']']:
					if '[' in data_hash[k1][x][2]:
						block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1]),int(data_hash[k1][x][1])+flank_length]
						block_b_temp=data_hash[k1][x][2].split('[')[1].split(':')
						block_b=[block_b_temp[0],int(block_b_temp[1]),int(block_b_temp[1])+flank_length]
						out.append( ['a','b^',block_a,block_b])
					elif ']' in data_hash[k1][x][2]:
						block_a=[data_hash[k1][x][0],int(data_hash[k1][x][1]),int(data_hash[k1][x][1])+flank_length]
						block_b_temp=data_hash[k1][x][2].split(']')[1].split(':')
						block_b=[block_b_temp[0],int(block_b_temp[1])-flank_length,int(block_b_temp[1])]
						out.append( ['a','b',block_a,block_b])
			return out
		def complementary(seq):
			seq2=[]
			for i in seq:
					if i in 'ATGCN':
							seq2.append('ATGCN'['TACGN'.index(i)])
					elif i in 'atgcn':
							seq2.append('atgcn'['tacgn'.index(i)])
			return ''.join(seq2)
		def reverse(seq):
			return seq[::-1]
		def chop_pacbio_read_long(bps):
			flank_length=0
			block_length={}
			len_cff=bps[-1]-bps[1]
			fbam=os.popen(r'''samtools view %s %s:%d-%d'''%(bam_in,bps[0],int(bps[1])-flank_length,int(bps[-1])+flank_length))
			out=[]
			out2=[]
			out3=[]
			for line in fbam:
					pbam=line.strip().split()
					#out3.append(int(pbam[3])-int(bps[1])+flank_length)
					if not pbam[0]=='@': 
						if int(pbam[3])<int(bps[1])-flank_length+1:
							align_info=cigar2alignstart(pbam[5],int(pbam[3]),bps,flank_length)
							align_start=align_info[0]
							miss_bp=align_info[1]
							align_pos=int(pbam[3])
							target_read=pbam[9][align_start:]
							if len(target_read)>flank_length+len_cff:
								out.append(target_read[:len_cff+2*flank_length])
								out2.append(miss_bp)
								out3.append(pbam[0])
								#out3.append(pbam[:9])
			fbam.close()
			return [out,out2,out3]
		def Pacbio_prodce_ref_alt_long(ref,new_structure_single,txt_file):
			#eg: new_structure_single=new_structure[0]
			x=new_structure_single
			struc_hash={}
			struc_hash[x[0][0]]=''
			struc_hash[x[0][0]+'^']=''
			fread1=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,x[2][0],x[2][1],x[2][2]))
			pin=fread1.readline().strip().split()
			for line in fread1:
				pin=line.strip().split()
				struc_hash[x[0][0]]+=pin[0]
			fread1.close()
			struc_hash[x[0][0]+'^']=reverse(complementary(struc_hash[x[0][0]]))
			struc_hash[x[1][0]]=''
			struc_hash[x[1][0]+'^']=''
			fread1=os.popen(r'''samtools faidx %s %s:%d-%d'''%(ref,x[3][0],x[3][1],x[3][2]))
			pin=fread1.readline().strip().split()
			for line in fread1:
				pin=line.strip().split()
				struc_hash[x[1][0]]+=pin[0]
			fread1.close()
			struc_hash[x[1][0]+'^']=reverse(complementary(struc_hash[x[1][0]]))
			fo1=open('.'.join(txt_file.split('.')[:-1])+'.alt.fa','w')
			print >>fo1, '_'.join([str(i) for i in [x[0],x[1]]+x[2]+x[3]])
			new_seq=''
			for y in x[:2]:
				new_seq+=struc_hash[y]
			for y in chop_x_for_fasta(new_seq):
				print >>fo1,y
			fo1.close()
			os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,x[2][0],x[2][1],x[2][1]+2*flank_length,'.'.join(txt_file.split('.')[:-1])+'.ref1.fa'))
			os.system(r'''samtools faidx %s %s:%d-%d > %s'''%(ref,x[3][0],x[3][1],x[3][1]+2*flank_length,'.'.join(txt_file.split('.')[:-1])+'.ref2.fa'))
			return ['.'.join(txt_file.split('.')[:-1])+'.alt.fa','.'.join(txt_file.split('.')[:-1])+'.ref1.fa','.'.join(txt_file.split('.')[:-1])+'.ref2.fa']
		def calcu_eu_dis_comp_long(SVelter_Name,info,new_structure):
			flank_length=500
			txt_rec_file_writing(out_path,SVelter_Name,info[2:])
			txt_file=out_path+SVelter_Name+'.txt'
			out=[[],[],[]]
			out_al_a=[[],[],[]]
			out_al_b=[[],[],[]]
			for x_structure in new_structure:
				fa_file_names=Pacbio_prodce_ref_alt_long(ref,x_structure,txt_file)
				ref_info1=[x_structure[2][0],x_structure[2][1],x_structure[2][1]+2*flank_length]
				ref_info2=[x_structure[3][0],x_structure[3][1],x_structure[3][1]+2*flank_length]
				delta=50
				rsquare_ref=[]
				rsquare_alt1=[]
				rsquare_alt2=[]
				if bps_check(ref_info1)==0:
					bps=ref_info1
					end_point=bps[-1]
					all_reads1=chop_pacbio_read_long(bps)
				if bps_check(ref_info2)==0:
					bps=ref_info2
					end_point=bps[-1]
					all_reads2=chop_pacbio_read_long(bps)
				read_hash=all_reads1[0]+all_reads2[0]
				miss_hash=all_reads1[1]+all_reads2[1]
				name_hash=name_hash_modify(all_reads1[2]+all_reads2[2])
				rec_len=0
				fa1=out_path+'temp.fa'
				Capitalize_ref(fa_file_names[0])
				Capitalize_ref(fa_file_names[1])
				Capitalize_ref(fa_file_names[2])
				miss_rec=-1
				if not read_hash==[]:
					if len(read_hash)>20:
							new_read_hash=random.sample(range(len(read_hash)),20)
							read_hash2=[read_hash[i] for i in new_read_hash]
							miss_hash2=[miss_hash[i] for i in new_read_hash]
							name_hash2=[name_hash[i] for i in new_read_hash]
							read_hash=read_hash2
							miss_hash=miss_hash2
							name_hash=name_hash2
					name_rec=['0','0']
					for x in read_hash:
						if len(x)>1000:
							x=x[:1000]
						miss_rec+=1
						y=x
						y2=miss_hash[miss_rec]
						y2_name=name_hash[miss_rec]
						fo=open(out_path+'temp.fa','w')
						print >>fo, '>temp'
						for z in chop_x_for_fasta(y):
								print >>fo, z
						fo.close()
						fo_ref=txt_file.replace('.txt','.dotplot.alt')
						fo1_alt=txt_file.replace('.txt','.dotplot.ref1')
						fo2_alt=txt_file.replace('.txt','.dotplot.ref2')
						os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,txt_file.replace('.txt','.alt.fa'),fo_ref))
						os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,txt_file.replace('.txt','.ref1.fa'),fo1_alt))
						os.system(r'''dotdata.py %d %s %s %s'''%(window_size_1,fa1,txt_file.replace('.txt','.ref2.fa'),fo2_alt))
						eu_dis_calcu_complex_long(fo_ref,fo1_alt,fo2_alt,rsquare_ref,rsquare_alt1,rsquare_alt2,y2,delta,info)
						#Smaller the dis is, better the prediction is
						if not len(rsquare_ref)==len(rsquare_alt1)==len(rsquare_alt2):
								min_len=min([len(rsquare_ref),len(rsquare_alt1),len(rsquare_alt2)])
								rsquare_ref=rsquare_ref[:min_len]
								rsquare_alt1=rsquare_alt1[:min_len]
								rsquare_alt2=rsquare_alt2[:min_len]
						if not rsquare_ref==[]:
							if rsquare_ref[-1]-max([abs(rsquare_alt1[-1]),abs(rsquare_alt2[-1])])<rec_len:
								rec_len=rsquare_ref[-1]-max([abs(rsquare_alt1[-1]),abs(rsquare_alt2[-1])])
								rec_start=y2
								os.system(r'''cp %s %s'''%(fo_ref+'.txt',fo_ref+'longest'))
								os.system(r'''cp %s %s'''%(fo1_alt+'.txt',fo1_alt+'longest'))
								os.system(r'''cp %s %s'''%(fo2_alt+'.txt',fo2_alt+'longest'))
								os.system(r'''cp %s %s'''%(fa1,txt_file.replace('.txt','.sample.fa')))	 
								name_rec[1]=y2
								name_rec[0]=y2_name
					if not os.path.isfile(fo_ref+'longest'):
						print 'Cannot Validate'
					else:
						test_fin=os.popen(r'''wc -l %s'''%(txt_file.replace('.txt','.alt.fa')))
						test_pin=test_fin.readline().strip().split()
						test_fin.close()
						if not test_pin[0]=='1':
							os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),txt_file.replace('.txt','.alt.fa'),fo_ref+'.'+name_rec[0]+'.png'))
							os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),txt_file.replace('.txt','.ref1.fa'),fo1_alt+'.'+name_rec[0]+'.png'))
							os.system(r'''dotplot.py %d %s %s %s'''%(window_size,txt_file.replace('.txt','.sample.fa'),txt_file.replace('.txt','.ref2.fa'),fo2_alt+'.'+name_rec[0]+'.png'))
							os.system(r'''mv %s %s'''%(fo_ref+'.'+name_rec[0]+'.png.txt',fo_ref+'.'+name_rec[0]+'.start.'+str(name_rec[1])))
							os.system(r'''mv %s %s'''%(fo1_alt+'.'+name_rec[0]+'.png.txt',fo1_alt+'.'+name_rec[0]+'.start.'+str(name_rec[1])))
							os.system(r'''mv %s %s'''%(fo2_alt+'.'+name_rec[0]+'.png.txt',fo2_alt+'.'+name_rec[0]+'.start.'+str(name_rec[1])))
				remove_files_short(txt_file)
				out1=[[1 for i in rsquare_ref],[abs(float(rsquare_alt1[i])/float(rsquare_ref[i])) for i in range(len(rsquare_ref)) if not rsquare_ref[i]==0],[abs(float(rsquare_alt2[i])/float(rsquare_ref[i])) for i in range(len(rsquare_ref)) if not rsquare_ref[i]==0]]
				out_al1=[[],[],[]]
				out_al2=[[],[],[]]
				for x in range(len(out1[2])):
						if out1[1][x]>out1[2][x]:
							out_al1[0].append(out1[0][x])
							out_al1[1].append(out1[1][x])
							out_al1[2].append(out1[2][x])
						else:
							out_al2[0].append(out1[0][x])
							out_al2[1].append(out1[1][x])
							out_al2[2].append(out1[2][x])
				out[0] +=out1[0]
				out_al_a[2]+=out_al1[2]
				out_al_b[1]+=out_al2[1]
			return [out_al_a[2]+out_al_b[1],out[0]]
		opts,args=getopt.getopt(sys.argv[2:],'i:',['vcf-file=','vali-file=','window-size=','reference=','bam-file=','output-prefix=','output-path=','sv=','path-lumpy=','path-pindel=','path-delly='])
		dict_opts=dict(opts)
		out_prefix=dict_opts['--output-prefix']
		if not '--output-path' in dict_opts.keys():
			out_path='/'.join(out_prefix.split('/')[:-1])
		else:
			out_path=dict_opts['--output-path']
		if not out_path[-1]=='/':
			out_path+='/'
		if not os.path.isdir(out_path):
			os.system(r'''mkdir %s'''%(out_path))
		vcf_in=dict_opts['--vcf-file']
		vcf_data_hash=read_csv_in(vcf_in)
		start=0
		if '--window-size' in dict_opts.keys():
			window_size=int(dict_opts['--window-size'])
		else:
			window_size=10
		bam_in=dict_opts['--bam-file']
		ref=dict_opts['--reference']
		min_length=50
		flank_length=100
		min_read_compare=20
		if '--path-delly' in dict_opts.keys():
			Path_Delly=dict_opts['--path-delly']
			Path_Delly=path_name_check(Path_Delly)
		else:
			Path_Delly='Null'
		if '--path-lumpy' in dict_opts.keys():
			Path_Lumpy=dict_opts['--path-lumpy']
			Path_Lumpy=path_name_check(Path_Lumpy)
		else:
			Path_Lumpy='Null'
		if '--path-pindel' in dict_opts.keys():
			Path_Pindel=dict_opts['--path-pindel']
			Path_Pindel=path_name_check(Path_Pindel)
		else:
			Path_Pindel='Null'
		SVelter_in=dict_opts['--vali-file']
		sample_name=SVelter_in.split('/')[-1].split('.')[0]
		out_path=path_name_check(out_path)
		SVelter_hash=csv_info_readin(SVelter_in)
		Delly_hash=Other_Algorithm_result_readin(Path_Delly)
		Lumpy_hash=Other_Algorithm_result_readin(Path_Lumpy)
		Pindel_hash=Other_Algorithm_result_readin(Path_Pindel)
		comparison=compare_hash()
		case_rec=0
		chromos=[]
		fref=open(dict_opts['--reference']+'.fai')
		for line in fref:
			pref=line.strip().split()
			chromos.append(pref[0])
		fref.close()
		fo_fina=open(SVelter_in+'.compare.stat','w')
		print >>fo_fina, ' '.join(['SVelter','Delly','Lumpy','Ref','info'])
		for k1 in comparison.keys():
			for k2 in comparison[k1]:
				flag_run=0
				for k3 in k2[0][0][5:]:
					if k3 in chromos:
						flag_run+=1
				if flag_run==0:
					case_rec+=1
					SVelter_Name='.'.join(['SVelter',sample_name,'_'.join(k2[0][0][4:])])
					Delly_Name='.'.join(['Delly',sample_name,'_'.join(k2[0][0][4:])])
					Lumpy_Name='.'.join(['Lumpy',sample_name,'_'.join(k2[0][0][4:])])
					Pindel_Name='.'.join(['Pindel',sample_name,'_'.join(k2[0][0][4:])])
					if not int(k2[0][0][-1])-int(k2[0][0][5])>10**4:
						SVelter_Info=calcu_eu_dis_algorithm(SVelter_Name,k2[0][0])
						fo=open(out_path+SVelter_Name+'.SVelter.stat','w')
						for x in range(len(SVelter_Info[0])):
							print>>fo, ' '.join([str(cm) for cm in [SVelter_Info[0][x],SVelter_Info[1][x]]])
						fo.close()
					else:
						if '_'.join(k2[0][0][4:]+['C']) in vcf_data_hash.keys():
							new_structure=produce_structure_to_vali(vcf_data_hash,'_'.join(k2[0][0][4:]+['C']))
						elif '_'.join(k2[0][0][4:]+['S']) in vcf_data_hash.keys():
							new_structure=produce_structure_to_vali(vcf_data_hash,'_'.join(k2[0][0][4:]+['S']))									
						else:
							new_structure=[]
						SVelter_Info=calcu_eu_dis_comp_long(SVelter_Name,k2[0][0],new_structure)
						fo=open(out_path+SVelter_Name+'.SVelter.stat','w')
						for x in range(len(SVelter_Info[0])):
							print>>fo, ' '.join([str(cm) for cm in [SVelter_Info[0][x],SVelter_Info[1][x]]])
						fo.close()
					if not Path_Delly=='Null':
						Delly_Info=multi_calcu_eu(Delly_Name,k2[1])
					if not Path_Lumpy=='Null':
						Lumpy_Info=multi_calcu_eu(Lumpy_Name,k2[2])
					if not Path_Pindel=='Null':
						Pindel_Info=multi_calcu_eu(Pindel_Name,k2[3])
					SVelter_Info2=heter_info_Filter(SVelter_Info)
					SVelter_Stat=numpy.mean(SVelter_Info2[1])
					Delly_Stat=1
					Lumpy_Stat=1
					Pindel_Stat=1
					if not Path_Delly=='Null':
						if not Delly_Info[0]==[]:
							Delly_Info2=heter_info_Filter(Delly_Info)
							Delly_Stat=numpy.mean(Delly_Info2[1])
					if not Path_Lumpy=='Null':
						if not Lumpy_Info[0]==[]:
							Lumpy_Info2=heter_info_Filter(Lumpy_Info)
							Lumpy_Stat=numpy.mean(Lumpy_Info2[1])
					if not Path_Pindel=='Null':
						if not Pindel_Info[0]==[]:
							Pindel_Info2=heter_info_Filter(Pindel_Info)
							Pindel_Stat=numpy.mean(Pindel_Info2[1])
					print >>fo_fina, ' '.join([str(i) for i in [1-SVelter_Stat,1-Delly_Stat,1-Lumpy_Stat,1-Pindel_Stat]+['_'.join([str(j) for j in [k1]+k2[0][0][:4]])]])
		fo_fina.close()
	if function_name=='svelter-to-vali':
		import os
		import sys
		import getopt
		import re
		opts,args=getopt.getopt(sys.argv[2:],'i:',['qc-structure=','input=','output=','reference=','vali-file=','output-prefix=','output-path='])
		dict_opts=dict(opts)
		input_file=dict_opts['--input']
		if not '--output' in dict_opts.keys():
			if input_file.split('.')[-1]=='svelter':
				output_file=input_file.replace('.svelter','.vali')
			else:
				output_file=input_file+'.vali'
		else:
			out_path=dict_opts['--output']
		if not '--qc-structure' in dict_opts.keys():
			qc_score=-20
		else:
			qc_score=float(dict_opts['--qc-structure'])
		fin=open(input_file)
		fo=open(output_file,'w')
		pin=fin.readline().strip().split()
		for line in fin:
			pin=line.strip().split()
			if float(pin[-1])>qc_score:
				print >>fo, '\t'.join([pin[-3],pin[-2]]+pin[-4].split(':'))
		fin.close()
	if function_name=='integrate-simple':
		import getopt
		import re
		import pickle
		import time
		import datetime
		import random
		def rsquare_file_analysis(rsquare_file):
			fin=open(rsquare_file)
			pin=fin.readline().strip().split()
			if not pin[0]=='alt' or not pin[1]=='ref':
				return 'Error'
			else:
				total_num=0
				pass_qc_mum=0
				for line in fin:
					pin=line.strip().split()
					total_num+=1
					if float(pin[0])>float(pin[1]):
						pass_qc_mum+=1
			fin.close()
			if not total_num==0:
				return float(pass_qc_mum)/float(total_num)
			else:
				return 'Nan'
		def txt_file_analysis(txt_file):
			fin=open(txt_file)
			pin=fin.readline().strip().split()
			fin.close()
			return pin
		opts,args=getopt.getopt(sys.argv[2:],'i:',['vcf-file=','vali-file=','window-size=','reference=','bam-file=','output-prefix=','output-path=','sv=','path-lumpy=','path-pindel=','path-delly='])
		dict_opts=dict(opts)
		out_prefix=dict_opts['--output-prefix']
		if not '--output-path' in dict_opts.keys():
			out_path='/'.join(out_prefix.split('/')[:-1])
		else:
			out_path=dict_opts['--output-path']
		if not out_path[-1]=='/':
			out_path+='/'
		fo=open(out_prefix+'.stat','w')
		print >>fo, ' '.join(['score','ref','alt','pos'])
		for k1 in os.listdir(out_path):
			if k1.split('.')[-1]=='rsquare':
				qc_score_structure=rsquare_file_analysis(out_path+k1)
				info=txt_file_analysis(out_path+'.'.join(k1.split('.')[:-1]+['txt']))
				print >>fo, ' '.join([str(i) for i in [qc_score_structure]+info[:2]+['_'.join(info[2:])]])
		fo.close()
	if function_name=='show-stats':
		import os
		import sys
		import getopt
		if len(sys.argv)==2:
			print 'Pacbio.Validate.py show-stats'
			print 'Parameters:'
			print '--input: input file with pacbio validation stat'
			print '--cutoff: minimum score required to be considered as validated. default: 0.3'
		else:
			opts,args=getopt.getopt(sys.argv[2:],'i:',['input=','cutoff='])
			dict_opts=dict(opts)
			if '--cutoff' in dict_opts.keys():
				cutoff=float(dict_opts['--cutoff'])
			else:
				cutoff=0.3
			filein=dict_opts['--input']
			fin=open(filein)
			total_num=0
			vali_num=0
			not_vali_num=0
			cannot_vali_num=0
			for line in fin:
				pin=line.strip().split()
				if not pin[0]=='score' and not pin[0]=='SVelter':
					total_num+=1
					if not pin[0]=='Nan':
						if float(pin[0])>cutoff:
							vali_num+=1
						else:
							not_vali_num+=1
					else:
						cannot_vali_num+=1
			fin.close()
			print ' '.join(['validated','Non-validated','Cannot-validate','total-num'])
			print ' '.join([str(i) for i in [vali_num,not_vali_num,cannot_vali_num,total_num]])
	if function_name=='svs-to-vali':
		import os
		import sys
		import getopt
		opts,args=getopt.getopt(sys.argv[2:],'i:',['sv-type=','svelter-file=','qc-structure='])
		dict_opts=dict(opts)
		qc_score=-20
		if '--qc-structure' in dict_opts.keys():
				qc_score=float(dict_opts['--qc-structure'])
		sv_type_hash={}
		fin=open(dict_opts['--sv-type'])
		for line in fin:
				pin=line.strip().split()
				if not pin[0] in sv_type_hash.keys():
						sv_type_hash[pin[0]]={}
				if not pin[1] in sv_type_hash[pin[0]].keys():
						sv_type_hash[pin[0]][pin[1]]=[]
		fin.close()
		if dict_opts['--svelter-file'].split('.')[-1]=='svelter':
				fin=open(dict_opts['--svelter-file'])
				pin=fin.readline().strip().split()
				for line in fin:
						pin=line.strip().split()
						if pin[-3] in sv_type_hash.keys():
								if pin[-2] in sv_type_hash[pin[-3]].keys():
										if float(pin[-1])>qc_score:
												sv_type_hash[pin[-3]][pin[-2]].append(pin[-4])
				fin.close()
		elif dict_opts['--svelter-file'].split('.')[-1]=='rec':
				fin=open(dict_opts['--svelter-file'])
				pin=fin.readline().strip().split()
				for line in fin:
						pin=line.strip().split()
						if pin[0] in sv_type_hash.keys():
								if pin[1] in sv_type_hash[pin[0]].keys():
										if float(pin[-1])>qc_score:
												sv_type_hash[pin[0]][pin[1]].append(':'.join(pin[2:-1]))
				fin.close()
		fo=open(dict_opts['--sv-type']+'.vali','w')
		for k1 in sorted(sv_type_hash.keys()):
				for k2 in sorted(sv_type_hash[k1].keys()):
						for k3 in sv_type_hash[k1][k2]:
							check=k3.split(':')
							flag_check=0
							for x in range(len(check))[1:-1]:
								if 'chr' in check[x] or 'M' in check[x] or 'Y' in check[x] or 'X' in check[x]:
									flag_check+=1
								elif 'chr' in check[x+1] or 'M' in check[x+1] or 'Y' in check[x+1] or 'X' in check[x+1]:
									flag_check+=1
								else:
									if int(check[x+1])<int(check[x]):
										flag_check+=1
							if flag_check==0:
								print >>fo, '\t'.join([k1,k2]+k3.split(':'))
		fo.close()
	if function_name=='split-vali':
		import os
		import sys
		import getopt
		opts,args=getopt.getopt(sys.argv[2:],'i:',['input=','size='])
		dict_opts=dict(opts)
		filein=dict_opts['--input']
		if '--size' in dict_opts.keys():
			sub_size=int(dict_opts['--size'])
		else:
			sub_size=100
		sub_rec=[[]]
		fin=open(filein)
		for line in fin:
			pin=line.strip().split()
			if len(sub_rec[-1])<sub_size:
				sub_rec[-1].append(pin)
			else:
				sub_rec.append([pin])
		fin.close()
		rec=0
		for x in sub_rec:
			rec+=1
			fin=open('.'.join(filein.split('.')[:-1])+'.'+str(rec)+'.vali','w')
			for y in x:
				print >>fin, '\t'.join(y)
			fin.close()




#Not in use
		def rsquare_calcu(fo_ref,rsquare_ref):
			fo1=open(fo_ref+'.txt')
			temp_data1=[[],[]]
			for line in fo1:
				po1=line.strip().split()
				temp_data1[0].append(int(po1[0]))
				temp_data1[1].append(int(po1[1]))
			fo1.close()
			slope, intercept, r_value, p_value, std_err = stats.linregress(temp_data1[0],temp_data1[1])
			rsquare_ref.append(r_value**2)
