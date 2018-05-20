# By NYU.
# Free program for academic purpose.
# This program is the main body of GaussianCpG.
# It will generate a bundle of files for CpG's locations, distances between CpG locations, statistics of distances (index of loc includes), predicted file.
# Its input parameters include the theshold of CpG island distance for clustering, debug switches, repeat or non-repeat search.

# Also, it can generate some other debug files. 
# showgreaterthreshold.sort.txt  --> all clustered CpG-rich areas
# gaussian.txt --> pre-generated Gaussian Filter.
# cpg_i_chr21.txt.init.txt --> debug the clustered CpG rich areas.
# cluster_sub.txt --> entropy test. Human whole genome's dinucleotide percentage transfer to entropy and test the cluster to see if they can beat the entropy.
# cpg_i_chr21.txt.result.txt --> 


# cpgreport: 
# http://imed.med.ucm.es/cgi-bin/emboss.pl?_action=input&_app=cpgreport


# Cpgcluster:
# http://bioinfo2.ugr.es/CpGislands/


# cpgProd:
#  http://pbil.univ-lyon1.fr/software/cpgprod_query.php


# cpgplot:
# http://www.hpa-bioinfotools.org.uk/pise/cpgplot.html


# human assembled chromosome
# ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/Assembled_chromosomes/

# human CpG island database
# http://data.microarrays.ca/cpg/

# mouse assembled chromosome
# ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/Assembled_chromosomes/seq/

# mouse CpG island database
# http://data.microarrays.ca/cpgmouse/index.htm

# refseq TSS promoter
# ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/protein/

# Epigenetics database resource
# http://epigenie.com/epigenetic-tools-and-databases/
# http://wanglab.ucsd.edu/star/epigram/


# python CpGs.py artificial_dna.txt cpg_loc_artifi.txt cpg_dis_loc_artifi.txt cpg_stat_dis_artifi.txt 1 1 1  1 95 cpg_i_artifi.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr21.dna.txt cpg_loc_chr21.txt cpg_dis_loc_chr21.txt cpg_stat_dis_chr21.txt 1 1 1  1 95 cpg_i_chr21.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr21.dna.txt cpg_loc_chr21.txt cpg_dis_loc_chr21.txt cpg_stat_dis_chr21.txt 1 1 1  1 95 cpg_i_chr21.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr21.dna.txt cpg_loc_chr21.txt cpg_dis_loc_chr21.txt cpg_stat_dis_chr21.txt 1 1 1  0 95 cpg_i_chr21.txt



# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/dataset_cpgprod_sub1.dna.txt cpg_loc_dataset_cpgprod_sub1.txt cpg_dis_loc_dataset_cpgprod_sub1.txt cpg_stat_dis_dataset_cpgprod_sub1.txt 1 1 1  1 95 cpg_i_dataset_cpgprod_sub1.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/dataset_cpgprod_sub2.dna.txt cpg_loc_dataset_cpgprod_sub2.txt cpg_dis_loc_dataset_cpgprod_sub2.txt cpg_stat_dis_dataset_cpgprod_sub2.txt 1 1 1  1 95 cpg_i_dataset_cpgprod_sub2.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/dataset_cpgprod_sub3.dna.txt cpg_loc_dataset_cpgprod_sub3.txt cpg_dis_loc_dataset_cpgprod_sub3.txt cpg_stat_dis_dataset_cpgprod_sub3.txt 1 1 1  1 95 cpg_i_dataset_cpgprod_sub3.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/dataset_cpgprod_sub4.dna.txt cpg_loc_dataset_cpgprod_sub4.txt cpg_dis_loc_dataset_cpgprod_sub4.txt cpg_stat_dis_dataset_cpgprod_sub4.txt 1 1 1  1 95 cpg_i_dataset_cpgprod_sub4.txt
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/dataset_cpgprod_sub5.dna.txt cpg_loc_dataset_cpgprod_sub5.txt cpg_dis_loc_dataset_cpgprod_sub5.txt cpg_stat_dis_dataset_cpgprod_sub5.txt 1 1 1  1 95 cpg_i_dataset_cpgprod_sub5.txt


# python CpGs.py mouse_extracted_seqs_not_fasta/chr1.txt cpg_loc_chr1.m.txt cpg_dis_loc_chr1.m.txt cpg_stat_dis_chr1.m.txt 1 1 1  1 95 cpg_i_chr1.m.txt
# python CpGs.py artificial_dna_mouse.txt cpg_loc_artifi.m.txt cpg_dis_loc_artifi.m.txt cpg_stat_dis_artifi.m.txt 1 1 1  1 95 cpg_i_artifi.m.txt
# python CpGs.py mouse_extracted_seqs_not_fasta/chr19.txt cpg_loc_chr19.m.txt cpg_dis_loc_chr19.m.txt cpg_stat_dis_chr19.m.txt 1 1 1  1 95 cpg_i_chr19.m.txt


# 07/24/2015

# study chr21 as the basic one.
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr21.dna.txt cpg_loc_chr21.txt cpg_dis_loc_chr21.txt cpg_stat_dis_chr21.txt 1 1 1  1 95 cpg_i_chr21.txt 0

# for cpg in repeat region only.
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr1.dna.txt cpg_loc_chr1.txt cpg_dis_loc_chr1.txt cpg_stat_dis_chr1.txt 1 1 1  1 95 cpg_i_chr1.txt 2
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr1.dna.txt cpg_loc_chr1.txt cpg_dis_loc_chr1.txt cpg_stat_dis_chr1.txt 1 1 1  1 95 cpg_i_chr1.txt 0
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr21_hg38.dna.txt cpg_loc_chr21.txt cpg_dis_loc_chr21.txt cpg_stat_dis_chr21.txt 1 1 1  1 95 cpg_i_chr21.txt 0
# python CpGs.py extracted\ dna\ sqs\ \(not\ fasta\)/chr22_hg38.dna.txt cpg_loc_chr22.txt cpg_dis_loc_chr22.txt cpg_stat_dis_chr22.txt 1 1 1  1 120 cpg_i_chr22.txt 0




import sys, os
import math


dna_len = 0;

min_region_len = 8;


#This function only reads the files extracted from RawFasta2CppReadablefile.py .
#Such files have no any punctions and special symbols, even no new line, only letters.
def read_dna_file(readfilename):
	

	with open(filename) as file:
		contents = file.read()
	
	#split function split the whole contents into pieces.
	entries = contents.split('>')[1:]
	
	#partition function partition each entry into 3 parts, the middle one is '\n'
	#partition function is different from split.
	partitioned_entries = [entry.partition('\n') for entry in entries]
	pairs = [(entry[0],entry[2]) for entry in partitioned_entries]
	pairs2 = [(pair[0], pair[1].replace('\n','')) for pair in pairs]
	result = [(pair[0].split('|'), pair[1]) for pair in pairs2]
	return result
	
	with open(readfilename) as file:
		for line in file: 
			if line.startswith('Input1'):
				
				parts = line.split(': (');
				part = parts[1].split(',');
				if part[0] == 'Input2':
					reverse_flag = 1;
	
			elif line.startswith('path '):
				
				parts = line.split('|');
				hitstr = parts[4].replace(' [', '');
				subhit = hitstr.split(' ');		#hit=20 gap=5584 err=1857
				myhit = subhit[0].split('=');
				hit = myhit[1];
				
				bigstr = parts[2].replace(' ','');
				bigstr = bigstr.replace('x:','');
				bigstr = bigstr.replace('y:','');
				part = bigstr.split(',');
				if part == ['']:
					continue;
				else: 
					firstpair = part[0].split('..');
					secondpair = part[1].split('..');
						
					if reverse_flag == 0:
						element.append([int(firstpair[0]),int(firstpair[1]),int(secondpair[0]),int(secondpair[1]), int(hit)]);
					else :
						element.append([int(secondpair[0]),int(secondpair[1]), int(firstpair[0]),int(firstpair[1]), int(hit)]);
	return element;



readfilename = sys.argv[1];	
outputfilename  = sys.argv[2];			#cpg location file records the locations of all CpG sites.
distancefilename = sys.argv[3];			#this file stores the distances between every two neighboring CpG sites, following the first outputfile.
statisticfilename = sys.argv[4];		#statistic is used to store the counts for specific distance, and the index in distance array.
outputfile_on = sys.argv[5];			
distancefile_on = sys.argv[6];
statisticfile_on = sys.argv[7];

threshold_on = sys.argv[8];
threshold = int(sys.argv[9]);

clusterfilename = sys.argv[10];

repeat_on = (sys.argv[11]);

#csize = int(sys.argv[4]);
#number = int(sys.argv[5]);

contents = [];
with open(readfilename) as file:
	contents = file.read()


print len(contents);




################ to find all the CGs in DNA genome.
last = '';
loci = [];		#loci is used to store all CpG's locations.
for i, c in enumerate(contents):
	if last != '':
		# this is for the statistical data for genomic analysis..
		if repeat_on == '1':
			if (last == 'C' or last == 'c') and (c =='G' or c =='g'):
				loci.append(i);
		
		#ignore the repeat segment.		
		if repeat_on == '0':
		#this is for prediction of CGI.
			if (last == 'C' ) and (c =='G' ):
				loci.append(i);
		
		#repeat area only or semi-repeat areas.
		if repeat_on == '2':
			if (c=='g' and last == 'c') or (c =='g' and last == 'C') or (c =='G' and last == 'c') :
				loci.append(i);

		
	#if c != 'N' and last == 'N' and i < 10000000:
	#	print "if c != 'N' and last == 'N' and i < 10000000:", i;
	last = c;


###############Output all loci of CpGs to file.
###############to write the loci of CgGs.
#	cpg_loc_chr21.txt 
fd = open("./intermediate/{}".format(outputfilename), 'w');
#fd.write("{}\n".format(len(contents)));

buf = "";
j = 0;

if outputfile_on == '1':
	for k in range(len(loci)):
		if j  < 4000:
			buf += "{}\n".format(loci[k]);
		else:
			fd.write(buf);
			buf = "";
			buf += "{}\n".format(loci[k]);
			j = 0;
		j += 1;	
	
fd.write(buf);

del buf;
fd.close();


print 'The number of all CpGs ->', len(loci);


##############to get the difference/distance, between loci.
diff = [];
last = 0;
#dif[1] is the index of loci array.
for i, loc in enumerate(loci):
	if i != 0:	
		dif = loc - last;
		diff.append([dif,i]);
	last = loc;

#diff contains all differences/distances between two neighbouring CpGs.
diff.sort();


#write these distances between two CpGs to file.
#eg. cpg_dis_loc_chr21.txt for (diff)
# all [distance, location]
disfile = open("./intermediate/{}".format(distancefilename), 'w');

buf = "";
j = 0;
if distancefile_on == '1':
	for dif in diff:
		if j  < 4000:
			buf += "{}\n".format(dif);
		else:
			disfile.write(buf);
			buf = "";
			buf += "{}\n".format(dif);
			j = 0;
		j += 1;		
disfile.write(buf);	
disfile.close();




#Obtaining the statistic for the purpose of clustering them.
#Note that: this aims to do the statistics about the distances and the # of distances.
#laststat [0] is the distance between CpGs. laststat[1] is the number of this distance, [2] is the first one's index in diff.
#statistic[][0] and [][1] are same as laststat[0] and [1]. statistic[][2] is the index of the first laststat[0] in diff.
#index_diff is a map whose key is the distance of CpGs, whose value is the index of diff corresponding to the first distance..
statistic = [];
laststat = [0,1,0];		

index_diff = {};



print "the number of all distances is:  ", len(diff);
#count the number of each distance. 

for i, dif in enumerate(diff):
	if (dif[0] != laststat[0]) and (laststat[0] != 0):
		statistic.append([laststat[0],laststat[1],laststat[2]]);
		index_diff[laststat[0]] = laststat[2];
		#print last;
		laststat[1] = 1;	# The counter is set to 1.
		laststat[2] = i;	# the first one's index.
	elif (dif[0] == laststat[0]) and (laststat[0] != 0):
		laststat[1] += 1;

	#disfile.write("{}\n".format(dif));
	laststat[0] = dif[0];
	
	#if it is the last one, must dump its memory.
	if i == len(diff)-1:
		statistic.append([laststat[0],laststat[1],laststat[2]]);
		index_diff[laststat[0]] = laststat[2];	
	

	
#print statistic;


#print index_diff;

#write the distance, the statistical counter and the first index to the statistic file.
# eg. cpg_stat_dis_chr21.txt
# [distance, counter, first index in (diff)]
statifile = open("./intermediate/{}".format(statisticfilename), 'w+');

counter = 0;
pflag50 = 0;
pflag8585 = 0;
pflag8595 = 0;
pflag8599 = 0;

pflag683 = 0;
pflag85 = 0;
pflag95 = 0;
pflag99 = 0;


buf = "";
j = 0;

if statisticfile_on == '1':
	for st in statistic:
		if j  < 4000:
			buf += "{}\t{}\t{}\n".format(st[0],st[1],st[2]);
		else:
			statifile.write(buf);
			buf = "";
			buf += "{}\t{}\t{}\n".format(st[0],st[1],st[2]);
			j = 0;
		j += 1;			

statifile.write(buf);
statifile.close();



#It aims to print out some important statistics of distance x and counters y.
#Percentage is important to identify its distribution function. 
#Note that: st[0] is the distance; st[1] is the # of this distance.
for st in statistic:
	#if j  < 4000:
	#	buf += "{}\t{}\n".format(st[0],st[1]);
	#else:
	#	statifile.write(buf);
	#	buf = "";
	#	j = 0;
	#j += 1;			
	#statifile.write("{}\t{}\n".format(st[0],st[1]));
	counter += st[1];
	if (st[0] 	 == 27):
		print "distance 27 accumulated %-->", float(counter)/float(len(diff))
	if (st[0] == 28):
		print "distance 28 accumulated %-->", float(counter)/float(len(diff))
		
	if (float(counter)/float(len(diff)) > 0.50) and (pflag50 == 0):
		print "0.50->distance ",st;
		pflag50 = 1;	
	if (float(counter)/float(len(diff)) > float(0.85*0.85)) and (pflag8585 == 0):
		print "0.85*0.85=",float(0.85*0.85),'->',st;
		pflag8585 = 1;
	if (float(counter)/float(len(diff)) > float(0.85*0.95)) and (pflag8595 == 0):
		print "0.85*0.95=",float(0.85*0.95),'->',st;
		pflag8595 = 1;
	if (float(counter)/float(len(diff)) > float(0.85*0.99)) and (pflag8599 == 0):
		print "0.85*0.99=",float(0.85*0.99),'->',st;
		pflag8599 = 1;
	if (float(counter)/float(len(diff)) > 0.683) and (pflag683 == 0):
		print "0.683->",st;
		pflag683 = 1;
	elif (float(counter)/float(len(diff)) > 0.85) and (pflag85 == 0):
		print "0.85->",st;
		pflag85 = 1;
	elif (float(counter)/float(len(diff)) > 0.95) and (pflag95 == 0):
		print "0.95->",st;
		pflag95 = 1;
	elif (float(counter)/float(len(diff)) > 0.99) and (pflag99 == 0):
		print "0.99->",st;
		pflag99 = 1;

#statifile.write(buf);
#statifile.close();

#os.close(outputfilename);
#read_dna_file(readfilename);

if threshold_on == '0':
	print "off threshold";
	exit();


#print index_diff;

#if int(threshold) > len(index_diff):
if index_diff.get(threshold) == None :	
	int_threshold = int(threshold)-1;
	while (index_diff.get(threshold) == None ) and (int_threshold > 0) :
		int_threshold -= 1;
		threshold = (int_threshold);
		#print int_threshold, threshold, index_diff.get(threshold);
		#print index_diff.keys();
	
	print "threshold changed ", threshold, "<-",len(index_diff)-1;

	

#index_diff is a map. the key is distance of 
#print threshold, len(index_diff), index_diff[int(threshold)];
print "index of threshold {} is {}".format(threshold, index_diff[int(threshold)]);#, index_diff[3], index_diff[4];

#any one who is greater than threshold should be in the following array.
diff_threshold = [];
for dif in diff[(index_diff[int(threshold)]):-1]:
	#dif[1] is the index of (loci) array. dif[0] is the distance between two loci.
	diff_threshold.append([dif[1],dif[0]]);
# The last one of diff[-1] is not added at all!
diff_threshold.append([diff[-1][1],diff[-1][0]]);


print "diff_threshold unsorted: ", diff_threshold[-1], "len: ", len(diff_threshold);


#diff_threshold is sorted in the order of index in loci array.
diff_threshold.sort();

showthreshold = open("./intermediate/showgreaterthreshold.sort.txt", 'w+');
for sh in diff_threshold:
	showthreshold.write("{}\n".format(sh));
showthreshold.close();

print "len(diff_threshold)",len(diff_threshold);	
#print diff_threshold;

#clusters fill with 
#[0]: CpG density.  #[1]: # of CpGs.  #[2]: length of CpG island region. 
#[3]: index of the first CpG of this island region in loci array.
#[4]: index of the last CpG of this island in loci array.
#[5]: the absolute start location of the first CpG in this island
#[6]: the absolute end location of this CpG island.
#[7]: the %G+C. #[8]: the number of G+C. 
#[9]: the energy density. #[10]: the total energy.

clusters = [];
last_loc = [0,0];

#diff_threshold contains the gap/distance/difference whose lengths are greater than threshold.
# each one's first element in diff_threshold holds the index of loci array. 
for loc in diff_threshold:
	# If current gap's index in loci array (loc[0]) is not the last gap's index + 1 in loci array (last_loc[0]+1),
	# 
	#if loc[0] > 13853 and loc[0] < 13859:
	#	print loc[0],loc[0]-1, loci[loc[0]-1], loci[last_loc[0]], "X ", last_loc[0], last_loc[0] + 1;
	if loc[0] != last_loc[0] + 1:
		clusters.append([float(loci[loc[0]-1]-loci[last_loc[0]]+2)/float(loc[0] -last_loc[0]), loc[0] -last_loc[0], \
		loci[loc[0]-1]-loci[last_loc[0]]+2, last_loc[0], loc[0]-1, loci[last_loc[0]], loci[loc[0]-1], 0 , 0 , 0, 0 ]);
		#clusters.append([loci[loc[0]-1]-loci[last_loc[0]], last_loc[0], loc[0]-1, loci[last_loc[0]], loci[loc[0]-1]]);
		#clusters.append([float(loci[sh[0]-1]-loci[last_sh[0]])/float(sh[0] -last_sh[0]), sh[0] -last_sh[0], loci[sh[0]-1]-loci[last_sh[0]], last_sh[0], sh[0]-1, loci[last_sh[0]], loci[sh[0]-1]]);
	last_loc = loc;

clusters.sort();
#clusters.reverse();

print "len(clusters)",len(clusters);

####################initialize the gs_matrix.
gs = [];		#read from gaosi.txt, same size to the filter's length
gs_matrix = [];
gsfile = open("gaussian.txt",'r');

for l in gsfile:
	gs.append(int(l));
	
gsfile.close();


tmp = [];
#initial the first row 
tmp.append(gs[0]);
for i in range(1,len(gs)):
	tmp.append(tmp[i-1]+gs[i]);
gs_matrix.append(tmp);

#after initialization of the first row, fill up the others of gs_matrix.
for i in range(1,len(gs)):	
	tmp = [];
	tmp.append(gs_matrix[0][i]);
	for j in range(1,len(gs)):
		tmp.append(tmp[j-1]+gs[j]);
	gs_matrix.append(tmp);


for m in gs_matrix:
	print m;
######################initialization end of Gaussian.

filter_GC = 0;
filter_len8 = 0;
filter_len50 = 0;
filter_0 = 0;
filter_0_val = float(256/float(2*len(gs)+2));		#filter for energy density. only 3000 left.
filter_den_cpg = 0;


refined_clusters = [];

refined_debug = open("refined_debug.txt", "w+");

#After calculating the clusters' raw CpGIslands, calculate the %G+C must greater than 50% for each island.
#And the energy density.
buf = "";
for i, cl in enumerate(clusters):
	#get the buf filled with the letters.
	if cl[2] >= min_region_len:
		buf = contents[cl[5]-1:cl[6]+1];
	else:
		if (min_region_len - cl[2]) < 0: 
			print "\n\nWRONG...............";
		#if the region is less than min_region_len's definition, fill it to length of (min_region_len).	
		buf = contents[(cl[5]- 1 - (min_region_len - cl[2] )/2): (cl[6] + (min_region_len - cl[2] )/2)];
	
	#each one count the %GC
	counter_gc = 0;
	for j in buf:
		if j == 'c' or j == 'C' or j == 'g' or j == 'G':
			counter_gc += 1;
	
	cl[7] = float( float(counter_gc) / float(len(buf)));
	cl[8] = counter_gc;
	
	#Gaussian function
	k_i = 0; k_j = 0;
	energy = 0;
	#find all CpGs in CGI.
	for k in range(cl[3],cl[4]+1):
		#k_i is row; k_j is column.
		if loci[k] + len(gs) - 1 <= loci[cl[4]]:
			k_j = len(gs)-1;
		else:
			k_j = loci[cl[4]] - loci[k];
		
		if loci[k] - (len(gs) - 1) >=  loci[cl[3]]:
			k_i = len(gs)-1;
		else:
			k_i = loci[k] - loci[cl[3]] + 1;
	
		if k_i > len(gs) or k_j > len(gs):
			print "out...", k, cl[3], cl[4];
		
		energy += gs_matrix[k_i][k_j];
	cl[9] = float(float(energy)/float(cl[2]));
	cl[10] = energy;	
	
	#after calcualting the %G+C and the energy density for each cluster, we need to inspect each cluster to determinate if its energy density is lower than
	# expectation. If so, it might be divided into subcluster with higher energy density.
	# It will calculate the energy at each point within this cluster. At each point its energy is affected by the neighboring CpGs..
	
	#cl[3]: index of the first CpG of this island region in loci array.
	#cl[4]: index of the last CpG of this island in loci array.
	#cl[5]: the absolute start location of the first CpG in this island
	#cl[6]: the absolute end location of this CpG island.	
	as0 = cl[5]- 1;	#absolute start of buf.
	b1 = cl[3]; b2 = cl[4]; #index for the cover scope of lower band and upper band.
	vbuf = [];

	if cl[9] < filter_0_val: #if less than the threshold, it must exist some inner-islands.
		vbuf = [0 for b in range(len(buf))];
		#calculate each one in vbuf.
		for b, v in enumerate(vbuf):
			e_provider = [];
			for potential in range(b1, b2+1):
				if loci[potential] < as0+b - (len(gs)-1):
					continue;
				elif loci[potential] >= as0+b - (len(gs)-1) and loci[potential] <= as0+b+(len(gs)-1):
					e_provider.append(potential);
				elif loci[potential] > as0+b + (len(gs)-1):
					break;

			#if i == 21239:
			#	print b;			
			for potential in e_provider:
				if 	abs(loci[potential]-(as0+b)) >= len(gs):
					continue;
				vbuf[b] += gs[abs(loci[potential]-(as0+b))];
		
		#refined_debug.write("{} {}\n {} \n {} \n".format(i, vbuf, buf, cl));
		if i == 21239:
			refined_debug.write("{} {}\n {} \n {} \n".format(i, vbuf, buf, cl));	
		#scan vbuf, for those cluster whose value above filter_0_val, put it to clusters.			
		
		
		
		
		#this one is trying to re-cluster this cluster.
		sb = 0 ; eb = 0; nb = 0; #sb:start point, eb: end point, nb: new one.
		for b, v in enumerate(vbuf):
			#core should be powerful. 2*threshold
			if v >= 2*filter_0_val and nb == 0:
				sb = b;
				eb = b;
				nb = 1;
			elif v < 2*filter_0_val and nb == 1:
				eb = b-1;
				if eb - sb > min_region_len:
					refined_clusters.append([as0+sb, as0+eb, cl[3],cl[4], cl[5], cl[6]]);
					#if continues a subcluster. it is assigned to the threshold.
					cl[9] = filter_0_val;
				nb = 0;
				sb = 0;
				eb = 0;
			elif v >= filter_0_val and nb == 1:
				eb = b;

		if nb == 1 and eb-sb > min_region_len:
			refined_clusters.append([as0+sb, as0+eb, cl[3],cl[4], cl[5], cl[6]]);	
			
			
			#	if loci[b1] + len(gs)-1 < as+b
			#		b1 += 1;
			#		continue;
			#	elif loci[b1] > as+b:
			#		break;
			#	elif loci[p] + len(gs)-1 >= as+b and loci[p] < as+b:
			#		v += gs[loci[p]-(as+b)];
			
			#for k in range(b1, cl[4]):
			#	if loci[b1] - (len(gs)-1) <= (as+b):
			#		b2 += 1; 

refined_debug.close();
	
print "len(refined_clusters) -> ", len(refined_clusters);	

refinedfile = open("refined.txt", 'w+');		 
for i in refined_clusters:
	refinedfile.write("{}\n".format(i));
refinedfile.close();


clustersfile = open("./intermediate/{}.init.txt".format(clusterfilename), 'w+');

clustersfile.write("[density, CpG#, span, start_ind, end_ind, start_loc, end_loc]\n");
for cl in clusters:
	#if the %G+C is less than 50%, drop it off.
	#Chr21 total #56***.
	#if cl[7] <= 0.5:	#27xxx.
	#	continue;
	#if cl[2] < min_region_len: #26xxx.
	#	continue;
	clustersfile.write("{}\n".format(cl));
	#The following is to output the CpG Island only.
	#clustersfile.write("{} {}\n".format(cl[5],cl[6]));



#statistic clusters is for counting the number of cpg islands in (clusters) which have the same density.
#last_cl[0] is the density of cpg. [1] is the counter.
last_cl = [0,1];
cluster_stat = [];
for cl in clusters:
	if int(cl[0]) != int(last_cl[0]): 
		cluster_stat.append([int(last_cl[0]), last_cl[1]]);
		last_cl = [int(cl[0]), 1];
	else:
		last_cl[1] += 1;
cluster_stat.append([int(last_cl[0]), last_cl[1]]);



#write the statistics to the file.
for st in cluster_stat:	
	#print st[0],"\t",st[1];	
	clustersfile.write("{}\t{}\n".format(st[0],st[1]));
	
clustersfile.close();



##################################Begin entropy is blow for things in each cluster (raw CGI).
human_di_prob = {'CG':1, 'GC':4.3, 'CC':4.7, 'GT':4.9, 'GG':5.0, 'AC':5.4, 'TC':5.7, 'GA':6.1,\
							'TA':6.7,'AG':7.0,'CT':7.1, 'TG':7.4,'CA':7.4, 'AT':8.1, 'AA':9.7,'TT':9.7};
hu_di_I = {'CG':0, 'GC':0, 'CC':0, 'GT':0, 'GG':0, 'AC':0, 'TC':0, 'GA':0,\
							'TA':0,'AG':0,'CT':0, 'TG':0,'CA':0, 'AT':0, 'AA':0,'TT':0};
human_E = 0;							
for di in human_di_prob.keys():
	hu_di_I[di] = -math.log(float(human_di_prob[di])/100,2);
	print di, hu_di_I[di];
	human_E += float(human_di_prob[di])*hu_di_I[di]/100;
print human_E;

buf = "";

CGI_len_threshold = 2*len(gs);
gap_len_threshold = 6*len(gs);

tmp_dbg_file = open("./intermediate/cluster_sub.txt", 'w+');
print "len(gs)",len(gs);
subclusters = [];
for i, cl in enumerate(clusters):
	if cl[2] <= CGI_len_threshold:
		continue;
	#get the buf filled with the letters.
	buf = contents[cl[5]-1:cl[6]+1];
	tmp_dbg_file.write("{}\n".format(cl));
	tmp_dbg_file.write("{}\n".format(buf));
	#subclusters = [];
	last_k = cl[3];
	start_k = cl[3];
	for k in range(cl[3]+1,cl[4]+1):
		buf = contents[loci[last_k]-1:loci[k]+1];
			 
		#if loci[k] - loci[last_k] > gap_len_threshold or ('a' in buf or 'c' in buf or 'g' in buf or 't' in buf):
		if loci[k] - loci[last_k] > gap_len_threshold :
			if start_k != last_k:
				buf = contents[loci[last_k]-1:loci[k]+1];
				if 'a' in buf or 'c' in buf or 'g' in buf or 't' in buf:
					if (loci[last_k] - loci[start_k]) > CGI_len_threshold:
						subclusters.append([float((loci[last_k] - loci[start_k]))/last_k-start_k+1,last_k-start_k,(loci[last_k] - loci[start_k]),start_k, last_k,loci[start_k],loci[last_k]]);
						buf = contents[loci[start_k]-1:loci[last_k]+1];
						tmp_dbg_file.write("\t{}_{} k={}- {} \t{}\n".format(loci[start_k],loci[last_k], k, loci[k], buf));
					else:
						start_k = k;
				else:
					di_num = {'CG':0, 'GC':0, 'CC':0, 'GT':0, 'GG':0, 'AC':0, 'TC':0, 'GA':0,\
							'TA':0,'AG':0,'CT':0, 'TG':0,'CA':0, 'AT':0, 'AA':0,'TT':0};
					this_E = 0;
					last_b = 0;
					twin = '';
					for b in range(1,len(buf)):
						twin = '{}{}'.format(buf[last_b],buf[b]);
						di_num[twin] += 1;
						last_b = b;

					for di in di_num.keys():
						this_E += float(di_num[di])*hu_di_I[di]/(len(buf)-1);
						
					if this_E < human_E:
						if (loci[last_k] - loci[start_k] )> CGI_len_threshold:
							subclusters.append([float((loci[last_k] - loci[start_k]))/last_k-start_k+1,last_k-start_k,(loci[last_k] - loci[start_k]),start_k, last_k,loci[start_k],loci[last_k]]);
							buf = contents[loci[start_k]-1:loci[last_k]+1];
							tmp_dbg_file.write("\tcluster @ {}_{},  discard gap {}_{}, cluster: \t{} \n".format(loci[start_k],loci[last_k], loci[last_k], loci[k], buf));
						else:
							#if length of CGI is less than expected,
							start_k = k;
						buf = contents[loci[last_k]-1:loci[k]+1];
						tmp_dbg_file.write("\t\t gap @ {}_{}  this_E {} \t{}\n".format( loci[last_k], loci[k], this_E, buf));
					else: #if this entropy is greater than human's entropy, nothing to do. even though
					# the gap is greater than 2*len(gs), it should keep running.start_k should not change.
						tmp_dbg_file.write("\t\t contain gap{}_{} because this_E {} good \t{}\n".format(loci[last_k], loci[k], this_E, buf));
						last_k = k;
						continue;
			start_k = k;
		#else loci[k] - loci[last_k] < 2*len(gs):			
			
		elif k == cl[4]:
			if start_k != k:
				if (loci[k] - loci[start_k] ) > CGI_len_threshold:
					subclusters.append([float((loci[k] - loci[start_k]))/k-start_k+1,k-start_k,(loci[k] - loci[start_k]),start_k, k,loci[start_k],loci[k]]);
					buf = contents[loci[start_k]-1:loci[k]+1];
					tmp_dbg_file.write("\tlast {}_{} :\t{}\n".format(loci[start_k],loci[k],buf));
			start_k = k;			
		last_k = k;
			
tmp_dbg_file.close();

subclusters.sort();

print len(subclusters);

##################################End entropy  ###############################



#clusters fill with 
#[0]: the energy density.   #[1]: CpG density.  #[2]: the %G+C.
#[3]: index of the first CpG of this island region in loci array.
#[4]: index of the last CpG of this island in loci array.
#[5]: the absolute start location of the first CpG in this island
#[6]: the absolute end location of this CpG island.
#[7]: number of CpGs. #[8]: #length of CpG island region.  

resultfile = open("./intermediate/{}.result.txt".format(clusterfilename), 'w+');
result_clusters = [];
clusters_len = [];

for cl in clusters:
	if cl[7] <= 0.5:	#27xxx.
		filter_GC += 1;
		continue;
	if cl[2] < 4 * min_region_len: #26xxx. filter 5*8=40
		filter_len8 += 1;
		continue;
	if 	float(cl[9])/float(1) < 1.8*filter_0_val:
		filter_0 += 1;
		continue;
	if cl[0] > 27:
		filter_den_cpg += 1;
		continue;
	result_clusters.append([cl[9], cl[0], cl[7], cl[3], cl[4], cl[5], cl[6],cl[1], cl[2], cl[8], cl[10]] )
	clusters_len.append([cl[2],cl[1], cl[9], cl[0], cl[7], cl[3], cl[4], cl[5], cl[6], cl[8], cl[10]]);

#result_clusters is applied by all sorts of filters.	
result_clusters.sort();
result_clusters.reverse();

#clusters_len is used for display of length.
clusters_len.sort();
clusters_len.reverse();

filter_repeat = 0;

last_res = result_clusters[0];
for res in result_clusters:	
	if last_res[0] == res[0] and last_res[1] == res[1] and last_res[2] == res[2]:
		filter_repeat += 1;
	#print st[0],"\t",st[1];	
	resultfile.write("{}\n".format(res));
	last_res = res;

resultfile.write("========================================================================\n\n\n");
	
for len_c in clusters_len:
	if len_c[0] < 50:
		filter_len50 += 1;
	resultfile.write("{}\n".format(len_c));

resultfile.close();


print "filter_GC=", filter_GC, "filter_0_val = ", filter_0_val, " filter_0=", filter_0, "filter_len8", filter_len8, "filter_den_cpg", filter_den_cpg, "filter_len50", filter_len50;
print "len(result_clusters)", len(result_clusters), " - repeat", filter_repeat-1, "  = ", len(result_clusters)-filter_repeat+1;

final_cgi = [];
for fnl in subclusters:
#for fnl in result_clusters:
	final_cgi.append([fnl[5],fnl[6]]);
final_cgi.sort();




finalfile = open("./final/{}".format(clusterfilename), 'w+');
for fnl in final_cgi:
	finalfile.write("{} {}\n".format(fnl[0],fnl[1]));
finalfile.close();

print "write in",len(result_clusters);


