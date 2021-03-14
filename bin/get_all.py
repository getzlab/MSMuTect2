#Pyhton managing file


def EM_ML (vec,num_allele,P) :
	#vec is a two dimenssional vectro continaing the histogram


	#The search for the alleles will be only from lenghts that are larger than 5 reads.
	#This mean that I need to have two lists of availbile lenghts. All of them, and only those that are larger than 5 reads.
	#vec=vec[:,np.nonzero(vec[0,:]<40 and np.nonzero(vec[0,:]>4)];
	#vec=vec[:,np.nonzero(vec[0,:]>4 )];
	exsiting_reads=vec[0,:];
	exsiting_reads_num=vec[1,:];

	Large_reads=exsiting_reads[np.where(exsiting_reads_num>5)[0]];
	Large_reads_location=np.where(exsiting_reads_num>5)[0];
	#print "Large_Reads",Large_reads, "num allele" , num_allele;

	num_allele_v=np.arange(num_allele);


	L_max=-1e9;

	for attempt in np.arange(10):

		Q_Theta_j=0*np.ones(num_allele);
		New_theta=0*np.ones(num_allele);
		Theta_new_temp=0*np.ones(Large_reads_location.size);

		if num_allele<1:
		#Fit only for one allele viw brute force

			return 1;
		else:
			A=0*np.ones(num_allele);
			f=0*np.ones(num_allele-1);
			f_new=0*np.ones(num_allele);

			#initial values for the maximization process

			#A=Large_reads[np.random.randint(0,Large_reads.size,num_allele)];
			A=Large_reads[np.random.permutation(Large_reads.size)[0:num_allele]]
			#print "initial A " ,A
			f=np.random.rand(num_allele-1);
			f=np.append(f,1-sum(f[:num_allele-1]));
			f=np.ones(num_allele)*(1.0/num_allele);
			Z_i_j=0*np.ones([44,num_allele])
			Q_Theta_j=-1e9*np.ones(num_allele);
			new_theta=0*np.ones(num_allele);

			change=1e6;
			prevL=1e6;
			it=0;
			while change>1e-5:
			#for it in np.arange(10):
				it=it+1;
				#print exsiting_reads," bla ",num_allele_v
				#print P;
				#print A;
				#print f;
				#print Large_reads ," and ",Large_reads_location;

				#Step number 0 for each iteration make the Z_i_j matrix
				for j in num_allele_v:
					for i in exsiting_reads:
					    Z_i_j[i,j]=P[A[j],i]*f[j]/sum(P[A[:],i]*f[:]+1e-10);
				#print Z_i_j;
				#print "old f ", f
				#Step number 1. From the Z_i_j's estimate the new frequencies.
				for j in num_allele_v:
					f_new[j]=sum(Z_i_j[exsiting_reads,j]*exsiting_reads_num)/sum(exsiting_reads_num);
					#print Z_i_j[exsiting_reads,j],exsiting_reads_num , sum(Z_i_j[exsiting_reads,j]*exsiting_reads_num) , sum(exsiting_reads_num)

				#Step number 2. Maximize the new Thetas
				for j in num_allele_v:
					for k in np.arange(Large_reads.size):
						Test_theta=Large_reads[k];
						#for i in exsiting_reads:
						Theta_new_temp[k]=sum(Z_i_j[exsiting_reads,j]*np.log(P[Test_theta,exsiting_reads]+1e-10)*exsiting_reads_num)
					New_theta[j]=Large_reads[Theta_new_temp.argmax()];
					#print Theta_new_temp;

				for j in num_allele_v:
					A[j]=New_theta[j];
					f[j]=f_new[j];

				Q_theta=0;
				Q_f=0;
				for j in num_allele_v:
					Q_theta=Q_theta+sum(Z_i_j[exsiting_reads,j]*np.log(P[A[j],exsiting_reads]+1e-10)*exsiting_reads_num)
					Q_f=Q_f+sum(Z_i_j[exsiting_reads,j]*np.log(f[j]+1e-10)*exsiting_reads_num);


				#Calcualte the likelihood
				#L(Theta|D)=P(D|Theta)=PI_i(d_i|Theta)=PI_i(SUM_j(f[j]*p(d_i|theta_j))).
				#In our case we can combine all the reads with the same number of repeats thus:
				L_k_log=0;
				for k in np.arange(exsiting_reads.size):
					L_k_log=L_k_log+exsiting_reads_num[k]*np.log(sum(f*P[A,exsiting_reads[k]])+1e-10 )
					#print 	f,P[A,exsiting_reads[k]],exsiting_reads_num[k],"Log L",L_k_log;
				#print "New_theta  is ", New_theta , " f_new is " , f_new, " Q_theta is ", Q_theta, " Q_f is ", Q_f, " Q is ", Q_theta+Q_f ,"Log L =" , L_k_log;
				if L_k_log>L_max:
					L_max=L_k_log;
					A_best=A;
					f_best=f;
				change=np.abs(prevL-L_k_log);
				prevL=L_k_log;
				#print "New_theta  is ", New_theta , " f_new is " , f_new, " Q_theta is ", Q_theta, " Q_f is ", Q_f, " Q is ", Q_theta+Q_f ,"Log L =" , L_k_log, "cha" ,change;


			#print L_max ,A_best,f_best, A, f;

	A_res=[None]*num_allele;
	f_res=[None]*num_allele;
	for i in np.arange(num_allele):
		A_res[i]=A_best[i];
		f_res[i]=f_best[i];
	return L_max ,A_res,f_res

def Finding_the_allele(vec,P):



	#print vec;
	s=vec[0,:].size

	exsiting_reads=vec[0,:];
	exsiting_reads_num=vec[1,:];



	Large_reads=exsiting_reads[np.where(exsiting_reads_num>5)[0]];
	Large_reads_location=np.where(exsiting_reads_num>5)[0];

	[L1,A1,f1]=EM_ML(vec,1,P);
	#print L1,A1,f1


	[L2,A2,f2]=EM_ML(vec,2,P);
	#print L2,A2,f2

	D12=-2*L1+2*L2;
	#print L2,A2,f2,D12
	res1=L1;res2=A1;res3=f1;
	if D12>0:
		P12_value= stats.chi2.pdf(D12,2);
		#print L2,A2,f2,D12,P12_value
		if   P12_value>0.05:
			res1=L1;res2=A1;res3=f1;
		elif Large_reads.size==2:
			res1=L2;res2=A2;res3=f2;
		else:
			[L3,A3,f3]=EM_ML(vec,3,P);
			#print L3,A3,f3
			D23=-2*L2+2*L3;
			P23_value= stats.chi2.pdf(D23,2);
			#print L2,A2,f2,D23,P23_value
			if   P23_value>0.05:
				res1=L2;res2=A2;res3=f2;
			elif  Large_reads.size==3:
				res1=L3;res2=A3;res3=f3;
			else:
				[L4,A4,f4]=EM_ML(vec,4,P);

				D34=-2*L3+2*L4;
				P34_value= stats.chi2.pdf(D34,2);
				#print L3,A3,f3,D34,P34_value
				if   P34_value>0.05:
					res1=L3;res2=A3;res3=f3;
				else:
					res1=L4;res2=A4;res3=f4;
	np.append(res1,res2);
	return [res1,res2[:],res3[:]];



import csv
import time
#import scipy as sci
import scipy.stats as stats
import numpy as np
import sys

start = time.time();
#Should upload each of the error tables
P=np.loadtxt(sys.argv[2], delimiter=',');

nl = "\n";

f = open(sys.argv[1], "r")
strin=sys.argv[1] + ".all"
#print strin;
f_w = open(strin, "w")
sto=1;
#for temp in fileinput.input(f):
count=0;
for temp in f:
#temp=f.readline();
#while sto==1:
	count=count+1;
	#print count;
	temp1=temp.split("\n")[0];
	temp2=temp1.split(", ");

	#print temp2;
	num_val=len(temp2);
	for i in np.arange(1,num_val):
		temp2[i]=int(temp2[i]);
	#del vec;
	#print temp2;
	vec=np.array([temp2[1:num_val:2],temp2[2:num_val:2]]);
	fi=np.nonzero(vec[0,:]<40)
	vec=vec[:,fi[0]]
	fi=np.nonzero(vec[0,:]>4)
	vec=vec[:,fi[0]]

	#print vec;
	exsiting_reads=vec[0,:];
	exsiting_reads_num=vec[1,:];


	Large_reads=exsiting_reads[np.where(exsiting_reads_num>5)[0]];

	#print Large_reads
	if Large_reads.size==1:
		A=[[0],Large_reads,[1]];
		#print A;
		#line=np.append(vec[0],A[0]);
		#line=np.append(line,A[1]);
		#line=np.append(line,A[2]);
		#print "vec[0] is", type(vec[0].tolist), "A[0] is", type(A[0]),"A[1] is", type(A[1]),"A[2] is", type(A[2]);
		#line=temp2+A[0]+A[1]+A[2]+[nl];
		#print type(line)
		#print line.size

		for pri in temp2:
			f_w.write("%s " % pri)
		f_w.write("-999 " )
		f_w.write("%f " % A[0][0])
		f_w.write("-999 " )
		for item in A[1]:
			f_w.write("%d " % item)
		f_w.write("-999 " )
		for item in A[2]:
			f_w.write("%f " % item)
		f_w.write("\n")
	elif Large_reads.size>1:
		A=Finding_the_allele(vec,P);
	#print A;
		#line=np.append(vec[0],A[0]);
		#line=np.append(line,A[1]);
		#line=np.append(line,A[2]);
		#print "vec[0] is", type(vec[0].tolist), "A[0] is", type(A[0]),"A[1] is", type(A[1]),"A[2] is", type(A[2]);
		#line=temp2+A[0]+A[1]+A[2]+[nl];
		#print type(line)
		#print line.size


		for pri in temp2:
			f_w.write("%s " % pri)
		f_w.write("-999 " )
		f_w.write("%f " % A[0])
		f_w.write("-999 " )
		for item in A[1]:
			f_w.write("%d " % item)
		f_w.write("-999 " )
		for item in A[2]:
			f_w.write("%f " % item)
		f_w.write("\n")


	#temp=f.readline();
	#if temp==None:
	#	sto=1;

end = time.time()
#print end-start
