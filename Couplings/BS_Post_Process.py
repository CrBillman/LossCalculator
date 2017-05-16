N_ATOMS=1008

BISC=open("bisc.out","r")
PPOUT=open("pp-bisc.out","wb")

dup_cut=0.0005
dup_dist_cut=0.01

BISC.readline()
Good_TLS=[]
Index=[]
count=0
while True:
	flag=None
	line=BISC.readline()
	if( not line):
		break
	temp=line.split()
	if((temp[4]=="1") or (temp[4]=="5")):
		if((float(temp[5])>0.0) and (float(temp[6])>0.0)):
			for i in xrange(0,len(Good_TLS)):
				TLS_Comp=Good_TLS[i]
				dist_cut=dup_dist_cut*float(temp[3])
				asym_diff=abs(float(temp[2])-float(TLS_Comp[2]))
				dist_diff=abs(float(temp[3])-float(TLS_Comp[3]))
				if((asym_diff<dup_cut) and (dist_diff<dist_cut)):
					flag="Duplicate"
					break
			if(flag!="Duplicate"):
				Good_TLS.append(temp)
				PPOUT.write(line)
				Index.append(temp[1])
	count=count+1
PPOUT.close()
BISC.close()

MIN=open("min_saddle","r")
MOUT=open("pp-min_saddle","wb")

for i in Index:
	lookout_str="         " + str(i) + "\n"
	if(len(lookout_str)>11):
		lookout_str=lookout_str[len(lookout_str)-11:]
	
	while True:
		line=MIN.readline()
		if(line==lookout_str):
			MOUT.write(line)
			for j in xrange(0,N_ATOMS*4):
				MOUT.write(MIN.readline())
			break
