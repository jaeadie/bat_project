output = open('primer3inputs.fasta','w')


def primer3(header,template, start, length):
  output.write("".join(["SEQUENCE_ID=",header]))
  output.write("\n")
  output.write("".join(["SEQUENCE_TEMPLATE=",template]))
  output.write("\n")
  output.write("".join(["SEQUENCE_TARGET=",start,",",length]))
  output.write("\n")
  output.write("PRIMER_TASK=pick_detection_primers")
  output.write("\n")
  output.write("PRIMER_PICK_LEFT_PRIMER=1")
  output.write("\n")
  output.write("PRIMER_PICK_INTERNAL_OLIGO=1")
  output.write("\n")
  output.write("PRIMER_PICK_RIGHT_PRIMER=12")
  output.write("\n")
  output.write("PRIMER_OPT_SIZE=18")
  output.write("\n")
  output.write("PRIMER_MIN_SIZE=15")
  output.write("\n")
  output.write("PRIMER_MAX_SIZE=21")
  output.write("\n")
  output.write("PRIMER_MAX_NS_ACCEPTED=1")
  output.write("\n")
  output.write("PRIMER_PRODUCT_SIZE_RANGE=110-200")
  output.write("\n")
  output.write("P3_FILE_FLAG=1")
  output.write("\n")
  output.write("SEQUENCE_INTERNAL_EXCLUDED_REGION=37,21")
  output.write("\n")
  output.write("PRIMER_EXPLAIN_FLAG=1")
  output.write("\n")
  output.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=/homes/dmamartin-x/software/primer3-2.3.6/src/primer3_config/")
  output.write("\n")
  output.write("=")
  output.write("\n")


import subprocess

fh = open('pipymmfp.txt')
#fh2 = open('pipymmfpoutput.txt', 'w')
#fh3 = open('pipycontigs.txt', 'w')
#fh4 = open('pipyssrstart.txt', 'w')
#fh5 = open('pipyssrend.txt', 'w')
#fh6 = open('pipyssrlen.txt', 'w')

name = 0
start = 3
end = 4
motif = 1
seq = 9
seq_len = 5
array = []
lcname=''
lcend=0
contigs_checked=0
candidates=0
array3 = []

for line in fh.readlines():
  array2 = line.split("\t")
  if array2[name]=='Seq_Name':
      continue
  contigs_checked += 1
  if lcname==array2[name]: 
     if int(array2[start])-lcend < 130:
         lcend=int(array2[end])
         continue
  if (len(array2[motif])>3 and len(array2[seq])>89 and len(array2[seq]) <151):

     temp_seq="pipy:"+array2[name][0:14]
     subp=subprocess.Popen(['seqret',temp_seq, '-auto', '-stdout'], stdout=subprocess.PIPE)
     (stdout, stderr)=subp.communicate()
     tseq=''.join(stdout.split("\n")[1:])
     
     primer3(array2[name],tseq, array2[start], array2[seq_len])
     
     #array.append(line)
     #array3.append(line[0:14])	
     #fh2.write(line)
     #fh3.write("".join(["Pipy:",line[0:14],"[",str(array2[start]),":",str(array2[seq_len]),"]"]))
     #fh3.write("\n")	
     #fh6.write("".join([str(array2[start]),",",array2[seq_len]]))
     #fh6.write("\n")    
     candidates += 1
  lcname=array2[name]
  lcend=int(array2[end])
#fh2.close()  
print("%s contigs checked. %s candidates retained\n"%(contigs_checked, candidates))
