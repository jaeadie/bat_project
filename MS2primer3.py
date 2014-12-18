# open output file for primer3 inputs
output = open('primer3inputs.fasta','w')

# the function to write each entry into the primer3 file with all required fields. 
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

# import the subprocess module for use. 
import subprocess

# open the input microsatellite data file and assign it a file handle
fh = open('pipymmfp.txt')

# assign each tab of the input file to what its contents are, eg. the start position of the microsatellite being tab 3.
name = 0
start = 3
end = 4
motif = 1
seq = 9
seq_len = 5
# create an empty array 
array = []
# set the last contig checkeds name.
lcname=''
#set the end position of th elast contig checked
lcend=0
# the number of contigs that have been checked
contigs_checked=0
# the number of candidates returned
candidates=0
#create another empty array
array3 = []
# read and loop over each line in the input file
for line in fh.readlines():
# for each line of the input file split it by the tabs of the table and store in an array.
  array2 = line.split("\t")
#if the current lines name starts with Seq_name continue to the next line
  if array2[name]=='Seq_Name':
      continue
    # add one to the number of contigs checked.
  contigs_checked += 1
  #check to see if the last contigs name is the same as the current one
  if lcname==array2[name]: 
    # if the contigs are the same, check to see if the current microsatellite starts within 130bp of the previous one
     if int(array2[start])-lcend < 130:
       # change the last contig end to the one that has just been read.
         lcend=int(array2[end])
         continue
       # if the length of the motif tab is >3 and the sequences length is between 90 and 150bp then continue. 
  if (len(array2[motif])>3 and len(array2[seq])>89 and len(array2[seq]) <151):
# set the current line to be temp_seq inputting the first 14 letters of the contig name
     temp_seq="pipy:"+array2[name][0:14]
     #perform a subprocess to feed the output of temp_seq to retrieve the contig of the current candidate marker.
     subp=subprocess.Popen(['seqret',temp_seq, '-auto', '-stdout'], stdout=subprocess.PIPE)
    # communicate the subporocess and perform its function with a standard output file.
     (stdout, stderr)=subp.communicate()
    # assign tseq as the contig sequence from the stdout. splitting the input by each line and skipping the first line, which will contain the contig name rather than sequence. 
     tseq=''.join(stdout.split("\n")[1:])
    #perform the primer3 function defined earlier with each parameter set for writting to the correct line. 
     primer3(array2[name],tseq, array2[start], array2[seq_len])
     
# add one to the candidate markers
     candidates += 1
     #set the last contigs name to the current lines
  lcname=array2[name]
  # set the last contigs end postition to the current lines
  lcend=int(array2[end])
 # print the number of contigs that have been checked and how many candidates there were
print("%s contigs checked. %s candidates retained\n"%(contigs_checked, candidates))
