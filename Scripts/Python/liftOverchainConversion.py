#! /usr/bin/env python

import sys, getopt

def main(argv):
   inputfile = ''
   ref1 = ''
   ref2 = ''
   outputfile = ''

   try:
      opts, args = getopt.getopt(argv,"i:1:2:o:",["in=","ref1=","ref2=","out="])
   except getopt.GetoptError:
      print('liftOverchainConversion.py -i <inputfile> -o <outputfile> -1 <ref1> -2 <ref2>');
      sys.exit(2);
   for opt, arg in opts:
      if opt == '-h':
         print('liftOverchainConversion.py -i <inputfile> -o <outputfile> -1 <ref1> -2 <ref2>');
         sys.exit();
      elif opt in ("-i", "--in"):
         inputfile = arg;
      elif opt in ("-1", "--ref1"):
         ref1 = arg;
      elif opt in ("-2", "--ref2"):
         ref2 = arg;
      elif opt in ("-o", "--out"):
         outputfile = arg;


   dict_r1 = {}
   with open(ref1,'r') as f:
      for line in f:
         (val,key) = line.split()
         dict_r1[key] = val

   dict_r2 = {}
   with open(ref2,'r') as f:
      for line in f:
         (val,key) = line.split()
         dict_r2[key] = val


   fout = open(outputfile,'w')

   with open(inputfile,'r') as f:
      for line in f:
         if line.startswith('chain'):
            line_split = line.split();
            chr_1 = dict_r1[line_split[2]]
            chr_2 = dict_r2[line_split[7]]
            line_split[2] = chr_1
            line_split[7] = chr_2
            line_joined = ' ' .join(line_split)
            fout.write(line_joined)
         else:
            fout.write(line);


   fout.close();

if __name__ == "__main__":
   main(sys.argv[1:])

