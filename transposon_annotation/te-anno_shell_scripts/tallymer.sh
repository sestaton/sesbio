#!/bin/bash

#-------------------------
# build the suffix array
#-------------------------
gt suffixerator -dna -pl -tis -suf -lcp -v -parts 4 -db read1.fna read2.fna -indexname reads

(...)

#------------------------------------
# call occratio; get range of k-mers
#------------------------------------
gt tallymer occratio -output unique nonunique -minmersize 10 -maxmersize 20 -esa reads

# distribution of unique mers
10 223755
11 373775
12 444083
13 465859
14 468735
15 465646
16 460791
17 455449
18 450049
19 444720
20 439532
# distribution of non unique mers (counting each non unique mer only once)
10 135526
11 92162
12 62347
13 49611
14 44618
15 42301
16 40867
17 39769
18 38815
19 37961
20 37166
# space peak in megabytes: 0.06
# mmap space peak in megabytes: 3.93

#------------------------------------------------
# call occratio; get ratio relative to the total
#------------------------------------------------
gt tallymer occratio -output unique relative -minmersize 10 -maxmersize 20 -esa reads

# distribution of unique mers
10 223755 0.623
11 373775 0.802
12 444083 0.877
13 465859 0.904
14 468735 0.913
15 465646 0.917
16 460791 0.919
17 455449 0.920
18 450049 0.921
19 444720 0.921
20 439532 0.922
# space peak in megabytes: 0.06
# mmap space peak in megabytes: 3.93

#---------------------------------------------
# -mersizes option restricts the calculations
#---------------------------------------------
gt tallymer occratio -output unique nonunique -mersizes 10 13 17 -esa reads

# distribution of unique mers
10 223755
13 465859
17 455449
# distribution of non unique mers (counting each non unique mer only once)
10 135526
13 49611
17 39769
# space peak in megabytes: 0.06
# mmap space peak in megabytes: 3.93

#------------------------------------------------------------------------------
# while occratio runs for a range of mer sizes, mkindex runs on a fixed length
#------------------------------------------------------------------------------
gt tallymer mkindex -mersize 19 -minocc 40 -esa reads

1 444720
2 30886
3 3909
4 1397
5 640
6 335
7 111
8 172
(...)

# construct mer buckets for prefixlength 4
# numofcodes = 256
# indexfilename = tyr-reads
# alphasize = 4
# mersize = 19
# numofmers = 3166
# merbytes = 5
# space peak in megabytes: 0.06
# mmap space peak in megabytes: 3.93

#-------------------------------------------------------------------------------------------------------------
# The program search now uses the index reads and matches all 19-mers of the input sequence U89959 against it
#-------------------------------------------------------------------------------------------------------------
gt tallymer search -output qseqnum qpos counts sequence -tyr tyr-reads -q U89959.fna

0 +5966 4 tcttcttcttcttcttctt
0 +17269 14 atatatatatatatatata
0 +17270 12 tatatatatatatatatat
0 +17271 14 atatatatatatatatata
0 +17272 12 tatatatatatatatatat
0 +71281 6 tcatcatcatcatcatcat
(...)