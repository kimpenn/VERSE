VERSE supports the following modes of RNA-Seq quantification:
    	1. Original featureCounts (Default)
    	2. HTSeq Union (-z 1)
    	3. HTSeq Intersection-strict (-z 2)
    	4. HTSeq Intersection-nonempty (-z 3)
    	5. VERSE Union-strict (-z 4)
	6. VERSE Cover-length (-z 5)

Supported Quantification Schemes:
	Hierarchical Assign: assign reads to feature types according to their priority.
	Independent Assign: assign reads to feature types independently in a single run.

*** Installation *** 

	You can type "make" to see instructions.
	For Linux OS:
		make -f Makefile.Linux
	For Mac OS, use command:
		make -f Makefile.MacOS

*** Usage ***

Please run "./verse" to see the details.

A sample command:
./verse -a testdata/test.gtf -t 'exon' -g gene_id -z 3 -s 1 -o testdata/intersection_nonempty.stranded.paired testdata/PE.sam

A sample hierarchical assign command:
./verse -a testdata/test.gtf -t 'exon;intron;xine' -g gene_id -z 3 -o testdata/intersection_nonempty.unstranded.paired.hierarchical testdata/PE.sam

***************

VERSE is developed based on the framework of featureCounts(SUBREAD), which is written by Drs Yang Liao and Wei Shi.

This work is supervised and generously supported by Dr. Stephen Fisher and Professor Junhyong Kim.

Please contact Qin Zhu if you have any questions.

***************


by Qin Zhu
Junhyong Kim Lab
University of Pennsylvania
2015

