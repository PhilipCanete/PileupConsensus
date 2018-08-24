# PileupConsensus

Generates consensus bases from pileup file as argument, using the following steps:

1. Remove non-mismatch characters and indels (^. $ . , +/-AGCT), updating coverage as appropriate
2. Remove bases that do not have a Phred quality over 40, updating coverage as appropriate
3. Default to reference base if coverage is below 5 or above 100
4. Accept mismatch if  percentage of most common mismatch is over 80% of coverage, else default to reference base
5. Write to a file whose name is the pileup argument file with a -results.txt extension
