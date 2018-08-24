import sys
import os
import re
import collections

def main():
    filepath = sys.argv[1]
    if not os.path.isfile(filepath):                                                                                     # File error handling
        print("File path {} does not exist. Exiting...".format(filepath))
        sys.exit()

    with open(filepath) as file:
        f = open(sys.argv[1] + ".results", "w+")
        for lineNumber, line in enumerate(file):
            lineDict = {}

            for index, key in enumerate(("RefChromosome", "RefPosition", "RefBase", "Coverage", "ReadBases", "ReadQuality")):   # Dictionary keys
                lineDict[key] = line.strip().split('\t')[index]                                                         # Values are lists of strings
            lineDict["Coverage"] = int(lineDict["Coverage"])

            mismatch_chars_only(lineDict)
            if len(lineDict["ReadQuality"]) > 0:
                    check_quality(lineDict)
                    if len(lineDict["ReadQuality"]) > 0:
                         if 100 > lineDict["Coverage"] >= 5:
                             if check_mismatch_percent(lineDict):
                                 lineDict["RefBase"] = check_mismatch_percent(lineDict)
            write_consensus(lineDict, lineDict["RefBase"],f)

    if f is not None:
        f.close()

def write_consensus(lineDict, consensusBase,f):
    f.write("\nReference Chromosome: " + lineDict["RefChromosome"])
    f.write("\nReference Position: " + lineDict["RefPosition"])
    f.write("\nConsensus: " + consensusBase)
    f.write("\n\n")


def mismatch_chars_only(lineDict):
    lineDict["ReadBases"] = re.sub("(\^.)|\$","",lineDict["ReadBases"])                                                 # Remove [^.] mapping quality\

    i = 0
    replacements = []
    listQuality = list(lineDict["ReadQuality"])
    for m in re.finditer("([+\-]\d+[ACGT]+)", lineDict["ReadBases"]):                                                   # Iterate over indel regex matches
        replacements.append(m.start()-i)
        i = m.end()-m.start() + i - 1                                                                                   # i is the offset between the index of ReadBases corresponding ReadQuality
        lineDict["Coverage"] -= 1
    for index in reversed(replacements):
        del listQuality[index]                                                                                          # Remove indel quality, iterate backwards
    lineDict["ReadBases"] = re.sub("[+\-]\d+[ACGT]+", "", lineDict["ReadBases"])                                        # Remove indel from bases

    replacements = []
    for m in re.finditer("[,.]", lineDict["ReadBases"]):                                                                # Iterate over matched bases
        replacements.append(m.start()-i)
    for index in reversed(replacements):
        del listQuality[index]                                                                                          # Remove matched base details
    lineDict["ReadBases"] = re.sub("[,.]", "", lineDict["ReadBases"])                                                   # Remove matched base

    lineDict["ReadQuality"] = listQuality


def check_quality(lineDict):

        replacement = []
        listBases = list(lineDict["ReadBases"])
        listQuality = list(lineDict["ReadQuality"])

        for index, q in enumerate(lineDict["ReadQuality"]):                                                             # Generate list of subquality indices
            if (ord(q) - 33) < 40:
                replacement.append(index)
                lineDict["Coverage"] -= 1                                                                               # Decrement for each subquality index found
        for index in reversed(replacement):                                                                             # Remove indicies from ReadBases and ReadQuality
            del listQuality[index]
            del listBases[index]

        lineDict["ReadBases"] = listBases
        lineDict["ReadQuality"] = listQuality


def check_mismatch_percent(lineDict):                                                                                   # Checks is number of mismatches is over 80%
    if (collections.Counter(lineDict["ReadBases"]).most_common(1)[0][1]) > (lineDict["Coverage"] * .8):
        return collections.Counter(lineDict["ReadBases"]).most_common(1)[0][0]                                          # Returns mismatch base

if __name__ == '__main__':
    main()