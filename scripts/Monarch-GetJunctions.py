#! /usr/bin/python3
import sys

# function to store the best_hit and search for the second_hit
def GetCircles(BlastDict, delta, max_circle_size, min_sec_hit_size, min_read_cov, blast_results, invert_strand, strandness):
    with open(blast_results+".bed", "w") as outbed:
        with open(blast_results+".info","w") as outinfo:
            # count how many circles are found with each parameter
            count = 0
            # searching for each half in each read in the dictionary
            for read, alignments in BlastDict.items():

                # the best alignment is the first one (and probably with the lowest evalue)
                best_hit = alignments[0]

                # extracting variables of the first alignment
                chromosome = best_hit['chr']
                qstart_best_hit = best_hit['qstart']
                qend_best_hit = best_hit['qend']
                sstart_best_hit = best_hit['sstart']
                send_best_hit = best_hit['send']
                best_hit_strand = best_hit['sstrand']
                read_size = int(read.split('__')[1])
                best_hit_size = qend_best_hit - qstart_best_hit + 1

                # check if the best hit has the min size read_size/2
                # continue only if this is met
                # if the best hit is already bad, do not even continue
                if best_hit_size >= read_size/2:

                    second_hit = " "
                    best_hit_half = "A or B"

                    # checking which half is the best_hit
                    if qend_best_hit in range(int(read_size/2), read_size - min_sec_hit_size):
                        best_hit_half = "A"
                    else:
                        best_hit_half = "B"

                    # searching for second_hit
                    # 1) getting all hits in the same strand as the best hit
                    #    if not, it gets 0 just sstart index does not change
                    all_sstart = [int(ali['sstart']) if ali['sstrand'] == best_hit_strand else 0 for ali in alignments ]
                    all_send = [int(ali['send']) if ali['sstrand'] == best_hit_strand else 0 for ali in alignments ]

                    # inverting strand to be written in .bed files
                    if invert_strand == 'YES' and strandness == 'YES':
                        if best_hit_strand == 'plus':
                            strand = "-"
                        else:
                            strand = "+"
                    # or keeping it the same
                    elif strandness == 'YES':
                        if best_hit_strand == 'plus':
                            strand = "+"
                        else:
                            strand = "-"

                    else:
                        strand = '+'

                    # 2) checking each case of "strand" and "half"
                    if best_hit_strand == 'plus' and best_hit_half == "A":

                        # allowed range for second_hit search (max_circ_size parameter)
                        from_chr = abs(send_best_hit - max_circle_size)
                        to_chr = sstart_best_hit - 1
                        sec_hit_idx = [all_sstart.index(x) for x in all_sstart if x in range(from_chr, to_chr + 1)]

                        # if there is an index, there is a second_hit
                        # do not allow multimapper reads - only one hit
                        if len(sec_hit_idx) == 1:
                            # getting variables of second_hit
                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_second_hit - qend_best_hit - 1
                            # constrains on second_hit
                            # minimum coverage
                            if (qend_second_hit - qstart_best_hit + 1) >= (min_read_cov*read_size):
                                # overlap or gap maximum size (delta parameter)
                                if abs(overlap) <= delta:
                                    # minimum size
                                    if (second_hit_size > min_sec_hit_size):

                                        # correcting circ coordinates based on delta
                                        # send_best_hit = send_best_hit + overlap
                                        # only when overlap <= 0!
                                        if overlap <= 0:
                                            send_best_hit = send_best_hit + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome, str(second_hit['sstart']),
                                                                    str(send_best_hit), read.split("__")[0],
                                                                    str(score), strand])+"\n")


                                        # writing output
                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'A', len(sec_hit_idx), overlap,
                                        qstart_best_hit, qend_best_hit, qstart_second_hit, qend_second_hit,
                                        sstart_best_hit, send_best_hit, sstart_second_hit, send_second_hit,'\n']]))


                                        count+=1

                    elif best_hit_strand == 'plus' and best_hit_half == "B":

                        # proceed as first block...
                        from_chr = send_best_hit + 1
                        to_chr = max_circle_size + sstart_best_hit
                        sec_hit_idx = [all_send.index(x) for x in all_send if x in range(from_chr, to_chr + 1)]

                        if len(sec_hit_idx) == 1:
                            
                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_best_hit - qend_second_hit - 1

                            if (qend_best_hit - qstart_second_hit + 1) >= (min_read_cov*read_size):

                                if abs(overlap) <= delta:

                                    if second_hit_size > min_sec_hit_size:
                                        
                                        if overlap <= 0:
                                            send_second_hit = second_hit['send'] + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome,  str(sstart_best_hit),
                                                                    str(send_second_hit), read.split("__")[0],
                                                                    str(score), strand])+"\n")

                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'B', len(sec_hit_idx), overlap,
                                        qstart_second_hit, qend_second_hit, qstart_best_hit, qend_best_hit,
                                        sstart_second_hit, send_second_hit, sstart_best_hit, send_best_hit,'\n']]))

                                        count+=1

                    elif best_hit_strand == 'minus' and best_hit_half == "A":

                        # proceed as first block...
                        from_chr = sstart_best_hit + 1
                        to_chr = send_best_hit + max_circle_size
                        sec_hit_idx = [all_sstart.index(x) for x in all_sstart if x in range(from_chr, to_chr + 1)]

                        if len(sec_hit_idx) == 1:

                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_second_hit - qend_best_hit - 1

                            if (qend_second_hit - qstart_best_hit + 1) >= (min_read_cov*read_size):

                                if abs(overlap) <= delta:

                                    if (second_hit_size > min_sec_hit_size):
                                        
                                        if overlap <= 0:
                                            sstart_second_hit = second_hit['sstart'] + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome,
                                                                    str(send_best_hit), str(sstart_second_hit),
                                                                    read.split("__")[0], str(score), strand])+"\n")

                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'A', len(sec_hit_idx), overlap,
                                        qstart_best_hit, qend_best_hit, qstart_second_hit, qend_second_hit,
                                        sstart_best_hit, send_best_hit, sstart_second_hit, send_second_hit, '\n']]))

                                        count+=1

                    elif best_hit_strand == 'minus' and best_hit_half == "B":

                        # proceed as first block...
                        from_chr = abs(max_circle_size - sstart_best_hit)
                        to_chr = send_best_hit - 1
                        sec_hit_idx = [all_send.index(x) for x in all_send if x in range(from_chr, to_chr + 1)]

                        if len(sec_hit_idx) == 1:

                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_best_hit - qend_second_hit - 1

                            if (qend_best_hit - qstart_second_hit + 1) >= (min_read_cov*read_size):

                                if abs(overlap) <= delta:

                                    if (second_hit_size > min_sec_hit_size):
                                        
                                        if overlap <= 0:
                                            sstart_best_hit = sstart_best_hit + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome,
                                                                    str(send_second_hit), str(sstart_best_hit),
                                                                    read.split("__")[0], str(score), strand])+"\n")

                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'B', len(sec_hit_idx), overlap,
                                        qstart_second_hit, qend_second_hit, qstart_best_hit, qend_best_hit,
                                        sstart_second_hit, send_second_hit, sstart_best_hit, send_best_hit,'\n']]))

                                        count+=1

def ParseBlastResults(file, max_target_seqs, max_mismatch):
    BlastDict = dict()
    # variables to help put only max_target_seqs top alignments in the dict
    read = ''
    count = max_target_seqs
    with open(file, 'r') as blast:
        for row in blast.readlines():
            row = row.strip().split('\t')
            # reading blast columns
            sseqid, qseqid, qlen, length, qstart, qend, sstart, send, sstrand, mismatch = row
            # dict keys are 'read-name__read-size' (2 underline)
            DictKey = qseqid+'__'+str(qlen)
            # alignment is also a dict
            ThisAlignment = {'qstart': int(qstart),
                             'qend': int(qend),
                             'sstart': int(sstart),
                             'send': int(send),
                             'sstrand': sstrand,
                             'chr': sseqid,
                             'mismatch': int(mismatch)}
            # reset count and read to go to the next read
            if read != DictKey:
                read = DictKey
                count = max_target_seqs
            # putting in the dict only top alignments
            # with maximum max_mismatch 
            if count > 0:
                # building dict
                # checks if it is the first alignment of the read. if not, 
                # use append()
                if DictKey not in BlastDict.keys():
                    # only add to dictionary if number of mismatches in
                    # ThisAlignment is lower than max_mismatch
                    if ThisAlignment['mismatch'] <= max_mismatch:
                        BlastDict[DictKey] = [ThisAlignment]
                        count -= 1
                else:
                    # only add to dictionary if number of mismatches in
                    # ThisAlignment is lower than max_mismatch
                    if ThisAlignment['mismatch'] <= max_mismatch:
                        BlastDict[DictKey].append(ThisAlignment)
                        count -= 1
    return BlastDict

def main():
    blast_results = sys.argv[1]
    delta = int(sys.argv[2])
    max_circle_size = int(sys.argv[3])
    min_sec_hit_size = int(sys.argv[4])
    min_read_cov = float(sys.argv[5])
    max_mismatch = int(sys.argv[6])
    invert_strand = sys.argv[7]
    strandness = sys.argv[8]
    # Parsing results
    BlastDict = ParseBlastResults(blast_results, 250, max_mismatch)
    # Get Circles
    GetCircles(BlastDict, delta, max_circle_size,
               min_sec_hit_size, min_read_cov, blast_results, invert_strand, strandness)

if __name__ == '__main__':
        main()
