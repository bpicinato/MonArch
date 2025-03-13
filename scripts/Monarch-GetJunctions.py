#! /usr/bin/python3
import sys

# v_9-6-21
# adições/correções:
# 1) output corrigido de acordo com overlap ou gap
#    - escolhi corrigir sempre a ponta que for o final do círculo
# 2) considera o delta como uma condição para achar second_hit
# 3) correção das condições com max_circle_size, considerando que a janela em
#    que a segunda metade vai ser procurada corresponde ao tamanho máximo do
#    círculo que queremos

# função que armazena o best_hit e procura pelo second_hit
def GetCircles(BlastDict, delta, max_circle_size, min_sec_hit_size, min_read_cov, blast_results, invert_strand, strandness):
    with open(blast_results+".bed", "w") as outbed:
        with open(blast_results+".info","w") as outinfo:
            # variavel para contar quantos círculos são achados com cada parâmetro
            count = 0
            # procurando pelas duas metades em cada uma das reads do dicionário
            for read, alignments in BlastDict.items():

                # o melhor alinhamento será o primeiro (provavelmente o maior e com o
                # menor evalue)
                best_hit = alignments[0]

                # extraindo variáveis do blast do melhor alinhamento
                chromosome = best_hit['chr']
                qstart_best_hit = best_hit['qstart']
                qend_best_hit = best_hit['qend']
                sstart_best_hit = best_hit['sstart']
                send_best_hit = best_hit['send']
                best_hit_strand = best_hit['sstrand']
                read_size = int(read.split('__')[1])
                best_hit_size = qend_best_hit - qstart_best_hit + 1

                # antes de tudo vamos pegar o best_hit e verificar se ele tem o tamanho
                # mínimo necessário read_size/2
                #só vai fazer tudo embaixo se essa primeira condição for satisfeita
                # alguma garantia de evitar fazer comparações desnecessárias se o melhor
                # hit já é ruim (caso em que a read não tem junções e alinha vários
                # pedaços pequenos no genoma)
                if best_hit_size >= read_size/2:

                    second_hit = " "
                    best_hit_half = "A or B"

                    # verificando qual lado é o best_hit
                    # notar que troquei para verificar o final do best_hit
                    if qend_best_hit in range(int(read_size/2), read_size - min_sec_hit_size):
                        best_hit_half = "A"
                    else:
                        best_hit_half = "B"

                    # buscando second_hit
                    # 1) fazendo lista de todos os starts e ends da read
                    #    caso tenha mesma fita que o best_hit. caso não
                    #    seja, será atribuido valor zero, apenas para que
                    #    o índice do sstart não se altere
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

                    # 2) vendo cada caso de 'fita' e 'metade':
                    if best_hit_strand == 'plus' and best_hit_half == "A":

                        # range permitido para a busca do second_hit (tamanho do circ)
                        from_chr = abs(send_best_hit - max_circle_size)
                        to_chr = sstart_best_hit - 1
                        sec_hit_idx = [all_sstart.index(x) for x in all_sstart if x in range(from_chr, to_chr + 1)]

                        # se achou um index, existe o second_hit (lista não vazia)
                        #if sec_hit_idx:
                        # testing not allowing multimappers
                        if len(sec_hit_idx) == 1:
                            # pegando variáveis importantes do second_hit
                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_second_hit - qend_best_hit - 1
                            # condiçoes que o second_hit deve obedecer
                            # cobertura mínima
                            if (qend_second_hit - qstart_best_hit + 1) >= (min_read_cov*read_size):
                                # overlap ou gap tem que ser menor ou igual a delta
                                if abs(overlap) <= delta:
                                    # tamanho mínimo
                                    if (second_hit_size > min_sec_hit_size):

                                        # corrigindo coordenadas do círculo com o delta
                                        # send_best_hit = send_best_hit + overlap
                                        # MUDANDO PARA ALTERAR APENAS QUANDO OVERLAP <= 0
                                        if overlap <= 0:
                                            send_best_hit = send_best_hit + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome, str(second_hit['sstart']),
                                                                    str(send_best_hit), read.split("__")[0],
                                                                    str(score), strand])+"\n")


                                        # sys.stderr.writeando output como antes
                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'A', len(sec_hit_idx), overlap,
                                        qstart_best_hit, qend_best_hit, qstart_second_hit, qend_second_hit,
                                        sstart_best_hit, send_best_hit, sstart_second_hit, send_second_hit,'\n']]))


                                        count+=1

                    elif best_hit_strand == 'plus' and best_hit_half == "B":

                        # range permitido para a busca do second_hit (tamanho do circ)
                        from_chr = send_best_hit + 1
                        to_chr = max_circle_size + sstart_best_hit
                        sec_hit_idx = [all_send.index(x) for x in all_send if x in range(from_chr, to_chr + 1)]

                        # se achou um index, existe o second_hit (lista não vazia)
                        #if sec_hit_idx:
                        # testing not allowing multimappers
                        if len(sec_hit_idx) == 1:
                            # pegando variáveis importantes do second_hit
                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_best_hit - qend_second_hit - 1
                            # condiçoes que o second_hit deve obedecer
                            # cobertura mínima
                            if (qend_best_hit - qstart_second_hit + 1) >= (min_read_cov*read_size):
                                # overlap ou gap tem que ser menor ou igual a delta
                                if abs(overlap) <= delta:
                                    # tamanho mínimo
                                    if second_hit_size > min_sec_hit_size:

                                        # corrigindo coordenadas do círculo com o delta
                                        # send_best_hit = send_best_hit + overlap
                                        # MUDANDO PARA ALTERAR APENAS QUANDO OVERLAP <= 0
                                        if overlap <= 0:
                                            send_second_hit = second_hit['send'] + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome,  str(sstart_best_hit),
                                                                    str(send_second_hit), read.split("__")[0],
                                                                    str(score), strand])+"\n")


                                        # sys.stderr.writeando output como antes
                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'B', len(sec_hit_idx), overlap,
                                        qstart_second_hit, qend_second_hit, qstart_best_hit, qend_best_hit,
                                        sstart_second_hit, send_second_hit, sstart_best_hit, send_best_hit,'\n']]))

                                        count+=1

                    elif best_hit_strand == 'minus' and best_hit_half == "A":

                        # range permitido para a busca do second_hit (tamanho do circ)
                        from_chr = sstart_best_hit + 1
                        to_chr = send_best_hit + max_circle_size
                        sec_hit_idx = [all_sstart.index(x) for x in all_sstart if x in range(from_chr, to_chr + 1)]

                        # se achou um index, existe o second_hit (lista nao vazia)
                        #if sec_hit_idx:
                        # testing not allowing multimappers
                        if len(sec_hit_idx) == 1:
                            # pegando variaveis importantes do second_hit
                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_second_hit - qend_best_hit - 1
                            # condiçoes que o second_hit deve obedecer
                            # cobertura mínima
                            if (qend_second_hit - qstart_best_hit + 1) >= (min_read_cov*read_size):
                                # overlap ou gap tem que ser menor ou igual a delta
                                if abs(overlap) <= delta:
                                    # tamanho minimo
                                    if (second_hit_size > min_sec_hit_size):

                                        # corrigindo coordenadas do círculo com o delta
                                        # send_best_hit = send_best_hit + overlap
                                        # MUDANDO PARA ALTERAR APENAS QUANDO OVERLAP <= 0
                                        if overlap <= 0:
                                            sstart_second_hit = second_hit['sstart'] + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome,
                                                                    str(send_best_hit), str(sstart_second_hit),
                                                                    read.split("__")[0], str(score), strand])+"\n")


                                        # sys.stderr.writeando output como antes
                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'A', len(sec_hit_idx), overlap,
                                        qstart_best_hit, qend_best_hit, qstart_second_hit, qend_second_hit,
                                        sstart_best_hit, send_best_hit, sstart_second_hit, send_second_hit, '\n']]))

                                        count+=1

                    elif best_hit_strand == 'minus' and best_hit_half == "B":

                        # range permitido para a busca do second_hit (tamanho do circ)
                        from_chr = abs(max_circle_size - sstart_best_hit)
                        to_chr = send_best_hit - 1
                        sec_hit_idx = [all_send.index(x) for x in all_send if x in range(from_chr, to_chr + 1)]

                        # se achou um index, existe o second_hit (lista nao vazia)
                        #if sec_hit_idx:
                        # testing not allowing multimappers
                        if len(sec_hit_idx) == 1:
                            # pegando variaveis importantes do second_hit
                            second_hit = alignments[sec_hit_idx[0]]
                            qend_second_hit = second_hit['qend']
                            qstart_second_hit = second_hit['qstart']
                            send_second_hit = second_hit['send']
                            sstart_second_hit = second_hit['sstart']
                            second_hit_size = qend_second_hit - qstart_second_hit + 1
                            overlap = qstart_best_hit - qend_second_hit - 1
                            # condiçoes que o second_hit deve obedecer
                            # cobertura mínima
                            if (qend_best_hit - qstart_second_hit + 1) >= (min_read_cov*read_size):
                                # overlap ou gap tem que ser menor ou igual a delta
                                if abs(overlap) <= delta:
                                    # tamanho minimo
                                    if (second_hit_size > min_sec_hit_size):

                                        # corrigindo coordenadas do círculo com o delta
                                        # send_best_hit = send_best_hit + overlap
                                        # MUDANDO PARA ALTERAR APENAS QUANDO OVERLAP <= 0
                                        if overlap <= 0:
                                            sstart_best_hit = sstart_best_hit + overlap
                                        if best_hit_strand == second_hit['sstrand']:
                                            score = 100 + overlap
                                        else:
                                            score = 50 + overlap
                                        outbed.write("\t".join([chromosome,
                                                                    str(send_second_hit), str(sstart_best_hit),
                                                                    read.split("__")[0], str(score), strand])+"\n")


                                            # sys.stderr.writeando output como antes
                                        outinfo.write('\t'.join([str(x) for x in [chromosome, read.split('__')[0],
                                        strand, 'B', len(sec_hit_idx), overlap,
                                        qstart_second_hit, qend_second_hit, qstart_best_hit, qend_best_hit,
                                        sstart_second_hit, send_second_hit, sstart_best_hit, send_best_hit,'\n']]))

                                        count+=1

def ParseBlastResults(file, max_target_seqs, max_mismatch):
    BlastDict = dict()
    # variáveis para ajudar a colocar só os max_target_seqs top alinhamentos em Dict
    read = ''
    count = max_target_seqs
    with open(file, 'r') as blast:
        for row in blast.readlines():
            row = row.strip().split('\t')
            # lendo colunas do blast
            sseqid, qseqid, qlen, length, qstart, qend, sstart, send, sstrand, mismatch = row
            # chaves do Dict serão 'nome-da-read__read-size' (2 underline)
            DictKey = qseqid+'__'+str(qlen)
            # esse alinhamento será tbm um dicionário
            # já defini aqui o que é int e o que é string
            ThisAlignment = {'qstart': int(qstart),
                             'qend': int(qend),
                             'sstart': int(sstart),
                             'send': int(send),
                             'sstrand': sstrand,
                             'chr': sseqid,
                             'mismatch': int(mismatch)}
            # reseta count e read para ir para a próxima read
            if read != DictKey:
                read = DictKey
                count = max_target_seqs
            # colocando no Dict apenas os max_target_seqs top alinhamentos
            # com no máximo max_mismatch alinhamentos
            if count > 0:
                # construindo Dict
                # verifica se é o primeiro alinhamento dessa read. se não for,
                # deve usar o append()
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
