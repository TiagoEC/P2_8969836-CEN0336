#!/usr/bin/env python3

#importar bibliotecas necessárias
import sys

class NotaFastaError(RuntimeError):
#classe de erro para arquivos que não são fasta
    def __init__(self, message):
        self.message = message

    def __str__(self): #quando transformar o objeto em string, retornar a mensagem de erro (ex: print(erro))
        return self.message

def read_fasta(filename):
#método para pegar a informação do arquivo fasta (em dicionários)
    try:
        with open(filename, 'r') as file: #abrir o arquivo
            if not(bool(file)): #se o arquivo estiver vazio, mostrar erro
                raise RuntimeError("Arquivo vazio")
            sequences = {} #dicionário para armazenar as sequências
            current_seq = "" #string para armazenar a sequência atual
            for line in file:
                if line.startswith('>'): #se a linha for um cabeçalho, armazenar a sequência anterior e começar uma nova
                    if current_seq:
                        sequences[current_header] = current_seq 
                        current_seq = ""
                    current_header = line[1:].split()[0] #o nome da sequência é o que vem depois do > até o primeiro espaço
                else: #caso não seja um cabeçalho, adicionar a linha à sequência atual
                    current_seq += line.strip() 
            sequences[current_header] = current_seq #armazenar a última sequência

    except FileNotFoundError: #mostrar erro caso o arquivo não seja encontrado
        print(f"Erro: Arquivo (Caminho: {filename}) não encontrado.")
        exit()
    except RuntimeError as err: #mostrar erro caso o arquivo esteja vazio
        print(err)
        exit()
    else: #retornar o dicionário com as sequências caso não haja erro
        return sequences

def find_longest_orf(seq):
#médodo que encontra a maior orf na sequence, considerando os frames 1, 2, 3 e -1, -2, -3
    start_codon = 'ATG' #armazenar start códon
    stop_codons = ['TAA', 'TAG', 'TGA'] #armazenar stop códons
    # inicializar variáveis
    longest_orf = "" 
    longest_frame = 0
    longest_start = 0
    longest_end = 0

    for f in range(6): #encontrar o maior ORF dentre os 6 quadros de leitura
        if f >= 3: # reverso complementar para os frames negativos
            frame = f-3
            sequence = seq.translate(str.maketrans('ATGC', 'TACG'))[::-1][frame:] #traduzir a sequência para o reverso complementar e pegar o frame
        else: # frames positivos
            frame = f
            sequence = seq[frame:]
        
        start = -1 #inicializar variável para armazenar o início do ORF, indicando que não há ORF no momento
        for i in range(0, len(sequence), 3): #percorrer a sequência de 3 em 3
            codon = sequence[i: i+3] #pegar o códon
            if codon == start_codon and start == -1: #se for um start códon e não estiver dentro de um ORF, começar a sequência
                start = i
            elif codon in stop_codons and start != -1: #se estiver dentro de um ORF e for códon de parada, terminar a sequência
                end = i + 3 #o final da sequência é o número do começo do códon mais 3 (o tamanho do códon, considerando que i começa em 0)
                orf = sequence[start: end] #pegar a sequência do ORF
                if len(orf) > len(longest_orf): #caso essa sequência seja maior que a anterior, atualizar os resultados
                    longest_orf = orf
                    longest_frame = f+1 #armazenar o frame, somando 1 para que o frame 0 seja o frame 1
                    longest_start = start + frame + 1 #o começo dessa ORF na sequência é o começo da sequência mais o frame mais 1
                    longest_end = end + frame #o final dessa ORF na sequência é o final da sequência mais o frame
                start = -1 #recomeçar a sequência
    return longest_orf, longest_frame, longest_start, longest_end

def translate_to_protein(seq):
#método que traduz a sequência de DNA para sequência de proteína
    translation_table = { #dicionário com as traduções
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
        'AAT':'N', 'AAC':'N',
        'GAT':'D', 'GAC':'D',
        'TGT':'C', 'TGC':'C',
        'CAA':'Q', 'CAG':'Q',
        'GAA':'E', 'GAG':'E',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
        'CAT':'H', 'CAC':'H',
        'ATT':'I', 'ATC':'I', 'ATA':'I',
        'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
        'AAA':'K', 'AAG':'K',
        'ATG':'M',
        'TTT':'F', 'TTC':'F',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'TGG':'W',
        'TAT':'Y', 'TAC':'Y',
        'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
        'TAA':'*', 'TGA':'*', 'TAG':'*'
    }
    protein = ""
    for i in range(0, len(seq), 3): #percorrer a sequência de 3 em 3
        codon = seq[i: i+3] #pegar o códon
        if codon in translation_table: #se o códon estiver no dicionário de tradução, adicionar o aminoácido à sequência de proteína
            protein += translation_table[codon]
        elif len(codon) != 3: #se o códon não tiver 3 nucleotídeos, não fazer nada (não deve acontecer, apenas por sanidade)
            pass
        else:
            protein += "?" #Não deve acontecer, mas para sanidade
    return protein

def main(): 
    #função principal
    #Complexidade: O(MN) onde M é o número de sequências e N é o tamanho médio das sequências
    try:    
        file = sys.argv[1]
        if not file.endswith(".fasta") and not file.endswith(".fa"):
            # se o arquivo não for um arquivo fasta, mostrar erro
            raise NotaFastaError("Erro: Arquivo não é um arquivo fasta.")
    except IndexError:
        # se o arquivo não for fornecido, mostrar erro
        print(f'Erro: É necessário um arquivo de entrada')
        exit()
    except NotaFastaError as err:
        print(err)
        exit()
    
    sequences = read_fasta(file) #pegar as sequências do arquivo pelo método read_fasta
    try:
        with open("ORF.fna", 'w') as fna, open("ORF.faa", 'w') as faa: #criar os arquivos de output
            for header, seq in sequences.items():
                orf, frame, start, end = find_longest_orf(seq) #pegar os resultados da função find_longest_orf
                peptide = translate_to_protein(orf) #traduzir a sequência do ORF para proteína pelo método translate_to_protein
                #salvar as informações nos arquivos
                fna.write(f'>{header}_frame{frame}_{start}_{end}\n{orf}\n') 
                faa.write(f'>{header}_frame{frame}_{start}_{end}\n{peptide}\n')
    except:
        # caso não seja possível criar os arquivos, mostrar erro
        print("Erro ao criar os arquivos de output.")

if __name__ == "__main__":
    #se o script for executado diretamente, rodar a função main
    main()