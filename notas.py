#!/usr/bin/env python3

# Variáveis
total = 0
contador_notas = 0

while contador_notas < 10: #Pegar exatamente 10 notas do usuário
    nota = float(input("Digite a nota: ")) #recebe a nota do usuário pelo terminal e converte para float
    total += nota #soma a nota ao total
    contador_notas += 1 #para não entrar em loop infinito
media = total / 10 #calcula a média

print(media)

# Windows:
# certutil -hashfile notas.py MD5
# Linux:
# md5sum notas.py