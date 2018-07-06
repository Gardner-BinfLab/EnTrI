import pandas as pd
from Bio import SeqIO

def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

excluded = ['M5005_Spy_0009', 'M5005_Spy_0015c', 'M5005_Spy_0067', 'M5005_Spy_0072c', 'M5005_Spy_0074',
            'M5005_Spy_0075', 'M5005_Spy_0076', 'M5005_Spy_0098', 'M5005_Spy_0144c', 'M5005_Spy_0154', 'M5005_Spy_0174',
            'M5005_Spy_0181', 'M5005_Spy_0183', 'M5005_Spy_0184', 'M5005_Spy_0189', 'M5005_Spy_0190', 'M5005_Spy_0210',
            'M5005_Spy_0211', 'M5005_Spy_0234c', 'M5005_Spy_0241', 'M5005_Spy_0254c', 'M5005_Spy_0255c',
            'M5005_Spy_0259', 'M5005_Spy_0268', 'M5005_Spy_0294c', 'M5005_Spy_0333', 'M5005_Spy_0346', 'M5005_Spy_0350',
            'M5005_Spy_0396c', 'M5005_Spy_0401c', 'M5005_Spy_0405', 'M5005_Spy_0412', 'M5005_Spy_0431',
            'M5005_Spy_0448', 'M5005_Spy_0481c', 'M5005_Spy_0492c', 'M5005_Spy_0494c', 'M5005_Spy_0507c',
            'M5005_Spy_0535', 'M5005_Spy_0650', 'M5005_Spy_0666c', 'M5005_Spy_0667c', 'M5005_Spy_0669',
            'M5005_Spy_0675', 'M5005_Spy_0742', 'M5005_Spy_0744', 'M5005_Spy_0775', 'M5005_Spy_0797', 'M5005_Spy_0799',
            'M5005_Spy_0801', 'M5005_Spy_0806', 'M5005_Spy_0819', 'M5005_Spy_0853', 'M5005_Spy_0877', 'M5005_Spy_0887c',
            'M5005_Spy_0888c', 'M5005_Spy_0904', 'M5005_Spy_0912', 'M5005_Spy_0969c', 'M5005_Spy_0979c',
            'M5005_Spy_1004c', 'M5005_Spy_1041c', 'M5005_Spy_1045c', 'M5005_Spy_1049c', 'M5005_Spy_1053',
            'M5005_Spy_1074c', 'M5005_Spy_1078c', 'M5005_Spy_1087c', 'M5005_Spy_1089c', 'M5005_Spy_1166c',
            'M5005_Spy_1168c', 'M5005_Spy_1202c', 'M5005_Spy_1207c', 'M5005_Spy_1211c', 'M5005_Spy_1263c',
            'M5005_Spy_1294c', 'M5005_Spy_1312', 'M5005_Spy_1322c', 'M5005_Spy_1324', 'M5005_Spy_1369c',
            'M5005_Spy_1414c', 'M5005_Spy_1420c', 'M5005_Spy_1446c', 'M5005_Spy_1456c', 'M5005_Spy_1459c',
            'M5005_Spy_1463c', 'M5005_Spy_1510c', 'M5005_Spy_1535c', 'M5005_Spy_1541c', 'M5005_Spy_1574c',
            'M5005_Spy_1588c', 'M5005_Spy_1604c', 'M5005_Spy_1619', 'M5005_Spy_1631c', 'M5005_Spy_1643',
            'M5005_Spy_1644c', 'M5005_Spy_1645', 'M5005_Spy_1665c', 'M5005_Spy_1667c', 'M5005_Spy_1713c',
            'M5005_Spy_1717', 'M5005_Spy_1736', 'M5005_Spy_1739', 'M5005_Spy_1752c', 'M5005_Spy_1766c',
            'M5005_Spy_1792c', 'M5005_Spy_1816', 'M5005_Spy_1832c', 'M5005_Spy_1859c']

fasta = read_fasta_sequences('../embl2peptides/34-Streptococcus-pyogenes-MGAS5448.fasta')

with open('../excluded/34-Streptococcus-pyogenes-MGAS5448.fasta', 'w') as tofile:
    for name in excluded:
        if name in fasta.keys():
            tofile.write(fasta[name].format('fasta'))