import subprocess

ena_speclist = [('wickerhamomyces_anomalus','GCA_001661255.1'),
                ('candida_tropicalis','GCA_000006335.3'),
                ('yHMPu5000034957_hanseniaspora_osmophila_160519','GCA_001747045.1')
               ]

fname_base = '/mnt/g/My Drive/Crick_LMS/external_data/genomes/diverse_strains/ena/'
#             '/mnt/g/My\ Drive/Crick_LMS/external_data/genomes/diverse_strains/ena/'

for (spec, gca_name) in ena_speclist: 
    fname_in = fname_base + spec + '/' + gca_name + '.embl'
    fname_out = fname_base + spec + '/cds_raw.fasta'
    extractfeat_cmd = ['extractfeat', 
                       '-sequence', fname_in,
                       '-outseq',fname_out,
                       '-osformat2', 'fasta',
                       '-type', 'CDS',
                       '-describe', 'db_xref'
                      ]
                       
    #"'db_xref|locus_tag|protein_id'"]
    #extractfeat_cmd = ['ls', fname_base]
    #['echo', fname_in+fname_out]
    subprocess.run(extractfeat_cmd)

    

    


