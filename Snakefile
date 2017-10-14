rule all:
    input:
        'res/2017-10-10.peak.tsv'

rule start_jupyter:
    shell: 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/'

rule extract_peaks:
    input:
        rawdata = 'raw/{run}.RI.lcs_export.txt'
    output:
        'res/{run}.peak.tsv'
    shell:
        'scripts/parse_lcs_export.py {input.rawdata} > {output}'
