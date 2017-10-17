# Config {{{{{{1
configfile: 'pipeline.config.yaml'
max_threads = 99
# END Config }}}}}
# Library {{{{{1
def mothur(command, **kwargs):
    return ( "mothur -q '#{command}(".format(command=command)
           + ", ".join(["{key}={val}".format(key=key, val=kwargs[key]) for key in kwargs])
           + ")'")

def unlines(*lines):
    return '\n'.join(lines)

# END Library }}}}}
# Shortcuts {{{{{1

rule all:
    input:
        'res/2017-10-10.peak.tsv'

rule start_jupyter:
    shell: 'jupyter notebook --config=nb/jupyter_notebook_config.py --notebook-dir=nb/'

# END Shortcuts }}}}}
# HPLC {{{{{1

rule extract_peak_data:
    input:
        script = 'scripts/parse_lcs_export.py',
        rawdata = 'raw/hplc/{run}.RI.lcs_export.txt'
    output:
        'res/{run}.peak.tsv'
    shell:
        '{input.script} {input.rawdata} > {output}'

# END HPLC }}}}}

# RRS {{{{{1

rule concat_fused_fn:
    input:
        expand('seq/split/{library_id}.rrs.fuse.fn',
               library_id=config['libraries'])
    output: 'seq/rrs.fuse.fn'
    shell:
        'cat {input} > {output}'

rule make_groups_file:
    input:
           splits = expand('seq/split/{library_id}.rrs.fuse.fn',
                           library_id=config['libraries']),
           script = 'scripts/get_groups.sh'
    output: 'res/rrs.fuse.groups'
    shell:
        '{input.script} {input.splits} > {output}'

rule link_rrs_data:
    input:
        lambda wildcards: \
            ancient(expand('raw/rrs/{library_slug}_L001_R{read_number}_001.fastq.gz',
                           library_slug=config["libraries"][wildcards.library_id],
                           read_number=wildcards.read_number))
    output:
        'seq/split/{library_id}.rrs.r{read_number}.fastq.gz'
    shell:
        'ln -rs {input} {output}'

rule fuse_paired_reads:
    input:
        r1 = 'seq/split/{library_id}.rrs.r1.fastq.gz',
        r2 = 'seq/split/{library_id}.rrs.r2.fastq.gz'
    output:
        'seq/split/{library_id}.rrs.fuse.fn'
    shadow: 'full'
    shell:
        unlines('gzip -dc {input.r1} > {wildcards.library_id}.r1.fastq',
                'gzip -dc {input.r2} > {wildcards.library_id}.r2.fastq',
                mothur('make.contigs',
                       ffastq='{wildcards.library_id}.r1.fastq',
                       rfastq='{wildcards.library_id}.r2.fastq',
                       trimoverlap='T'),
                'mv {wildcards.library_id}.r1.trim.contigs.fasta seq/split/{wildcards.library_id}.rrs.fuse.fn'
               )

rule screen_seqs:
    input:
        seqs = 'seq/{prefix}.fn',
        groups = 'res/{prefix}.groups'
    output:
        seqs = 'seq/{prefix}.screen1.fn',
        groups = 'res/{prefix}.screen1.groups'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('screen.seqs',
                       fasta='{input.seqs}',
                       group='{input.groups}',
                       maxambig=str(0),
                       minlength=str(180),
                       maxlength=str(295),
                       processors='{threads}'),
                'mv seq/{wildcards.prefix}.good.fn {output.seqs}',
                'mv seq/{wildcards.prefix}.good.groups {output.groups}'
               )

rule uniq_seqs_fn:
    input:
        seqs = 'seq/{prefix}.fn',
    output:
        seqs = 'seq/{prefix}.uniq.fn',
        names = 'res/{prefix}.uniq.names'
    shadow: 'full'
    shell:
        unlines(mothur('unique.seqs', fasta='{input.seqs}'),
                'mv seq/{wildcards.prefix}.unique.fn {output.seqs}',
                'mv seq/{wildcards.prefix}.names {output.names}',
               )

rule uniq_seqs_afn:
    input:
        seqs = 'seq/{prefix}.afn',
    output:
        seqs = 'seq/{prefix}.uniq.afn',
        names = 'res/{prefix}.uniq.names'
    shadow: 'full'
    shell:
        unlines(mothur('unique.seqs', fasta='{input.seqs}'),
                'mv seq/{wildcards.prefix}.unique.afn {output.seqs}',
                'mv seq/{wildcards.prefix}.names {output.names}',
               )

rule count_groups_from_scratch:
    input:
        names = 'res/{prefix}.fuse.screen1.uniq.names',
        groups = 'res/{prefix}.fuse.screen1.groups'
    output:
        count_table = 'res/{prefix}.fuse.screen1.uniq.count_table'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('count.seqs',
                       name='{input.names}',
                       group='{input.groups}',
                       processors='{threads}'),
               )

rule uniq_seqs_groups_afn:
    input:
        seqs = 'seq/{prefix}.afn',
    output:
        seqs = 'seq/{prefix}.uniq.afn',
        names = 'res/{prefix}.uniq.names'
    shadow: 'full'
    shell:
        unlines(mothur('uniq.seqs', fasta='{output.seqs}'),
                'mv seq/{wildcards.prefix}.unique.afn {output.seqs}',
                'mv seq/{wildcards.prefix}.names {output.names}',
               )

rule align_rrs:
    input:
        seqs = 'seq/{prefix}.fn',
        ref = 'ref/silva.nr.pcr_v4.afn'
    output:
        align = 'seq/{prefix}.align.afn'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('align.seqs', fasta='{input.seqs}', reference='{input.ref}', processors='{threads}'),
                'mv seq/{wildcards.prefix}.align seq/{wildcards.prefix}.align.afn')


rule count_groups_adhoc:
    input:
        'res/{prefix}.groups'
    output:
        'res/{prefix}.groups.tally.tsv'
    shell:
        "cut -f2 {input}"
        " | uniq -c"
        " | awk -v OFS='\t' '{{print $2, $1}}'"
        " > {output}"

rule blast_spike_seq:
    input:
        spike_seqs = 'meta/spike.fn',
        reads = 'seq/rrs.fuse.fn'
    output:
        'res/rrs.fuse.blastn-spike.tsv'
    params:
        ident_thresh = 95
    shell:
        'blastn -query {input.reads}'
        ' -subject {input.spike_seqs}'
        ' -max_target_seqs 1'
        ' -perc_identity {params.ident_thresh}'
        ' -outfmt 6 -out {output}'
