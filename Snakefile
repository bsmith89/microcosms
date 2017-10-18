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

just_symlink = 'ln -rs {input} {output}'

# END Library }}}}}
# Shortcuts {{{{{1

rule all:
    input:
        'res/2017-10-10.peak.tsv',
        'res/rrs.fuse.screen1.uniq.align.screen2.press.screen3.screen4.clust.reps.blastn-sra.tsv'

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
        'seq/{prefix}.align.afn'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('align.seqs',
                       fasta='{input.seqs}',
                       reference='{input.ref}',
                       processors='{threads}'),
                'mv seq/{wildcards.prefix}.align {output}')

rule screen_aligned_seqs:
    input:
        seqs = 'seq/{prefix}.align.afn',
        counts = 'res/{prefix}.count_table'
    output:
        seqs = 'seq/{prefix}.align.screen2.afn',
        counts = 'res/{prefix}.align.screen2.count_table'
    shadow: 'full'
    shell:
        unlines(mothur('screen.seqs',
                       fasta='{input.seqs}',
                       start=str(3100),
                       count='{input.counts}',
                       end=str(10600),
                       maxhomop=str(8)),
                'mv seq/{wildcards.prefix}.align.good.afn {output.seqs}',
                'mv seq/{wildcards.prefix}.good.count_table {output.counts}'
               )

rule press_alignment:
    input: 'seq/{prefix}.afn'
    output: 'seq/{prefix}.press.afn'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('filter.seqs',
                       fasta='{input}',
                       vertical='T',
                       trump='.',
                       processors='{threads}'),
                'mv seq/{wildcards.prefix}.filter.fasta {output}')

# Changing the alignment of sequences (pressing or aligning) doesn't change the
# counts, so we just symlink the file to keep naming conventions.
rule symlink_identical_counts:
    input: 'res/{prefix}.count_table'
    wildcard_constraints:
        last_action = '(press|align)'
    output: 'res/{prefix}.{last_action}.count_table'
    shell:
        just_symlink

rule screen_chimeras:
    input:
        seqs = 'seq/{prefix}.afn',
        counts = 'res/{prefix}.count_table'
    output:
        seqs = 'seq/{prefix}.screen3.afn',
        counts = 'res/{prefix}.screen3.count_table'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('chimera.uchime',
                      fasta='{input.seqs}',
                      count='{input.counts}',
                      dereplicate='T',
                      processors='{threads}'),
                mothur('remove.seqs',
                       fasta='{input.seqs}',
                       accnos='seq/{wildcards.prefix}.denovo.uchime.accnos'),
                'mv seq/{wildcards.prefix}.pick.afn {output.seqs}',
                'mv seq/{wildcards.prefix}.denovo.uchime.pick.count_table {output.counts}'
                )

# fuse.screen1.uniq.align.screen2.press.uniq.screen3.screen4

rule classify_seqs:
    input:
        seqs = 'seq/{prefix}.afn',
        counts = 'res/{prefix}.count_table',
        ref_seqs = 'ref/silva.nr.pcr_v4.fn',
        ref_tax = 'ref/silva.nr.tax',
    output:
        'res/{prefix}.tax'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(mothur('classify.seqs',
                       fasta='{input.seqs}',
                       count='{input.counts}',
                       reference='{input.ref_seqs}',
                       taxonomy='{input.ref_tax}',
                       method='wang',
                       cutoff=str(0),
                       processors='{threads}'),
                'mv seq/{wildcards.prefix}.nr.wang.taxonomy {output}'
               )

rule screen_taxa:
    input:
        seqs = 'seq/{prefix}.afn',
        counts = 'res/{prefix}.count_table',
        tax = 'res/{prefix}.tax'
    output:
        seqs = 'seq/{prefix}.screen4.afn',
        counts = 'res/{prefix}.screen4.count_table',
        tax = 'res/{prefix}.screen4.tax'
    shadow: 'full'
    shell:
        unlines(mothur('remove.lineage',
                       fasta='{input.seqs}',
                       count='{input.counts}',
                       taxonomy='{input.tax}',
                       taxon='Chloroplast-Mitochondria-Archaea-Eukaryota-Unknown-Sphingopyxis'),
                'mv seq/{wildcards.prefix}.pick.afn {output.seqs}',
                'mv res/{wildcards.prefix}.pick.count_table {output.counts}',
                'mv res/{wildcards.prefix}.pick.tax {output.tax}'
               )

# Need to break ambiguity for producing *.screen4.tax, and screen_taxa
# doesn't require re-classifying anything.
ruleorder: screen_taxa > classify_seqs

rule cluster_otus:
    input:
        seqs = 'seq/{prefix}.afn',
        counts = 'res/{prefix}.count_table',
        tax = 'res/{prefix}.tax'
    output: 'res/{prefix}.clust.otus'
    shadow: 'full'
    threads: max_threads
    shell:
        unlines(
        mothur('cluster.split',
               fasta='{input.seqs}',
               count='{input.counts}',
               taxonomy='{input.tax}',
               method='opti',
               splitmethod='classify',
               taxlevel=str(2),
               cutoff='0.03',
               processors='{threads}'),
        'mv seq/{wildcards.prefix}.opti_mcc.unique_list.list {output}',
               )

rule make_shared:
    input:
        otus = 'res/{prefix}.clust.otus',
        counts = 'res/{prefix}.count_table'
    output:
        'res/{prefix}.clust.shared'
    shell:
        unlines(
        mothur('make.shared',
        list='{input.otus}',
        count='{input.counts}',
        label='0.03'),
        )

rule classify_otus:
    input:
        otus = 'res/{prefix}.clust.otus',
        counts = 'res/{prefix}.count_table',
        tax = 'res/{prefix}.tax'
    output:
        'res/{prefix}.clust.tax'
    shell:
        unlines(
        mothur('classify.otu',
        list='{input.otus}',
        count='{input.counts}',
        taxonomy='{input.tax}',
        label='0.03'),
        'mv res/{wildcards.prefix}.clust.0.03.cons.taxonomy {output}'
        )

rule get_otu_reps:
    input:
        otus = 'res/{prefix}.clust.otus',
        counts = 'res/{prefix}.count_table',
        seqs = 'seq/{prefix}.afn'
    output:
        'seq/{prefix}.clust.reps.afn'
    shell:
        unlines(
        mothur('get.oturep',
        method='abundance',
        list='{input.otus}',
        count='{input.counts}',
        fasta='{input.seqs}'),
        r"sed 's:>[^ ]\+\s\+\(Otu[0-9]\+\)|.*$:>\1:' < res/{wildcards.prefix}.clust.0.03.rep.fasta > {output}"
        )

# END RRS }}}}}}
# Misc. {{{{{1
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

rule best_reference_hit_afn:
    input:
        subject = 'ref/{reference}.fn',
        query = 'seq/{prefix}.fn'
    output:
        'res/{prefix}.blastn-{reference}.tsv'
    shell:
        'blastn -query {input.query}'
        ' -subject {input.subject}'
        ' -max_target_seqs 1'
        ' -outfmt 6 -out {output}'

rule unalign_nucl_seqs:
    input:
        '{prefix}.afn'
    output:
        '{prefix}.fn'
    run:
        from Bio.SeqIO import parse, write

        def unalign(rec):
            rec.seq = rec.seq.ungap('-').ungap('.')
            return rec

        write((unalign(rec) for rec in parse(input[0], 'fasta')), output[0], 'fasta')
