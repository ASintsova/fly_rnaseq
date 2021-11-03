from pathlib import Path
import sys

def getFastq1(wildcards): # todo generalize to other data structures
    search_str = f'*{config["fq_fwd"]}'
    try:
        return str(list(Path(DATADIR).joinpath(wildcards.sample).rglob(search_str))[0])
    except IndexError:
        print(Path(DATADIR).joinpath(wildcards.sample))
        print(search_str)
        sys.exit(1)


if config['qc']:

    rule qc:
        input:
            fq1 = getFastq1,
            adapters = Path(config['adapters']),
            phix = Path(config['phix'])
        output:
            fq1_clean = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
            adapter_matched = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.adapter.matched.fq.gz',
            adapter_stats = OUTDIR /'clean_reads/{sample}/{sample}.adapter.stats',
            phix_matched = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.phix.matched.fq.gz',
            phix_stats = OUTDIR /'clean_reads/{sample}/{sample}.phix.stats',
            qc_failed = OUTDIR /'clean_reads/{sample}/removedreads/{sample}.qc.failed.fq.gz',
            qc_stats = OUTDIR /'clean_reads/{sample}/{sample}.qc.stats',
            marker = touch(OUTDIR /'clean_reads/{sample}/{sample}.qc.done')
        params:
            trimq = config['trimq'],
            maq = config['mapq'],
            minlen = config['minlen'],
            qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qout',
            qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qerr',
            scratch = 500,
            mem = 8000,
            time = 235
        conda:
            "envs/qc.yaml"
        log:
            log = OUTDIR /'logs/{sample}.qc.log'
        threads:
            8
        shell:
            "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t " 
            "in={input.fq1} out=stdout.fq outm={output.adapter_matched} "
            "refstats={output.adapter_stats} statscolumns=5 "
            "overwrite=t ref={input.adapters} "
            "ktrim=r k=23 mink=11 hdist=1  2>> {log.log} | "
            "bbduk.sh -Xmx1G usejni=t pigz=t bgzip=f "
            "overwrite=t interleaved=f in=stdin.fq out=stdout.fq "
            "outm={output.phix_matched} ref={input.phix} k=31 hdist=1 "
            "refstats={output.phix_stats} statscolumns=5 2>> {log.log}| "
            "bbduk.sh -Xmx1G pigz=t bgzip=f usejni=t  "
            "overwrite=t interleaved=f in=stdin.fq fastawrap=10000 "
            "out={output.fq1_clean} outm={output.qc_failed} "
            "minlength={params.minlen} qtrim=rl maq={params.maq} maxns=1  "
            "stats={output.qc_stats} statscolumns=5 "
            "trimq={params.trimq}  2>> {log.log};"

else:
    rule qc:
        input:
            fq1 = getFastq1,
        output:
            fq1_clean = OUTDIR /'clean_reads/{sample}/{sample}.1.fq.gz',
            marker = touch(OUTDIR /'clean_reads/{sample}/{sample}.qc.done')
        params:
            qoutfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qout',
            qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.qc.qerr',
            scratch = 500,
            mem = 8000,
            time = 235
        conda:
            "envs/qc.yaml"
        log:
            log = OUTDIR /'logs/qc/{sample}.qc.log'
        threads:
            8
        shell:
            "cp {input.fq1} {output.fq1_clean} "


rule fastqc_before:
    input: fq1 = getFastq1,

    output:
        marker = touch(OUTDIR/'fastqc/before/{sample}.fastqc.done'),
        #fqc = OUTDIR/'fastqc/before/{sample}_fastqc.html'
    params:
        outDir = OUTDIR/'fastqc/before',
        scratch = 1000,
        time = 100,
        mem = 8000,
        qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc.before.qerr',
        qoutfile =  lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc.before.qout'
    conda:
        'envs/qc.yaml'
    threads:
        4
    log:
        log = OUTDIR/'logs/qc/{sample}.fastqc.before.log'
    shell:
        '''
        fastqc {input.fq1}  -o {params.outDir} &> {log}
        
        '''

rule fastqc_after:
    input: fq1 = OUTDIR/'clean_reads/{sample}/{sample}.1.fq.gz',
    output:
        marker = touch(OUTDIR/'fastqc/after/{sample}.fastqc.done'),
        fqc = OUTDIR/'fastqc/after/{sample}.1_fastqc.html',

    params:
        outDir = OUTDIR/'fastqc/after',
        scratch = 1000,
        time = 100,
        mem = 8000,
        qerrfile = lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc.after.qerr',
        qoutfile =  lambda wildcards: OUTDIR /f'logs/qc/{wildcards.sample}.fastqc.after.qout'
    conda:
        'envs/qc.yaml'
    threads:
        4
    log:
        log = OUTDIR/'logs/qc/{sample}.fastqc.after.log'
    shell:
        '''
        fastqc {input.fq1}  -o {params.outDir} &> {log}
        '''

def multiqc_input():
    if config['fastqc'] == 'before':
        return [OUTDIR/f'fastqc/before/{sample}.fastqc.done' for sample in SAMPLES]
    elif config['fastqc'] == 'after':
        return [OUTDIR/f'fastqc/after/{sample}.fastqc.done' for sample in SAMPLES]
    elif config['fastqc'] == 'both':
        return [OUTDIR/f'fastqc/after/{sample}.fastqc.done' for sample in SAMPLES]+[OUTDIR/f'fastqc/before/{sample}.fastqc.done' for sample in SAMPLES]
    else:
        return ''


def multiqc_output():
    if config['fastqc'] == 'before':
        return OUTDIR/f'fastqc/before.multiqc.done'
    elif config['fastqc'] == 'after':
        return OUTDIR/f'fastqc/after.multiqc.done'
    elif config['fastqc'] == 'both':
        return OUTDIR/f'fastqc/both.multiqc.done'
    else:
        return ''

# # todo NOT WORKING right now
rule multiqc:
    input:
        multiqc_input()
    output:
        OUTDIR/'fastqc/multiqc_report.html',
        touch(multiqc_output())
    params:
        outdir=OUTDIR/'fastqc',
        scratch = 1000,
        time = 100,
        mem = 4000,
        qerrfile = lambda wildcards: OUTDIR /f'logs/qc/multiqc.qerr',
        qoutfile =  lambda wildcards: OUTDIR /f'logs/qc/multiqc.qout'
    conda:
        'envs/qc.yaml'
    shell: "multiqc {params.outdir} -o {params.outdir}"


