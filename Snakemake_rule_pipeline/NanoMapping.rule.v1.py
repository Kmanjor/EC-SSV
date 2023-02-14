################################# Quality control ######################
rule NanoFilt:
    input:
        fastq = Nano_PATH + "/raw/{sample}.fastq.gz",
    output:
        fastq = Nano_PATH + "/clean/{sample}.fastq.gz",
    threads:
        THREADS
    log:
        Nano_PATH + "/log/NanoFilt_{sample}.log"
    params:
        readtype = config["readtype"],
        minQuality = config["minQuality"],
        minLength = config["minLength"],
        headcrop = config["headcrop"],
        tailcrop = config["tailcrop"],
    run:
        cmd = "gunzip -c %s | NanoFilt --readtype %s --quality  %s --length %s --headcrop %s --tailcrop %s | gzip > %s 2>%s" % (input.fastq, params.readtype, params.minQuality, params.minLength, params.headcrop, params.tailcrop, output.fastq, log)
        print(cmd)
        os.system(cmd)


################################## mapping using minimap2 ######################
rule minimap2Align:
    input:
        x = Nano_PATH + "/clean/{sample}.fastq.gz",
    output:
        sam = temp(Nano_PATH + "/mapping/clean_minimap2/bam/{sample}.sam"),
    threads:
        THREADS * ThreadFold
    params:
        RefGenome = config["RefGenome"],
    log:
        Nano_PATH + "/log/clean_minimap2Align_{sample}.log"
    run:
        ### map-ont, map-pb, asm20, asm5
        shell("minimap2 --MD -a -x map-ont -t {threads} {params.RefGenome} {input.x} > {output.sam} 2>{log}")


rule SAM2BAM:
    input:
        sam = Nano_PATH + "/mapping/clean_minimap2/bam/{sample}.sam",
    output:
        bam = Nano_PATH + "/mapping/clean_minimap2/bam/{sample}.bam",
        bai = Nano_PATH + "/mapping/clean_minimap2/bam/{sample}.bam.bai",
    threads:
        THREADS * ThreadFold
    log:
        Nano_PATH + "/log/sam2bam_clean_minimap2_{sample}.log"
    run:
        shell("samtools view -@ {threads} -bS {input.sam} | samtools sort - -m 5G -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")


################################## mapping using ngmlr ######################
rule ngmlrAlign:
    input:
        fastq = Nano_PATH + "/clean/{sample}.fastq.gz",
    output:
        sam = temp(Nano_PATH + "/mapping/clean_ngmlr/{sample}.sam"),
    threads:
        THREADS * ThreadFold
    params:
        RefGenome = config["RefGenome"],
        ngmlr = config["ngmlr"],
    log:
        Nano_PATH + "/log/ngmlrAlign_{sample}.log"
    run:
        ### --rg-id MD
        shell("{params.ngmlr} --presets ont --rg-id MD --rg-sm {wildcards.sample} -t {threads} -r {params.RefGenome} -q {input.fastq} -o {output.sam} >{log} 2>&1")

rule SAM2BAM2:
    input:
        sam = Nano_PATH + "/mapping/clean_ngmlr/{sample}.sam",
    output:
        bam = Nano_PATH + "/mapping/clean_ngmlr/bam/{sample}.bam",
        bai = Nano_PATH + "/mapping/clean_ngmlr/bam/{sample}.bam.bai",
    threads:
        THREADS * ThreadFold
    log:
        Nano_PATH + "/log/sam2bam_clean_ngmlr_{sample}.log"
    run:
        shell("samtools view -@ {threads} -bS {input.sam} | samtools sort - -m 4G -@ {threads} -o {output.bam} > {log} 2>&1")
        shell("samtools index {output.bam}")
