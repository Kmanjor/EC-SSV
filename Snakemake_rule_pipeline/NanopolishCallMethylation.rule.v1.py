rule call_methylation:
    input:
        bam = Nano_PATH + "/Methylation/{sample}/output.sorted.bam",
        fastq = Nano_PATH + "/Methylation/{sample}/merge.fastq"
    output:
        tsv = Nano_PATH + "/Methylation/{sample}/methylation_calls.{chr}.tsv"
    params:
        ref = config["RefGenome"],
        nanopolish = config["nanopolish"]
    threads: 
        THREADS * ThreadFold
    log:
        Nano_PATH + "/Methylation/{sample}/methylation_calls.{chr}.log"
    run:
        cmd = "%s call-methylation --progress --verbose -t %s -r %s -b %s -g %s -w %s > %s 2>%s" % (params.nanopolish, threads, input.fastq, input.bam, params.ref, wildcards.chr, output.tsv, log)
        print(cmd)
        os.system(cmd)


rule frequency_methylation:
    input:
        tsv = Nano_PATH + "/Methylation/{sample}/methylation_calls.{chr}.tsv"
    output:
        tsv = Nano_PATH + "/Methylation/{sample}/methylation_frequency.{chr}.tsv"
    params:
        nanopolish_fre = config["nanopolish_fre"]
    threads: 
        THREADS
    log:
        Nano_PATH + "/Methylation/{sample}/methylation_frequency.{chr}.log"
    run:
        cmd = "python %s %s > %s 2>%s" % (params.nanopolish_fre, input.tsv, output.tsv, log)
        print(cmd)
        os.system(cmd)


rule mtsv2bedGraph:
    input:
        tsv = Nano_PATH + "/Methylation/{sample}/methylation_calls.{chr}.tsv"
    output:
        bed = Nano_PATH + "/Methylation/{sample}/methylation.{chr}.bed.gz",
        tbi = Nano_PATH + "/Methylation/{sample}/methylation.{chr}.bed.gz.tbi"
    params:
        ref = config["RefGenome"],
        mtsv2bedGraph = config["mtsv2bedGraph"],
        bgzip = config["bgzip"],
        tabix = config["tabix"]
    threads: 
        THREADS
    log:
        Nano_PATH + "/Methylation/{sample}/mtsv2bedGraph.{chr}.log"
    run:
        cmd = "python  %s -i %s -g %s  | sort -k1,1 -k2,2n | %s > %s && %s -p bed %s 2>%s" % (params.mtsv2bedGraph, input.tsv, params.ref, params.bgzip, output.bed, params.tabix, output.bed, log)
        print(cmd)
        os.system(cmd)


rule seperate_bam:
    input:
        bam = Nano_PATH + "/Methylation/{sample}/output.sorted.bam",
    output:
        bam = Nano_PATH + "/Methylation/{sample}/output.sorted.{chr}.bam",
        bai = Nano_PATH + "/Methylation/{sample}/output.sorted.{chr}.bam.bai"
    params:
        samtools = config["samtools"]
    threads: 
        THREADS * ThreadFold
    log:
        Nano_PATH + "/Methylation/{sample}/seperate_bam.{chr}.log"
    run:
        cmd = "%s view -h -@ %s -bS %s %s -o %s && %s index %s 2>%s" % (params.samtools, threads, input.bam, wildcards.chr, output.bam, params.samtools, output.bam, log)
        print(cmd)
        os.system(cmd)


rule convert_bam_for_methylation:
    input:
        bam = Nano_PATH + "/Methylation/{sample}/output.sorted.{chr}.bam",
        bed = Nano_PATH + "/Methylation/{sample}/methylation.{chr}.bed.gz"
    output:
        bam = Nano_PATH + "/Methylation/{sample}/converted.{chr}.bam",
        bai = Nano_PATH + "/Methylation/{sample}/converted.{chr}.bam.bai"
    params:
        ref = config["RefGenome"],
        convert_bam_for_methylation = config["convert_bam_for_methylation"],
        samtools = config["samtools"]
    threads: 
        THREADS * ThreadFold
    log:
        Nano_PATH + "/Methylation/{sample}/convert_bam_for_methylation.{chr}.log"
    run:
        cmd = "python %s --verbose -t %s -b %s -c %s  -f %s | %s sort -o %s && %s index %s 2>%s" % (params.convert_bam_for_methylation, threads, input.bam, input.bed, params.ref, params.samtools, output.bam, params.samtools, output.bam, log)
        print(cmd)
        os.system(cmd)