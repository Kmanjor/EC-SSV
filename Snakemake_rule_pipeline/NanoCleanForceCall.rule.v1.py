##################################### force call SV using Sniffles (ngmlr) ###################################
rule ForceCallSniffles:
    input:
        bam = expand(Nano_PATH + "/mapping/clean_{tool}/bam/{sample}.bam",sample=SAMPLES,tool=TOOLS),
    output:
        tumor_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_sniffles/{pair}-tumor_{tool}.sniffles.vcf",
        normal_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_sniffles/{pair}-normal_{tool}.sniffles.vcf",
    threads:
        THREADS * ThreadFold
    params:
        min_support = config["min_support"],
        min_length  = config["min_length"],
        metafile = config["metafile"],
        sniffles = config["sniffles"],
    log:
        Nano_PATH + "/log/ForceCallSnifflesNano_{tool}_{pair}.log"
    run:
        normalFile, tumorFile = meta_target_file(input.bam, params.metafile, wildcards.pair)
        shell("{params.sniffles} --min_support {params.min_support} --min_length {params.min_length} -t {threads} -m {tumorFile} -v {output.tumor_vcf} --report_BND --genotype >{log} 2>&1")
        shell("{params.sniffles} --min_support {params.min_support} --min_length {params.min_length} -t {threads} -m {normalFile} -v {output.normal_vcf} --report_BND --genotype --Ivcf {output.tumor_vcf} >>{log} 2>&1")


rule combine_sniffles_forcecall:
    input:
        tumor_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_sniffles/{pair}-tumor_{tool}.sniffles.vcf",
        normal_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_sniffles/{pair}-normal_{tool}.sniffles.vcf",
    output:
        vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_sniffles/{pair}_{tool}.sniffles.vcf",
    threads:
        THREADS * ThreadFold
    params:
        combine_sniffles = MY_SRC_DIR + "/combine_sniffles_forcecall_output.py",
    log:
        Nano_PATH + "/log/combine_sniffles_forcecall_output_{pair}.log"
    run:
        shell("{params.combine_sniffles} -normal {input.normal_vcf} -tumor {input.tumor_vcf} -out {output.vcf} >{log} 2>&1")


##################################### force call SV using cuteSV (ngmlr) ###################################
rule ForceCallcuteSV:
    input:
        bam = expand(Nano_PATH + "/mapping/clean_{tool}/bam/{sample}.bam",sample=SAMPLES,tool=TOOLS),
    output:
        tumor_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_cuteSV/{pair}-tumor_{tool}.cuteSV.vcf",
        normal_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_cuteSV/{pair}-normal_{tool}.cuteSV.vcf",
    threads:
        THREADS * ThreadFold
    params:
        metafile = config["metafile"],
        RefGenome = config['RefGenome'],
        tmp_dir = Nano_PATH + "/clean_SVCall/force_call/{tool}_cuteSV/{pair}_temp_cuteSV",
        cuteSV = config['cuteSV'],
        min_support = config["min_support"],
        min_length  = config["min_length"],
    log:
        Nano_PATH + "/log/ForceCallcuteSVNano_{tool}_{pair}.log"
    run:
        shell("mkdir {params.tmp_dir}")
        normalFile, tumorFile = meta_target_file(input.bam, params.metafile, wildcards.pair)
        shell("{params.cuteSV} {tumorFile} {params.RefGenome} {output.tumor_vcf} {params.tmp_dir} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --min_support {params.min_support} --min_size {params.min_length} --genotype -t {threads} --max_size 3000000000 >{log} 2>&1")
        shell("{params.cuteSV} {normalFile} {params.RefGenome} {output.normal_vcf} {params.tmp_dir} --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 --min_support {params.min_support} --min_size {params.min_length} --genotype -t {threads} -Ivcf {output.tumor_vcf} --max_size 3000000000 >>{log} 2>&1")


rule combine_cuteSV_forcecall:
    input:
        tumor_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_cuteSV/{pair}-tumor_{tool}.cuteSV.vcf",
        normal_vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_cuteSV/{pair}-normal_{tool}.cuteSV.vcf",
    output:
        vcf = Nano_PATH + "/clean_SVCall/force_call/{tool}_cuteSV/{pair}_{tool}.cuteSV.vcf",
    threads:
        THREADS * ThreadFold
    params:
        combine_cuteSV = MY_SRC_DIR + "/combine_cuteSV_forcecall_output.py",
    log:
        Nano_PATH + "/log/combine_cuteSV_forcecall_output_{pair}.log"
    run:
        shell("{params.combine_cuteSV} -normal {input.normal_vcf} -tumor {input.tumor_vcf} -out {output.vcf} >{log} 2>&1")


##################################### force call SV filter  ###################################
rule forcecall_filt1plus:
    input:
        minimap2_sniffles_vcf = Nano_PATH + "/clean_SVCall/force_call/minimap2_sniffles/{pair}_minimap2.sniffles.vcf",
        minimap2_cuteSV_vcf = Nano_PATH + "/clean_SVCall/force_call/minimap2_cuteSV/{pair}_minimap2.cuteSV.vcf",
        ngmlr_sniffles_vcf = Nano_PATH + "/clean_SVCall/force_call/ngmlr_sniffles/{pair}_ngmlr.sniffles.vcf",
        ngmlr_cuteSV_vcf = Nano_PATH + "/clean_SVCall/force_call/ngmlr_cuteSV/{pair}_ngmlr.cuteSV.vcf",
    output:
        minimap2_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_minimap2.sniffles.filt1plus.vcf",
        minimap2_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_minimap2.cuteSV.filt1plus.vcf",
        ngmlr_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_ngmlr.sniffles.filt1plus.vcf",
        ngmlr_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_ngmlr.cuteSV.filt1plus.vcf",
    threads:
        THREADS * ThreadFold
    params:
        somaticSV_filt1plus = MY_SRC_DIR + "/somaticSV_filt1plus.py",
    log:
        Nano_PATH + "/log/somaticSV_filt1plus_{pair}.log"
    run:
        shell("python {params.somaticSV_filt1plus} -raw {input.minimap2_sniffles_vcf} -new {output.minimap2_sniffles_vcf} 2>{log}")  
        shell("python {params.somaticSV_filt1plus} -raw {input.minimap2_cuteSV_vcf} -new {output.minimap2_cuteSV_vcf} 2>>{log}")  
        shell("python {params.somaticSV_filt1plus} -raw {input.ngmlr_sniffles_vcf} -new {output.ngmlr_sniffles_vcf} 2>>{log}")
        shell("python {params.somaticSV_filt1plus} -raw {input.ngmlr_cuteSV_vcf} -new {output.ngmlr_cuteSV_vcf} 2>>{log}")


rule nanomonsv_filt1plus:
    input:
        vcf = Nano_PATH + "/SV_plus/nanomonsv_plus/{pair}-tumor.nanomonsv.result.ln.vcf",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_minimap2.nanomonsv.filt1plus.vcf",
    threads:
        THREADS
    params:
        somaticSV_nanomonsv_filt1plus = MY_SRC_DIR + "/somaticSV_nanomonsv_filt1plus.py",
    log:
        Nano_PATH + "/log/somaticSV_nanomonsv_filt1plus_{pair}.log"
    run:
        shell("python {params.somaticSV_nanomonsv_filt1plus} -raw {input.vcf} -new {output.vcf} 2>{log}")