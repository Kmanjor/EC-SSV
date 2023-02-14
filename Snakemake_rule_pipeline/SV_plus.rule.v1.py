
##################################### force call SV filter (minimap2 & ngmlr) ###################################
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


################################# reformat for jasmine ######################
rule Reformat:
    input:
        minimap2_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_minimap2.sniffles.filt1plus.vcf",
        minimap2_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_minimap2.cuteSV.filt1plus.vcf",
        ngmlr_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_ngmlr.sniffles.filt1plus.vcf",
        ngmlr_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_ngmlr.cuteSV.filt1plus.vcf",
        vcf = Nano_PATH + "/SV_plus/filt1plus/{pair}_minimap2.nanomonsv.filt1plus.vcf",
    output:
        minimap2_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_minimap2.sniffles.filt1plus.vcf",
        minimap2_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_minimap2.cuteSV.filt1plus.vcf",
        ngmlr_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_ngmlr.sniffles.filt1plus.vcf",
        ngmlr_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_ngmlr.cuteSV.filt1plus.vcf",
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_minimap2.nanomonsv.filt1plus.vcf",
    threads:
        THREADS
    params:
        reformat_sniffles_for_jasmine = MY_SRC_DIR + "/reformat_sniffles_for_jasmine.py",
        reformat_cuteSV_for_jasmine = MY_SRC_DIR + "/reformat_cuteSV_for_jasmine.py",
        reformat_nanomonsv_for_jasmine = MY_SRC_DIR + "/reformat_nanomonsv_for_jasmine.py",
    log:
        Nano_PATH + "/log/Reformat_for_jasmine_{pair}.log"
    run:
        shell("python {params.reformat_sniffles_for_jasmine} -raw {input.minimap2_sniffles_vcf} -new {output.minimap2_sniffles_vcf} 2>{log}")
        shell("python {params.reformat_cuteSV_for_jasmine} -raw {input.minimap2_cuteSV_vcf} -new {output.minimap2_cuteSV_vcf} 2>{log}")
        shell("python {params.reformat_sniffles_for_jasmine} -raw {input.ngmlr_sniffles_vcf} -new {output.ngmlr_sniffles_vcf} 2>>{log}")        
        shell("python {params.reformat_cuteSV_for_jasmine} -raw {input.ngmlr_cuteSV_vcf} -new {output.ngmlr_cuteSV_vcf} 2>>{log}")
        shell("python {params.reformat_nanomonsv_for_jasmine} -raw {input.vcf} -new {output.vcf} 2>>{log}")


rule jasmine_5merge:
    input:
        minimap2_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_minimap2.sniffles.filt1plus.vcf",
        minimap2_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_minimap2.cuteSV.filt1plus.vcf",
        ngmlr_sniffles_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_ngmlr.sniffles.filt1plus.vcf",
        ngmlr_cuteSV_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_ngmlr.cuteSV.filt1plus.vcf",
        nanomonsv_vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_format/{pair}_minimap2.nanomonsv.filt1plus.vcf",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/{pair}.filt1plus.5merge.vcf",
    threads:
        THREADS * ThreadFold
    params:
        RefGenome = config["RefGenome"],
        samtools = config["samtools"],
        tmp_dir = Nano_PATH + "/jasmine/tmp"
    log:
        Nano_PATH + "/log/jasmine_5_{pair}.log"
    run:
        files = ','.join([input.minimap2_sniffles_vcf,input.minimap2_cuteSV_vcf,input.ngmlr_sniffles_vcf,input.ngmlr_cuteSV_vcf,input.nanomonsv_vcf])
        cmd = "source activate jasminesv && jasmine file_list=%s out_file=%s --comma_filelist genome_file=%s samtools_path=%s threads=%s out_dir=%s --dup_to_ins --output_genotypes --ignore_strand 2>%s" % (files, output.vcf, params.RefGenome, params.samtools, threads, params.tmp_dir, log)
        print(cmd)
        os.system(cmd)

rule jasmine_5merge_three:
    input:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/{pair}.filt1plus.5merge.vcf",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/{pair}.filt1plus.5merge_three.vcf",
    threads:
        THREADS * ThreadFold
    params:
        RefGenome = config["RefGenome"],
        samtools = config["samtools"],
        tmp_dir = Nano_PATH + "/jasmine/tmp"
    log:
        Nano_PATH + "/log/jasmine_5_{pair}.log"
    run:
        shell("grep -v 'SUPP=1' {input.vcf} | grep -v 'SUPP=2' > {output.vcf} 2>{log}")  


rule jasmine_5merge_three_lengthflt:
    input:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/{pair}.filt1plus.5merge_three.vcf",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/{pair}.filt1plus.5merge_three.flt.vcf",
    threads:
        THREADS
    params:
        jasmine_length_filt = MY_SRC_DIR + "/jasmine_length_filt.py",
    log:
        Nano_PATH + "/log/jasmine_5merge_three_lengthflt_{pair}.log"
    run:
        shell("python {params.jasmine_length_filt} -raw {input.vcf} -new {output.vcf} 2>{log}")


rule jasmine_5merge_three_lengthflt_2bed:
    input:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/{pair}.filt1plus.5merge_three.flt.vcf",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/bed/{pair}.filt1plus.5merge_three.flt.bed",
    threads:
        THREADS
    params:
        jasmine2bed4depth = MY_SRC_DIR + "/jasmine2bed4depth.py",
    log:
        Nano_PATH + "/log/jasmine_5merge_three_lengthflt_2bed_{pair}.log"
    run:
        shell("python {params.jasmine2bed4depth} -raw {input.vcf} -new {output.vcf} 2>{log}")

rule jasmine_5merge_three_lengthflt_depth:
    input:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/bed/{pair}.filt1plus.5merge_three.flt.bed",
        bam = Nano_PATH + "/mapping/clean_minimap2/bam/{pair}-normal.bam",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/bed/{pair}.filt1plus.5merge_three.flt.depth",
    threads:
        THREADS
    params:
        samtools = config["samtools"],
    log:
        Nano_PATH + "/log/jasmine_5merge_three_lengthflt_2depth_{pair}.log"
    run:
        shell("samtools depth -b {input.vcf} {input.bam} > {output.vcf} 2>{log}")
        

rule jasmine_5merge_three_lengthflt_depthfilt:
    input:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/{pair}.filt1plus.5merge_three.flt.vcf",
        depth = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/bed/{pair}.filt1plus.5merge_three.flt.depth",
    output:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/depth_flt/{pair}.filt1plus.5merge_three.flt.depth.vcf",
    threads:
        THREADS
    params:
        jasmine_depth_filt = MY_SRC_DIR + "/depth_filt.py",
    log:
        Nano_PATH + "/log/jasmine_5merge_three_lengthflt_depthfilt_{pair}.log"
    run:
        shell("python {params.jasmine_depth_filt} -raw {input.vcf} -depth {input.depth} -new {output.vcf} 2>{log}")

     
rule jasmine_5merge_three_lengthflt_depthfilt_getigv:
    input:
        vcf = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/depth_flt/{pair}.filt1plus.5merge_three.flt.depth.vcf",
        t_bam = Nano_PATH + "/mapping/clean_minimap2/bam/{pair}-tumor.bam",
        n_bam = Nano_PATH + "/mapping/clean_minimap2/bam/{pair}-normal.bam",
    output:
        point = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/depth_flt/IGVsnapshot/{pair}_point_igv.batch",
        sv = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/depth_flt/IGVsnapshot/{pair}_sv_igv.batch",
    threads:
        THREADS
    params:
        RefGenome = config["RefGenome"],
        IGV_image = Nano_PATH + "/SV_plus/filt1plus/jasmine_merge/3in5/lenth_flt/depth_flt/IGV_image",
        jasmine2IGV = MY_SRC_DIR + "/jasmine2IGV.py",
    log:
        Nano_PATH + "/log/jasmine_5merge_three_lengthflt_depthfilt_{pair}.log"
    run:
        bam_files = ','.join([input.t_bam,input.n_bam])
        shell("python {params.jasmine2IGV} -b {input.vcf} -p {output.point} -o {output.sv} -d {params.IGV_image} -g {params.RefGenome} -m {bam_files} -s {wildcards.pair} -l 800 2>{log}")