nextflow.enable.types = true

process blastn {
    container "${projectDir}/images/powb-blast.sif"
    
    input:
        cols: String
        record( query: Path, asmid: String, subj: Path )
    
    output:
        record( asmid: asmid, blout: file("blastn_${asmid}.tsv") )
    
    script:
        """gunzip -c $subj |                    \
           blastn -evalue  0.01                 \
                  -outfmt  "6 ${cols}"          \
                  -query   $query               \
                  -subject -                    \
                  -out     blastn_${asmid}.tsv
        """
}

process filtres {
    container "${projectDir}/images/powb-python.sif"
    
    input:
        cols: String
        qcovs: Integer
        record( asmid: String, blout: Path )
    
    output:
        file("filter_${asmid}.tsv")
    
    script:
        """filter.py                       \
           --qcovs   ${qcovs}              \
           --cols   "${cols}"              \
           --input   $blout                \
           --output  filter_${asmid}.tsv
        """
}

process mergeres {
    input:
        paths: List<Path>
    
    output:
        file('merged.tsv')
    
    exec:
        def merged = task.workDir.resolve('merged.tsv')
        
        def header = paths.first().readLines().first()
        merged << header + '\n'
        
        paths.each { path ->
            def lines = path.readLines()
            if ( lines.size() > 1 ) {
                merged << lines.tail().join('\n') + '\n'
            }
        }
}

workflow {
    main:
        def cols = 'qseqid qlen sseqid length qcovs pident evalue sseq'
        
        col_ch = channel.value(cols)
        cov_ch = channel.value(params.qcovs)
        q_ch   = channel.fromPath(params.query)
        s_ch   = channel.fromPath(params.genomesDir + "/*.fna.gz")
        
        asm_s_ch = s_ch.map { path ->
            def asmid = (path.name =~ /(.+)\.fna/)[0][1]
            record( asmid: asmid, subj: path )
        }
        in_ch = q_ch.combine(asm_s_ch)
        .map { query, rec ->
            record( query: query, asmid: rec.asmid, subj: rec.subj )
        }
        
        blout_ch = blastn(col_ch, in_ch).filter { rec -> rec.blout.size() > 0 }
        fout_ch  = filtres(col_ch, cov_ch, blout_ch)
        
        merge_ch = mergeres( fout_ch.collect().map{ it -> it.toList() } )
}

