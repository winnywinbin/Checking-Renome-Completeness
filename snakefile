configfile: "config.yaml"
singularity: "docker://snakemake/snakemake"
singularity: "40f9294e32eb0e6a0a3c357efd188f593f50c7bc095ccbc9040bb6369ce2e9b7"
rule all:
        input:
            "completeness.txt",
            "tree_met_data.png"

rule tRNA_scan:
        input:
                datas =expand("{file}", file=config["files"], allow_missing=True)
        params:
                userinput=config["soort"]
        output:
                dataEEN =temp("trnaSCANtussen_output.txt"),
                data = "trnaSCAN_output.txt"
        conda:
                "conda.yaml"
        shell:"""
                tRNAscan-SE {input.datas} {params.userinput} -o {output.dataEEN};
                cat trnaSCANtussen_output.txt | tail -n+4 | awk -F"\t" '{{print $1"\t"substr($1,1,9)"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' > {output.data}
              """

rule data_filteren:
        input:
                data = "trnaSCAN_output.txt"
        output:
                filtered_data="FilteredData.txt"
        conda:
                "conda.yaml"
        script:
                "filter.py"

rule r_op_filtered_data:
          input:
                    filtered_data = "FilteredData.txt"
          output:
                    data_counts = "data_counts.csv",
                    fisher_output_R = "data_fisher_results.txt"
            #conda:
            #        "conda.yaml"
          script:
                    "tRNAfisher.R"
                    
rule fisher_naar_excel:
        input:
                results="data_fisher_results.txt"
        output:
                excel="Table.xlsx"
        conda:
                "conda.yaml"
        script:
                "Results_naar_tabel.py"
                
rule cosine_similarity:
        input:
                data="FilteredData.txt"
        output:
                co_excel="gmatrix.xlsx"
        conda:
                "conda.yaml"
        script:
                "cossim.R"
                
rule hamming_distance:
        input:
                data="FilteredData.txt"
        output:
                ha_excel="hmatrix.xlsx"
        conda:
                "conda.yaml"
        script:
                "hamming.R"

rule GtRNA_filteren:
        input:
                data=expand("{GtRNA}.txt", GtRNA=config["GtRNAdbsets"], allow_missing=True)
        output:
                filtered_data=expand("Filtered_{GtRNA}.txt", GtRNA=config["GtRNAdbsets"])
        conda:
                "conda.yaml"
        script:
                "filter_GtRNA.py"

rule NMDS_plot_voorbereiding:
        input:
                filtered_archaea_data = "Filtered_archaea_GtRNAdb.txt",
                filtered_bacteria_data = "Filtered_bacteria_GtRNAdb.txt",
                filtered_eukaryotes_data = "Filtered_eukaryota_GtRNAdb.txt"
        output:
                filtered_bacteria_data = temp("Filtered_bac_GtRNAdb.txt"),
                filtered_eukaryotes_data = temp("Filtered_euk_GtRNAdb.txt"),
                Filtered_GtRNAdbsets = "Filtered_GtRNAdbsets.txt"
        shell:
                """
                cat {input.filtered_bacteria_data} | tail -n+2 > {output.filtered_bacteria_data};
                cat {input.filtered_eukaryotes_data} | tail -n+2 > {output.filtered_eukaryotes_data};
                cat {input.filtered_archaea_data} {output.filtered_bacteria_data} {output.filtered_eukaryotes_data} >> {output.Filtered_GtRNAdbsets}
                """

rule NMDS_plot:
        input:
                domains_combined = "Filtered_GtRNAdbsets.txt"
        output:
                nmdsplot = "nmdsplot.html"
        conda:
                "conda.yaml"
        script:
                "NMDS.R"

rule histogram_anticodons:
        input:
                filtered_archaea_data = "Filtered_archaea_GtRNAdb.txt"
        output:
                histogram = "histogram.png"
        script:
                "histogram.R"
                
rule heatmaps:
        input:
                filtered_archaea_data = "Filtered_archaea_GtRNAdb.txt"
        output:
                HeatmapAnticodonsNormalized = "HeatmapAnticodonsNormalized.html"
        conda:
                "conda.yaml"
        script:
                "Heatmap.R"

rule GC_percentage:
        input:
                expand("{file}", file=config["files"], allow_missing=True)
        output:
                "rapport_gc_perc.txt"
        conda:
                "conda.yaml"
        shell:"""
                bash genome_gc.sh
            """

rule completeness:
        input:
                gc_percentage="rapport_gc_perc.txt",
                filtered_archaea_data="Filtered_archaea_GtRNAdb.txt",
                filtered_data="FilteredData.txt"
        output:
                completeness = "completeness.txt"
        script:
                "completeness.R"
            

rule ncbi_genomes:
        output:
                ass_summary = "assembly_summary.txt",
                ftpdirpaths1 = temp("ftpdirpaths1"),
                ftpdirpaths = "ftpdirpaths",
                https="https.txt"
        params:
                userinput=config["soort_naam"]
        conda:
                "conda.yaml"
        shell:"""
                curl -O ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/{params.userinput}/{output.ass_summary};
                awk -F "\t" '$12=="Complete Genome" && $11=="latest"{{print $20}}' {output.ass_summary} > {output.ftpdirpaths1};
                awk 'BEGIN{{FS=OFS="/";filesuffix="genomic.fna.gz"}}{{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}}' {output.ftpdirpaths1} > {output.ftpdirpaths};
                sed -r "s/^(.{{0}})(.{{3}})/\https/" {output.ftpdirpaths} > https.txt;
                mkdir -p all/;
                cat https.txt | parallel -j 20 --verbose --progress "cd all/ && curl -O {{}}";
                """
                
rule unzip:
        input:
                downloaded = "assembly_summary.txt"
        output:
                ziped = "zipdone.txt"
        shell:"""
                gunzip all/*;
                touch zipdone.txt
                """

                
rule genomes_tRNAscan:
    input:
            "zipdone.txt",
            input_files = expand("{file}", file=config["files"], allow_missing=True)
    output:
            genomes_tRNA = "output_genomes_tRNA.txt"
    params:
            userinput=config["soort"]
    conda:
            "conda.yaml"
    shell:"""
            rm zipdone.txt;
            cp {input.input_files} all/;
            tRNAscan-SE all/* {params.userinput} -o {output.genomes_tRNA}
            """
            
            
rule boom_genereren:
    input:
            x = "output_genomes_tRNA.txt"
    output:
            accession = "accessions_genomes.txt",
            fasta_complete = "fasta_compleet.fa",
            domain_aln = "domain.aln",
            domain_ph = "domain.ph"
    shell:"""
            bash boom_code.sh
            """
            
            
rule safe_data:
    input:
            x="accessions_genomes.txt"
    output:
            tdata = "tdata.csv",
            clean_data = "tdata_clean.csv"
    shell:"""
            bash data_opslaan.sh
            """
            
rule tree:
    input:
            domain_ph = "domain.ph",
            clean_data = "tdata_clean.csv"
    output:
            tree = "tree_met_data.png"
    script:
            "tree_script.R"
            
            

rule NMDSwithsubmitted:
    input:
            domains_combined = "Filtered_GtRNAdbsets.txt",
            data="FilteredData.txt"
    output:
            nmds = "nmdsplotwithsubmitted.html"
    script:
            "NMDSwithsubmitted.R"
