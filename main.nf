#!/usr/bin/env nextflow

RES_DIR = params.resultsDir


process p01_process_data {
    conda "envs/process_data.yml"
    cache 'lenient'
    publishDir "$RES_DIR/01_process_data", mode: params.publishDirMode

    input:
        file "lib/*" from Channel.fromPath("lib/scio.py")
        // it would be better to include all input files explicitly.
        // However not that easy due to the fact that file identifier
        // is its parent folder.
        file 'data' from Channel.fromPath("data")
        file 'tables/*' from Channel.fromPath("tables/vanderburg_01_samples.csv")
        file 'process_data.py' from Channel.fromPath("analyses/01_process_counts.py")

    output:
        file "adata.h5ad" into process_data_adata

    """
    python process_data.py
    """
}


process p02_filter_data {
    def id = "02_filter_data"
    conda "envs/run_notebook.yml"
    cpus = 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'tables/*' from Channel.fromPath(
            "tables/{mitochondrial_genes,biomart,ribosomal_genes}.tsv"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from process_data_adata

    output:
        file "adata.h5ad" into filter_data_adata
        file "${id}.html" into filter_data_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="input_file=input_adata.h5ad output_file=adata.h5ad"
    """
}


process p03_correct_data {
    def id = "03_correct_data"
    conda "envs/run_notebook.yml"
    cpus = 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/{jupytertools,scio,scpp}.py").collect()
        file 'tables/*' from Channel.fromPath(
            "tables/{biomart,cell_cycle_regev}.tsv"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from filter_data_adata

    output:
        file "adata.h5ad" into correct_data_adata
        file "${id}.html" into correct_data_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="input_file=input_adata.h5ad output_file=adata.h5ad"
    """
}


process p04_annotate_cell_types {
    def id = "04_annotate_cell_types"
    conda "envs/run_notebook.yml"
    cpus = 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'tables/*' from Channel.fromPath(
            "tables/cell_type_markers.csv"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from correct_data_adata

    output:
        file "adata.h5ad" into annotate_cell_types_adata
        file "${id}.html" into annotate_cell_types_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="input_file=input_adata.h5ad output_file=adata.h5ad"
    """

}


process p05_prepare_de_analysis {
    def id = "05_prepare_de_analysis"
    conda "envs/run_notebook.yml"
    cpus 32
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'tables/*' from Channel.fromPath(
            "tables/cell_type_markers.csv"
        ).collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from annotate_cell_types_adata

    output:
        file "*.rda" into prepare_de_analysis_rdatas
        file "adata.h5ad" into prepare_de_analysis_adata,
            prepare_de_analysis_adata_2,
            prepare_de_analysis_adata_3,
            prepare_de_analysis_adata_4,
            prepare_de_analysis_adata_5,
            prepare_de_analysis_adata_6,
            prepare_de_analysis_adata_7
        file "adata_obs.tsv" into prepare_de_analysis_adata_obs
        file "${id}.html" into prepare_de_analysis_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --engine=papermill \
        --cpus=${task.cpus} \
        --params="input_file=input_adata.h5ad \
                  output_file=adata.h5ad \
                  output_file_obs=adata_obs.tsv \
                  cpus=${task.cpus} \
                  results_dir=."
    """
}


process p05_run_de_analysis {
    def id = "05-2_run_de_analysis"
    conda "envs/run_de.yml"
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    //found to have max performance for the problem.
    //with more CPUs the communication overhead between
    //BLAS workers slows down the analysis.
    cpus 6

    input:
        file input_data from prepare_de_analysis_rdatas.flatten()

    output:
        file "${input_data}.res.tsv" into run_de_analysis_results,
            run_de_analysis_results_2,
            run_de_analysis_results_3,
            run_de_analysis_results_4

    """
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
            MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus} \
            MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus} \
            MKL_THREADING_LAYER=GNU
    run_de.R ${input_data} ${input_data}.res.tsv --cpus=${task.cpus}
    """
}


process p10_analysis_t_nk {
    def id = "10_analysis_t_nk"
    conda "envs/run_notebook.yml"
    cpus 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from prepare_de_analysis_adata
        file 'input_adata_obs.tsv' from prepare_de_analysis_adata_obs
        file "*" from run_de_analysis_results.collect()

    output:
        file "${id}.html" into analysis_t_nk_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="input_adata=input_adata.h5ad \
                  input_obs=input_adata_obs.tsv \
                  input_de_res_dir=. \
                  cpus=${task.cpus}"
    """
}


process p50_prepare_analysis_cd39 {
    def id = "50_prepare_analysis_cd39"
    conda "envs/run_notebook.yml"
    cpus 8
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/jupytertools.py")
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from prepare_de_analysis_adata_7

    output:
        file "*.rda" into prepare_analysis_cd39_rdatas
        file "adata.h5ad" into prepare_analysis_cd39_adata,
            prepare_analysis_cd39_adata_2
        file "${id}.html" into prepare_analysis_cd39_html

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --engine=papermill \
        --cpus=${task.cpus} \
        --params="input_adata=input_adata.h5ad \
                  output_adata=adata.h5ad \
                  cpus=${task.cpus} \
                  results_dir=."
    """
}


process p50_run_de_analysis_cd39 {
    def id = "50-2_run_de_analysis_cd39"
    conda "envs/run_de.yml"
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    cpus 6

    input:
        file input_data from prepare_analysis_cd39_rdatas.flatten()

    output:
        file "${input_data}.res.tsv" into run_de_analysis_cd39_results

    """
    export OPENBLAS_NUM_THREADS=${task.cpus} OMP_NUM_THREADS=${task.cpus} \
            MKL_NUM_THREADS=${task.cpus} OMP_NUM_cpus=${task.cpus} \
            MKL_NUM_cpus=${task.cpus} OPENBLAS_NUM_cpus=${task.cpus} \
            MKL_THREADING_LAYER=GNU
    run_de.R ${input_data} ${input_data}.res.tsv --cpus=${task.cpus}
    """
}


process p51_analysis_cd39 {
    def id = "51_analysis_cd39"
    conda "envs/run_notebook3.yml"
    cpus 16
    publishDir "$RES_DIR/$id", mode: params.publishDirMode

    input:
        file 'lib/*' from Channel.fromPath("lib/{jupytertools,scenictools}.py").collect()
        file 'notebook.Rmd' from Channel.fromPath("analyses/${id}.Rmd")
        file 'input_adata.h5ad' from prepare_analysis_cd39_adata
        file 'adata_scenic.h5ad' from adata_scenic_4
        file "res_de_cd39/*" from run_de_analysis_cd39_results.collect()

    output:
        file "${id}.html" into analysis_cd39_html
        file "${id}.zip" into analysis_cd39_tables

    """
    reportsrender notebook.Rmd \
        ${id}.html \
        --cpus=${task.cpus} \
        --params="input_file=input_adata.h5ad \
                  adata_scenic=adata_scenic.h5ad \
                  input_de_res_dir_cd39=res_de_cd39 \
                  cpus=${task.cpus}"
    zip ${id}.zip *.xls*
    """
}



process deploy {
    conda "envs/reportsrender_index.yml"
    publishDir "${params.deployDir}", mode: "copy"
    executor "local"

    input:
        file "input/*" from Channel.from().mix(
            filter_data_html,
            correct_data_html,
            annotate_cell_types_html,
            prepare_de_analysis_html,
            analysis_t_nk_html,
            analysis_cd39_html,
            prepare_analysis_cd39_html,
            analysis_cd39_tables,
        ).collect()

    output:
        file "*.html"
        file "*.zip"

    """
    cp input/*.{html,zip} .
    """
}

