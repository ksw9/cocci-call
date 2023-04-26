process SummarizeRun {

  // Parse logs from TrimGalore, Kraken, and BWA
  
  label 'slurm'

  publishDir "${projectDir}/results", mode: "copy", pattern: "pipeline_run_summary_*.tsv"

  input:
  path reads_list
  path trimming_reports
  path kraken_reports
  path mapping_reports_immitis
  path coverage_stats_immitis
  path dup_metrics_immitis
  path mapping_reports_posadasii
  path coverage_stats_posadasii
  path dup_metrics_posadasii

  output:
  path "pipeline_run_summary_*.tsv"

  """
  python ${projectDir}/scripts/make_run_summary.py --reads_list_file ${reads_list}
  """

}