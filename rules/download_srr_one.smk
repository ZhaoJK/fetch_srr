""" Download one SRR files"""

rule download_one_srr:
	"""Create fastqc report"""
	input: 
		run_id = srr_acc
	output: 
        '{input.run_id}_.fastq.gz'
    params:
        work_dir = config['WORK']['working_dir']
    shell:
		"""
            cd  {params.work_dir}
            kingfisher get -p  {input.run_id}  -m aws-http  aws-cp -f fastq.gz
		"""