""" Download one SRR list """

rule download_srr_list:
	"""Create fastqc report"""
	input: 
		bioproj_id = bioproject_acc
	output: 
        'logs/{sample}_R1_fastqc.html'
    params:
        work_dir = config['WORK']['working_dir']
		"""
            cd  work_dir
            kingfisher get -p  {input.run_id}  -m aws-http  aws-cp -f fastq.gz
		"""