rule pull_kingfisher:
    output:
        "./kingfisher_3.1.sif"
    shell:
        """
        if [ ! -f {output} ]; then
            apptainer pull --name {output} docker://wwood/kingfisher:0.3.1
        fi
        """
