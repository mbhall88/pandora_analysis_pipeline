def get_technology_param(wildcards):
    if wildcards.technology=="illumina":
        return "--illumina"
    else:
        return ""
