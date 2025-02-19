## Using singularity containers

The available singularity containers are documented here: https://github.com/BF528/pipeline_containers

You may use these containers in your modules by specifying:

```
container 'ghcr.io/bf528/<name-of-tool>:latest'
```

For this project, you will need to make use of the following tools / containers:

ghcr.io/bf528/fastqc:latest
ghcr.io/bf528/star:latest
ghcr.io/bf528/verse:latest
ghcr.io/bf528/pandas:latest
ghcr.io/bf528/multiqc:latest
