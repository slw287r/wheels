# CollectSequencingArtifactMetrics

## How to install

* This tools requires <a href="https://github.com/samtools/samtools.git" target="_blank">samtools</a> and <a href="https://zlib.net" target="_blank">zlib</a> to compile. Please follow the link and install samtools first.

* Successful installation of samtools results in the following files (presumably in `/usr/local`):
    ```
    /usr/local/include/zlib.h
    /usr/local/include/bam.h
    /usr/local/include/htslib/*.h
    /usr/local/lib/libz.a
    /usr/local/lib/libbam.a
    /usr/local/lib/libhts.a
    ```

* Compile csam

	run `make` or
    ```
    gcc -o csam csam.c -lhts -lbam -lz
    ```
    or
    ```
    gcc -I/usr/local/include -L/usr/local/lib -o csam csam.c -lhts -lbam -lz
    ```
    incase `gcc` can't find the headers or libs

## How to use it

```
csam -i <bam> -r <ref> -p <prefix>
```

## Output files

```
${prefix}.pre_adapter_detail_metrics.txt
${prefix}.bait_bias_detail_metrics.txt
```

## Caveats

CONTEXT_SIZEs other than 1 are not yet supported in this implementation.

## Reference

<a href="https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.12.0/picard_analysis_artifacts_CollectSequencingArtifactMetrics.php" target="_blank">Picard CollectSequencingArtifactMetrics</a>
