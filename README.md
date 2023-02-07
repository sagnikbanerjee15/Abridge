# ABRIDGE

Advancement in technology has enabled sequencing machines to produce vast amounts of genetic data, causing an increase in storage demands. Most genomic software utilizes read alignments for several purposes including transcriptome assembly and gene count estimation. Herein we present, ABRIDGE, a state-of-the-art compressor for SAM alignment files offering users both lossless and lossy compression options. This reference-based  file compressor achieves the best compression ratio among all compression software ensuring lower space demand and faster file transmission. Central to the software is a novel algorithm that retains non-redundant information. This new approach has allowed ABRIDGE to achieve a compression 16% higher than the second-best compressor for RNA-Seq reads and over 35% for DNA-Seq reads.  ABRIDGE also offers users the option to randomly access location without having to decompress the entire file. ABRIDGE is distributed under MIT license and can be obtained from GitHub and docker hub.  We anticipate that the user community will adopt ABRIDGE within their existing pipeline encouraging further research in this domain.

If you use this software in your research please cite the paper: 

# Installation

`abridge` is released through GitHub. Download the latest release from GitHub

```bash
wget https://github.com/sagnikbanerjee15/Abridge/archive/refs/tags/ABRIDGE_v1.0.0.tar.gz
tar -xvzf ABRIDGE_v1.0.0.tar.gz
cd ABRIDGE_v1.0.0
echo "export PATH=\$PATH:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

Download `abridge` directly from GitHub using `git clone`. Please note that this version is developmental and might contain bugs


```bash
git clone --recurse-submodules https://github.com/sagnikbanerjee15/Abridge.git
cd Abridge
make htslib
make install
echo "export PATH=\$PATH:$(pwd)/bin" >> ~/.bashrc
source ~/.bashrc
```

`abridge` will execute either `docker` or `singularity` depending on which software is available. 

# Generating alignments

Alignments can be generated using a number of different software. For more details on how to generate alignments please refer to https://github.com/alexdobin/STAR or https://github.com/DaehwanKimLab/hisat2

## Inputs to `abridge`

```bash
abridge --help

usage: abridge [-h] [-o OUTPUT_DIRECTORY]
               (-isam INPUTSAMFILENAMES | -iabr INPUTABRFILENAMES) -g GENOME
               [-cmp] [-dcmp] [-r] [-H] [-l LEVEL] [-ss] [-igqual]
               [-qual QUALITY] [-gsc] [-gm] [-gs] [-gu] [-sq] [-aq] [-q]
               [-n CPU] [-run_diagnostics] [-p POSITIONS] [-rp READ_PREFIX]
               [--keep_intermediate_error_files]
               [--error_directory ERROR_DIRECTORY]

Compress alignments for storage, decompress from compressed file, view alignments from random locations and generate coverages

optional arguments:
  -h, --help            show this help message and exit
  -isam INPUTSAMFILENAMES, --inputsamfilenames INPUTSAMFILENAMES
                        Enter the name of the alignment file you wish to compress. Alignments in SAM format only is expected. Ensure that the file is sorted by coordinate. Also, files must have the header section with the reference information available. You can compress only one file at a time.
  -iabr INPUTABRFILENAMES, --inputabrfilenames INPUTABRFILENAMES
                        Enter the name of the compressed alignment files you wish to merge. These files must be compressed using abridge. You can decompress only one file at a time.
  -cmp, --compress      Set this option if you wish to compress the alignment file
  -dcmp, --decompress   Set this option if you wish to decompress the alignment file
  -r, --random          Retrieve alignments from random locations
  -H, --header          Print only the header of reference sequences during decompression

Required arguments:
  -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                        Enter the name of the output directory. If nothing is specified then the compressed file will be put in the same location as the input samfile
  -g GENOME, --genome GENOME
                        Enter a single fasta file for the reference

Optional arguments:
  -l LEVEL, --level LEVEL
                        This can accept an integer from the set (1,2,3). If level is set to 1 then abridge will perform the fastest but low compression. abridge will use brotli to compress. Decompression will be fast. Setting level to 2 will prompt abridge to perform the medium level compression using 7z. Compression will take time but decompression will be fast. If level is set to 3 then abridge will perform the best compression using 7paq. Both compression and decompression will take average time to complete
  -ss, --ignore_scores  Request abrigde to store the quality scores and the alignment score
  -igqual, --ignore_quality_scores
                        Ignore all quality scores
  -qual QUALITY, --quality QUALITY
                        Enter dummy quality scores while decompressing
  -gsc, --ignore_soft_clippings
                        No soft clippings will be stored. Read will be trimmed down to only the portion which matched to nucleotides in the reference
  -gm, --ignore_mismatches
                        All mismatches will be ignored
  -gs, --ignore_sequence
                        No nucleotide sequence will be produced during decompression
  -gu, --ignore_unmapped_reads
                        Request abridge to discard all reads that are unmapped
  -sq, --save_all_quality_scores
                        Request abridge to save all quality scores
  -aq, --save_exact_quality_scores
                        Adjust quality scores for matched bases to achieve better encoding. For more details please check ...
  -q, --quiet           Prevent abridge from printing any log information. By default logging is enables
  -n CPU, --cpu CPU     Enter the number of CPU cores to be used. This option will be used during compression or decompression.
  -run_diagnostics, --run_diagnostics
                        abridge will run diagnostics on the cigar compression and decompression. It will exit on discovery of any discrepancies
  -p POSITIONS, --positions POSITIONS
                        Enter the position as chromosome:start-end from which reads will be retrieved
  -rp READ_PREFIX, --read_prefix READ_PREFIX
                        Enter a read prefix for decompression - valid only for random access
  --keep_intermediate_error_files, -kief
                        Set this argument if you wish to preserve the intermediate error files to assess time and memory usage. Default behaviour is to delete those
  --error_directory ERROR_DIRECTORY, -edir ERROR_DIRECTORY
                        Enter a directory where all error files will be stored. If nothing is specified then error files will be stored in the output directory
```

### Testing abridge

Run the following command if you wish to test abridge with the provided examples

```bash
CPU=96 # Enter the number of CPUs
run_abridge --test --output_directory $PWD/abridge_tests --cpu $CPU
```

This command will execute abridge with all possible parameters and with two different input files. There are more than 2000 commands that will execute. So please give it plenty of time and as many cores you can spare. To check if everything has been executed properly issue the command `ls -lh $PWD/abridge_tests/*error`. All the error files should be empty if the runs were successful.

# Compress

## Compress a single RNA-Seq file

`abridge` is developed to run on a single CPU. You can run multiple samples together by providing a single CPU core to each. Please note that `abridge` can handle a single genome file. If you have the sequences of the chromosomes on separate files, please merge them into one. During compression `abridge` stores the parameters. These parameters are reused during decompression.

### Lossless compression

```bash
abridge \
--compress \
--output_directory $PWD/abridge_compressed \
--inputsamfilenames <Name of input SAM file> \
--genome <Name of the genome file> \
--level 2 \
--save_all_quality_scores \
--save_exact_quality_scores \
--container_name compress \
1> $PWD/abridge_compressed.output \
2> $PWD/abridge_compressed.error 
```

In the lossless mode, `abridge` will store the entire information.

### Lossy compression

`abridge` offers a range of choices for lossy compression. We have presented a few use cases below.

1. Compressing a SAM file that will be later used for genome guided transcriptome assembly. There is no need to store soft-clips, mismatches, mapping score and quality scores.

   ```bash
   run_abridge \
   --compress \
   --output_directory $PWD/abridge_compressed \
   --inputsamfilenames <Name of input SAM file> \
   --genome <Name of the genome file> \
   --level 2 \
   --ignore_soft_clippings \
   --ignore_mismatches \
   --ignore_scores \
   --ignore_quality_scores \
   --container_name compress \
   1> $PWD/abridge_compressed.output \
   2> $PWD/abridge_compressed.error 
   ```

   

2. Compressing a SAM file which will be used for calling SNPs. For this case, most quality scores can be ignored except for those nucleotides that were a mismatch to the reference. 

   ```bash
   run_abridge \
   --compress \
   --output_directory $PWD/abridge_compressed \
   --inputsamfilenames <Name of input SAM file> \
   --genome <Name of the genome file> \
   --level 2 \
   --container_name compress \
   1> $PWD/abridge_compressed.output \
   2> $PWD/abridge_compressed.error 
   ```

   

# Decompress

Decompression is performed by providing the compressed file and the genome file as input. There is no need to explicitly state the parameters with which the compression was performed since `abridge` stores that information during compression.

```bash
run_abridge \
--decompress \
--output_directory $PWD/abridge_decompressed \
--inputabrfilenames <Name of the abridge compressed file> \
--genome <Name of the genome file> \
--container_name decompress \
1> $PWD/abridge_decompressed.output \
2> $PWD/abridge_decompressed.error 
```



# Frequently asked questions (FAQ)

1. I have several RNA-Seq samples aligned to a reference which is no longer available. What should I do?

   You must have the reference file available. Without the reference file `abdridge` will not be able to compress or decompress files. Please make sure you privoded the same reference both for compression and also for decompression.

2. I aligned my RNA-Seq samples to a reference without the MD tag. How do I generate MD tag for `abrigde` compression?

   `MD` flags can be generated using the following command:

   ```bash
   samtools calmd -bAr aln.bam > aln_baq.bam
   ```

   This command will add the **NM** and **MD** tags at the same time. For more information please check http://www.htslib.org/doc/samtools-calmd.html

3. How do I know that the compression is lossless?

   You can decompress the file to verify that the compression was indeed lossless. It should produce all the alignments along with the quality scores.

4. I need to view only the references and not alignments. What command should I execute?

   You need to execute `abridge` with the `header` argument

   ```bash
   run_abridge \
   --header \
   --inputabrfilenames <Name of the abridge compressed file> 
   ```

5. How do I retreive alignments to a particular chromosome?

   You can derive alignments to a particular chromosome by using the `random` command

   ```bash
   run_abridge \
   --random \
   --positions chr1:1-100000000 \
   --output_directory $PWD/abridge_decompressed \
   --inputabrfilenames <Name of the abridge compressed file> \
   --genome <Name of the genome file> \
   1> $PWD/abridge_decompressed.output \
   2> $PWD/abridge_decompressed.error 
   ```

   To extract all alignments to a particular chromosome, provide an arbitrarily large number - higher that the size of the chromosome.

6. I have only bamfiles but `abridge` needs samfiles to compress. What should I do?

   Convert the bamfile to samfiles using the following command

   ```bash
   samtools view -@ <number of CPU> <bamfilename> > <samfilename>
   ```

7. How do I retrieve multiple locations using random access?

   Option to retrieve alignments from multiple locations will be provided in the future



# Future upgrades

Here is a list of future upgrades to `abridge`

1. Retrieve alignments from multiple locations from compressed file

# License

License information can be accessed from https://github.com/sagnikbanerjee15/Abridge/blob/main/LICENSE

# Contact

Please list all issues at https://github.com/sagnikbanerjee15/Abridge/issues

For all other queries please email sagnikbanerjee15@gmail.com 
