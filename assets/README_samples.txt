Renamed target-demultiplexed fastq files per Nathan naming convention are saved at:
/home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI

Sample set1: SIKI1 (2022_08)
    SIKI1,2,3: 3 replicates across 5 targets: ep300a, ep300b, etv2, h3h3d, hey2
    SIKI5,6,7: 3 replicates for a single target: ep300b 

1-1: Raw from Nathan:
/Users/kaihu/Dropbox/Kai_Hu_MCCB_Nathan_Lawson/data/SIKI1_2022_08_04/raw/PacBio/2022_8_4_SIKI_firstRun/demultiplexed_samples

1-2: Demultiplexed by targets (SIKI1,2,3):
/Users/kaihu/Projects/Nathan/work/sikiproject/sikipipe/test/demultiplex/data/res/ep300a 
/Users/kaihu/Projects/Nathan/work/sikiproject/sikipipe/test/demultiplex/data/res/ep300b
/Users/kaihu/Projects/Nathan/work/sikiproject/sikipipe/test/demultiplex/data/res/etv2
/Users/kaihu/Projects/Nathan/work/sikiproject/sikipipe/test/demultiplex/data/res/h3f3d
/Users/kaihu/Projects/Nathan/work/sikiproject/sikipipe/test/demultiplex/data/res/hey2
/Users/kaihu/Projects/Nathan/work/sikiproject/sikipipe/test/demultiplex/data/res/ep300b_siki567 (not being used by Nathan, but still include in the SRA submission folder)

1-3: Processed by sikipipe:

- Ref: 
/Users/kaihu/Dropbox/Kai_Hu_MCCB_Nathan_Lawson/2022_08_CRISPR_Nathan_Lawson_shared_Kai/Kai_work/0_notes.docx
/Users/kaihu/Dropbox/Kai_Hu_MCCB_Nathan_Lawson/data/SIKI1_2022_08_04


Sample set2: SIKI2 (2023_11)





## Prepare a single folder on SCI that contain renamed files for both SIKI1 and SIKI2
cd /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI1_2022_08/demultiplex_targets
cd ep300a
for file in bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].fastq.gz; do
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/ep300a_1_${file//--/_}" | bash
done

cd ../ep300b
for file in bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].fastq.gz; do
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/ep300b_1_${file//--/_}" | bash
done

cd ../etv2
for file in bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].fastq.gz; do
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/etv2_1_${file//--/_}" | bash
done

cd ../hey2
for file in bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].fastq.gz; do
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/hey2_1_${file//--/_}" | bash
done

cd ../h3f3d
for file in bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].fastq.gz; do
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/h3f3d_1_${file//--/_}" | bash
done

cd ../ep300b_siki567
for file in bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].fastq.gz; do
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/ep300b_1_${file//--/_}" | bash
done

cd /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI2_2023_11/demultiplex_targets
cd ep300a 
for file in demultiplex.bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].hifi_reads.fastq.gz; do
    new_file="${file#demultiplex.}"
    new_file="${new_file//hifi_reads./}"
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/ep300a_2_${new_file//--/_}" | bash
done

cd ../h3f3d
for file in demultiplex.bc[0-9][0-9][0-9][0-9]--bc[0-9][0-9][0-9][0-9].hifi_reads.fastq.gz; do
    new_file="${file#demultiplex.}"
    new_file="${new_file//hifi_reads./}"
    echo "cp $file /home/kai.hu-umw/pi/nathan.lawson-umw/data/SIKI/h3f3d_2_${new_file//--/_}" | bash
done