../run_SPADES.sh \
    -i SRR8959866_sampled_1000read.fastq.gz \
    -o SRR8959866_test_ont \
    -p SRR8959866_test_ont \
    -d gottcha2_database_test/gottcha_db.species.fna \
    -t 4 \
    --js-external \
    --ont-error-rate 0.03 \
    --ont

../run_SPADES.sh \
    -i SRR1553609_NC_002549_1.mapped.fastq.gz \
    -o SRR1553609_clinical_test \
    -p SRR1553609_clinical_test \
    -d gottcha2_database_test/gottcha_db.species.fna \
    -t 4 \
	--js-external

../run_SPADES.sh \
    -1 SRR12689945_1_1k.fastq.gz \
    -2 SRR12689945_2_1k.fastq.gz \
    -o SRR12689945_test_PE \
    -p SRR12689945_test_PE \
    -d gottcha2_database_test/gottcha_db.species.fna \
    --min-depth 5 \
    -t 4 \
	--js-external
