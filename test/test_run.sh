rm -rf *_output;

../run_SPADES.sh \
    -i SRR8959866_sampled_1000read.fastq.gz \
    -o SRR8959866_test_ont_output \
    -p SRR8959866_test_ont_output \
    -d gottcha2_database_test/gottcha_db.species.fna \
    -t 4 \
    --js-external \
    --ont-error-rate 0.03 \
    --ont

../run_SPADES.sh \
    -i SRR1553609_NC_002549_1.mapped.fastq.gz \
    -o SRR1553609_clinical_test_output \
    -p SRR1553609_clinical_test_output \
    -d gottcha2_database_test/gottcha_db.species.fna \
    -t 4 \
	--js-external

../run_SPADES.sh \
    -1 SRR12689945_1_1k.fastq.gz \
    -2 SRR12689945_2_1k.fastq.gz \
    -o SRR12689945_test_PE_output \
    -p SRR12689945_test_PE_output \
    -d gottcha2_database_test/gottcha_db.species.fna \
    --min-depth 5 \
    -t 4 \
	--js-external

../run_SPADES.sh \
    -i SRR8959866_sampled_1000read.fastq.gz \
    -o SRR8959866_test_ont_sqc_output \
    -p SRR8959866_test_ont_sqc_output \
    -d gottcha2_database_test/gottcha_db.species.fna \
    -t 4 \
    --ont-error-rate 0.03 \
    --ont

../run_SPADES.sh \
    -1 SRR12689945_1_1k.fastq.gz \
    -2 SRR12689945_2_1k.fastq.gz \
    -o SRR12689945_test_PE_sqc_output \
    -p SRR12689945_test_PE_sqc_output \
    --skip-qc \
    -d gottcha2_database_test/gottcha_db.species.fna \
    --min-depth 5 \
    -t 4

mkdir -p SRR8959866_stream_input/batch_001 SRR8959866_stream_input/batch_002
cp -p SRR8959866_sampled_1000read.fastq.gz SRR8959866_stream_input/batch_001/
cp -p SRR8959866_sampled_1000read.fastq.gz SRR8959866_stream_input/batch_002/

../stream_spades.py \
    --input-dir SRR8959866_stream_input \
    --outdir SRR8959866_test_stream_output \
    --prefix SRR8959866_test_stream_output \
    --db-path gottcha2_database_test/gottcha_db.species.fna \
    --once \
    --settle-seconds 0