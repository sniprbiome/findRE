mkdir ref
curl ftp://ftp.neb.com/pub/rebase/protein_gold_seqs.txt -o ref/REBASE_protein_gold_seqs.txt
./rebase2fasta ref/REBASE_protein_gold_seqs.txt > ref/REBASE_protein_gold_seqs.fa

curl ftp://ftp.neb.com/pub/rebase/protein_seqs.txt -o ref/REBASE_protein_seqs.txt
./rebase2fasta ref/REBASE_protein_seqs.txt > ref/REBASE_protein_seqs.fa
