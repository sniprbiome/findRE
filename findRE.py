#!/usr/bin/env python3
import sys
import os
import io
import subprocess
from collections import defaultdict

# BLAST cut offs
MINIMUM_PIDENT = 80;
MINIMUM_COVERAGE = 0.8;

def main():
    SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))

    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("USAGE: findRE.py genome.fasta [all]")
        exit(1)
        
    input_genome = sys.argv[1]

    # Make sure input file exists
    if not os.path.isfile(input_genome):
        print(f"File does not exist: {input_genome}")
        exit(1)

    # Select which REBASE database to use
    if len(sys.argv) > 2 and sys.argv[2] == "all":
        rebase_fasta = SCRIPT_DIR+"/ref/REBASE_protein_seqs.fa"
        rebase_metadata = parse_rebase_metadata(SCRIPT_DIR+"/ref/REBASE_protein_seqs.txt")
    else:
        rebase_fasta = SCRIPT_DIR+"/ref/REBASE_protein_gold_seqs.fa"
        rebase_metadata = parse_rebase_metadata(SCRIPT_DIR+"/ref/REBASE_protein_gold_seqs.txt")

    # Build blast database for genome
    input_genome_blastdb = build_blast_database(input_genome)

    # Find matches to REFBASE database
    best_hits_per_region = defaultdict(list)
    run_tblastn(
        rebase_fasta,
        input_genome_blastdb,
        best_hits_per_region
    )

    # Print output table
    summarize_results(
        best_hits_per_region,
        rebase_metadata,
        input_genome_blastdb
    )
    
    # Clean up temporary files
    for ext in ["nhr", "nin", "nog", "nsd", "nsi", "nsq"]:
        os.remove(input_genome_blastdb + "." + ext)


def parse_rebase_metadata(infile):
    '''Parses the rebase annotation file'''

    rebase_data = defaultdict(dict)
    with open(infile, 'r') as rebase_file:
        for line in rebase_file:
            if line.startswith(">"):
                header_parts = line[1:].strip().split("\t")
                rebase_id = header_parts[0]
                for keyvalue_pair in header_parts:
                    parts = keyvalue_pair.split(":")
                    key = parts[0].strip()
                    value = ":".join(parts[1:]).strip()
                    rebase_data[rebase_id][key] = value
                    
    return rebase_data


def build_blast_database(input_fasta):
    '''Create a nucleotide blast database from a fasta file'''
    
    tmp_blast_db = os.getcwd()+"/tmpblast_"+os.path.basename(input_fasta)
    subprocess.run(
        [
            "makeblastdb",
            "-in", input_fasta,
            "-dbtype", "nucl",
            "-parse_seqids",
            "-logfile", "/dev/null",
            "-out", tmp_blast_db
        ]
    )
    return tmp_blast_db


def run_tblastn(rebase_fasta, genome, best_hits_per_region):
    '''Runs tblastn and collects the best hit for every region of the genome'''
    
    blast_header_str = "qaccver saccver evalue bitscore pident length qstart qend sstart send slen qlen";
    blast_header = blast_header_str.split(" ");

    proc = subprocess.Popen(
        [
            "tblastn",
            "-query", rebase_fasta,
            "-db", genome,
            "-evalue", "1e-50",
            "-outfmt", "6 "+blast_header_str
        ], stdout=subprocess.PIPE
    )

    for line in io.TextIOWrapper(proc.stdout, encoding="utf-8"):
        hit_parts = line.strip().split("\t")
        hit = dict(zip(blast_header, hit_parts))

        if float(hit["pident"]) < MINIMUM_PIDENT or float(hit["length"])/float(hit["qlen"]) < MINIMUM_COVERAGE:
            continue

        collect_best_hits(best_hits_per_region, hit)


def collect_best_hits(best_hits, hit):
    '''
    Adds 'hit' to 'best_hits' if it doesn't overlap a previous hit, or has a
    higher bitscore than the overlapping hit (which is then replaced)
    '''

    contig = hit.get("saccver")
    if contig not in best_hits:
        best_hits[contig].append(hit)
        return

    for idx, saved in enumerate(best_hits[contig]):
        if hit["sstart"] < saved["send"] and hit["send"] > saved["sstart"] and hit["bitscore"] > saved["bitscore"]:
            best_hits[contig][idx] = hit

            
def extract_subsequence(genome_blast_db, contig, start, end):
    '''
    Extracts a subsequence from a blast database, given the contig + 
    start and end coordinate. Picks the minus strand if start > end
    '''

    from_pos = min(start, end)
    to_pos = max(start, end)
    strand = "plus" if start < end else "minus"
    
    proc = subprocess.run(
        [
            "blastdbcmd",
            "-db", genome_blast_db,
            "-entry", contig,
            "-range", f"{from_pos}-{to_pos}",
            "-strand", strand,
            "-outfmt", "%s"
        ],
        capture_output = True,
        text = True
    )

    return str(proc.stdout).strip()

def summarize_results(hits, rebase_data, genome):
    '''Print output table'''
    
    print("\t".join(
        [
            "REBASE ID",
            "UniProt ID",
            "Genome location",
            "Protein identity %",
            "Coverage %",
            "Organism",
            "Enzyme type",
            "Recognition sequence",
            "Matching sequence"
        ]
    ))
        
    for contig, hits_in_contig in hits.items():
        for hit in hits_in_contig:
            rebase = rebase_data.get(hit["qaccver"])
            coverage = 100 * float(hit["length"]) / float(hit["qlen"])
            sequence = extract_subsequence(genome, hit["saccver"], int(hit["sstart"]), int(hit["send"]))
            print(
                "\t".join([
                    rebase.get("REBASE", "-"),
                    rebase.get("UniProt", "-"),
                    hit["saccver"] + ":" + hit["sstart"] + "-" + hit["send"],
                    hit["pident"],
                    str(round(coverage, 2)),
                    rebase.get("OrgName", "-"),
                    rebase.get("EnzType", "-"),
                    rebase.get("RecSeq", "-"),
                    sequence
                ])
            )

if __name__ == "__main__":
    main()
