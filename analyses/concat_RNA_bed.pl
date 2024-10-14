#!/usr/bin/perl

# Concatenates gene expression data from RNA-seq samples into one master file

use File::Basename;

if ($#ARGV < 2){
    print stderr "USAGE: concat_RNA_bed.pl <bed file directory> <master bed file> <script directory>\n";
    exit(1);
}

$bed_dir = $ARGV[0];
$bed_file = $ARGV[1];
$scr_dir = $ARGV[2];

$tmp = "$bed_dir/tmp/";

if(!-d "$bed_dir"){
    die("$bed_dir does not exist\n");
}

if(!-d "$tmp"){
    `mkdir $tmp`;
}

opendir(DIR, $bed_dir);
@files = grep(/\.gene.count.bed$/, readdir(DIR));
closedir(DIR);
foreach $file (@files){
    $sample = basename($file, ".gene.count.bed");
    `cut -f7 "$bed_dir/$file" > "$tmp/$sample.count.bed.txt"`;
}

$file = @files[0];
`cut -f1-6 "$bed_dir/$file" > "$tmp/0.count.bed.txt"`;
`paste $tmp/* > $bed_file`;
`bash "$scr_dir/rename_chr.sh" $bed_file`;
`module load bedops/2.4.30; head -n 1 $bed_file > temp.txt; sort-bed $bed_file >> temp.txt && mv temp.txt $bed_file`;

`rm -r $tmp`;

exit 0;