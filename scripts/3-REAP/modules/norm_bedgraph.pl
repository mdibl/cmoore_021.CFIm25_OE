#!/usr/bin/perl
#
##Normalized value of bedgraph file by total mapped read number (counts per 1 million reads)
use Getopt::Std;
my %args;
getopt("tldim", \%args);

my $tlRnum=$args{t}; #total read number, 
my $read_len=$args{l}; #read average length, invalid if -t is set
my $ind=$args{d}; #input directory, eg: "/Home/cocochen/analyze/dep_seq/othProject/LVM_seq/3tophat_Map/"
my $bedgraph_fnames_str=$args{i}; # input directory
my $fnames_2_minus_value=$args{m}; #bedgraph which will be convert all values to -1*original value

my @bedgraph_fnames=split(/ /,$bedgraph_fnames_str); #eg. "sham_a.m.bedgraph sham_a.p.bedgraph"
my %minus_val_files;
foreach my $bedgraph_f(split(/ /,$fnames_2_minus_value)){
	$minus_val_files{$bedgraph_f}=1;
}

my $help= qq(eg. norm_bedgraph.pl -t 9898819 -i "sham_a.m.bedgraph sham_a.p.bedgraph" -m "sham_a.m.bedgraph");
#command eg: norm_bedgraph.pl -l 70 -d "/Home/cocochen/analyze/dep_seq/othProject/LVM_seq/3tophat_Map/" -i "sham_a.m.bedgraph sham_a.p.bedgraph" -m "sham_a.m.bedgraph"
#provide -t is better than provide read average length

die "error, read length (-l) or total read number (-t) not defined\n$help\n" if (!$tlRnum && !$read_len);
die "error, input bedgraph file(s) (-i) not defined\n$help\n" if (!$bedgraph_fnames_str);


if(!$tlRnum){
	my $total_coverage=0;
	foreach my $bedgraph_f(@bedgraph_fnames){
		open (INPUT,"$ind$bedgraph_f") || die "error open bedgraph input file $ind$bedgraph_f\n";
		print "#open $ind$bedgraph_f\n";
		while(<INPUT>){
			chomp;
			next if !$_;
			my ($contig,$start,$end,$val)=split(/\t/); #start is 0-based, end is 1-based
			$total_coverage+=($end-$start)*$val;
		}
		close INPUT;
	}
	$tlRnum=sprintf("%.0f",$total_coverage/$read_len);
}

print "#Total mapped reads number = $tlRnum; ";
print "Normalized bedgraph coverage to read counts per million total read:\n";

foreach my $bedgraph_f(@bedgraph_fnames){
	my $scale_factor=1/$tlRnum * 1000000 * ($minus_val_files{$bedgraph_f}?-1:1);
	open (INPUT,"$ind$bedgraph_f") || die "error $ind$bedgraph_f\n";
	open (OUTPUT,">$ind$bedgraph_f.normalized") || die "error write $ind$bedgraph_f.normalized\n";
	print "#write $ind$bedgraph_f.normalized\n";
	while(<INPUT>){
		chomp;
		next if !$_;
		my ($contig,$start,$end,$val)=split(/\t/); #start is 0-based, end is 1-based
		print OUTPUT "$contig	$start	$end	". sprintf("%.1f", $val*$scale_factor). "\n";
	}
	close INPUT;
	close OUTPUT;
}
