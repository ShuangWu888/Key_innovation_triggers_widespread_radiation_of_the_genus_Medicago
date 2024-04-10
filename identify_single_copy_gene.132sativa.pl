use strict;
use warnings;

my %cov=&sub_read_matrix("04.merge.gene_cov.132sativa.pl.coveper.out");
my %dp=&sub_read_matrix("04.merge.gene_cov.132sativa.pl.dpper.out");
my %len=&read_bed_len("Msa.gff.CDS.bed");

open (O,">05.cov0.5.dpper0.4-1.6.gene.pl.132sativa.out");
my @ind=sort keys %{$cov{MsaG000001}};#my$a=@ind;print "$a\n";exit;
print O "GeneID\tGene_len\tAlign_len\t",join("\t",@ind),"\n";
for my $k (sort keys %cov){
    next if $len{$k}<900;
    #my ($check1,$check2)=(0,0);
    my $line;#my @a;
    for my $ind (sort @ind){
        #$check1++ if $cov{$k}{$ind} >= 0.9;
        #$check2++ if (($dp{$k}{$ind} > 0.4) && ($dp{$k}{$ind} < 1.6));
	if ($cov{$k}{$ind} >= 0.5 && ($dp{$k}{$ind} > 0.4) && ($dp{$k}{$ind} < 1.6)){
	    #push(@a,$ind);
	    $line .= "$cov{$k}{$ind};$dp{$k}{$ind}\t";
	}
	else{$line .= "notpass\t";}
    }#if( @a==85){print "@a\n";exit;}
    #next if $check1 == 114;
    #next if $check2 == 114;#scalar(@ind)*0.8;
    print  O "$k\t$len{$k}\t$line\n";
}
close O;

#sub read_ortho{
 #   my %r;
  #  my ($in)=@_;
   # open (F,"$in")||die"$!";
    #while (<F>) {
     #   chomp;
      #  my @a=split(/\s+/,$_);
       # next if /^Cluster\s+/;
        #next if $a[2] < 300;
        #next if $a[2]/$a[1] < 0.5;
        #/Ore\|(\w+)/;
        #$r{$1}="$a[1]\t$a[2]";
    #}
    #close F;
    #return %r;
#}
sub sub_read_matrix{
    my %r;
    my ($in)=@_;
    my @id;
    open (F,"$in")||die"$!";
    while (<F>) {
        chomp;
	my @a=split(/\s+/,$_);
        if (/^\s+/){
            @id=@a;
            next;
        }
        for (my $i=1;$i<@a;$i++){
            $r{$a[0]}{$id[$i]}=$a[$i];
        }
    }
    close F;
    return %r;
}
sub read_name_info{
    my %r;
    my $nameinfo="Sample_location_132_info.txt";
    open (F,"$nameinfo")||die"$!";
    while (<F>) {
        chomp;
        /^(\S+)\s+(\S+)\s+(\S+)/;
        $r{$2}=$1;
    }
    close F;
    return %r;
}
sub read_bed_len{
    my %r;
    my ($tmpin)=@_;
    open (F,"$tmpin")||die"$!";
    while (<F>) {
	chomp;
	my @a=split(/\s+/,$_);
	$r{$a[3]} += $a[2]-$a[1]+1;
    }
    close F;
    return %r;
}
