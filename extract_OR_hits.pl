##Note that this assumes a blast output file and reference genome file uses as input!
my $blastfile=$ARGV[0];#take in files from command line
my $genomefile=$ARGV[1];
chomp $blastfile;chomp $genomefile;
my $outfile=$blastfile."_ORloci.fa";#create output file
open(RES, ">>$outfile");
open(GENOMEIN, "$genomefile");
open(GENOMEOUT, ">tmp_genome.fa");#remove whitespaces from current genome 
while(<GENOMEIN>){
    if($_=~m/\>/){
	print GENOMEOUT "\n".$_;
    }
    else{
	$_=~s/\s//g;
	print GENOMEOUT $_;
    }
}
close GENOMEIN;
open(IN, "$blastfile");
my @array=<IN>;
@array=split(/\>/, join('',@array));#loop through blast hits
close IN;
for(my $i=1;$i<scalar(@array);$i++){
    my $blasthit=$array[$i];
    my $header="";
    if($blasthit=~m/([\S]+)/){
	$header=">".$1;
    }
    my @hits=split(/[sS]core\s\=/,$blasthit);#loop through all subhits in blast
    foreach my $subhit(@hits){
	my $rev_comp=0;
        my $qpos1=0; my $qpos2=0; my $spos1=0; my $spos2=0;
	if($subhit=~m/Sbjct([0-9]+)/){
    $pos1=$1;
}
	if($subhit=~m/Query[\s]+([0-9]+)/){
    $qpos1=$1;
}
	while($subhit=~s/Query[\s]+[0-9]+[\s]+[A-Za-z\-]+[\s]+([0-9]+)//){
    $qpos2=$1;
}
	if($subhit=~m/Sbjct[\s]+([0-9]+)/){
    $spos1=$1;
}
	while($subhit=~s/Sbjct[\s]+[0-9]+[\s]+[A-Za-z\-]+[\s]+([0-9]+)//){#get all start and end coordinates for hits
	    $spos2=$1;
	}
	my $stringlen=0;
	if($spos1<$spos2){
	    $stringlen=$spos2-$spos1;
	}
	else{
	    $stringlen=$spos1-$spos2;
	}
	my $newpos1=0; my $newpos2=0;
	if($spos2<$spos1){#if hit is reverse compliment, make not and swap coordinates
	    my $tmp=$spos1;
	    $spos1=$spos2;
	    $spos2=$tmp;
	    $rev_comp=1;
	}
	$newpos1=$spos1-500;#go 500bp up and downstream
	$newpos2=$spos2+500;
	if($newpos1<=0){
	    $newpos1=0;
	}
	my $newlength=$newpos2-$newpos1;
	my $contig="";
	open(IN2, "tmp_genome.fa");#use coordinates to extract OR loci sequence
	my $switch=0;
	while(<IN2>){
	    my $seqheader="";
	    if($_=~m/(\>[\S]+)/){
		$seqheader=$1;
		if($seqheader eq $header){
		    print $seqheader." equals ".$header."\n";
		    $switch=1;
		}
	    }
	    else{
		$contig=$_;
		if($switch==0){
		}
		else{
		    $switch=0;
		    print $newpos1." ".$newlength."\n";
#		    print $contig;
		    my $sequence=substr($contig, $newpos1,$newlength);#if the hit is a reverse complement, mody the string sequence accordingly 
		    if($rev_comp==1){
			my $revcomp= reverse $sequence;
			$revcomp =~ tr/ATGCNatgcn/TACGNtacgn/;
			$sequence=$revcomp;
		    }
		    print RES ">".$header."_".$newpos1."_".$newpos2."\n".$sequence."\n";
		    $sequence="";
		}
	    }
	}
    }
}
