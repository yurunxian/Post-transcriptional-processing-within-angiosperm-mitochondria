use strict;

my @gene = ("atp1","atp4","atp6","atp8","ccmB","ccmC","ccmFC","ccmFN","cob","cox1","cox2","cox3","matR","mttB","nad1","nad2","nad3","nad4","nad4L","nad5","nad6","nad7","nad9");
my @species = ("Amborella","Nymphaea","Lirio","AriG","AriF","AriM","Thottea","Asarum","Saruma","allium","phoenix","zea","oryza","nelumbo","malus","arabidopsis","helianthus","nicotiana");
my %species = ("Amborella"=>3,"Nymphaea"=>4,"Lirio"=>5,"Saruma"=>6,"Asarum"=>7,"Thottea"=>9,"AriF"=>10,"AriG"=>11,"AriM"=>12,"allium"=>13,"phoenix"=>14,"zea"=>15,"oryza"=>16,"nelumbo"=>17,"arabidopsis"=>18,"malus"=>19,"nicotiana"=>21,"helianthus"=>23);
my %species_exp_all = ("Ginkgo"=>0,"Amborella"=>0,"Nymphaea"=>0,"Lirio"=>0,"Saruma"=>0,"Asarum"=>0,"Thottea"=>0,"AriF"=>0,"AriG"=>0,"AriM"=>0,"allium"=>0,"phoenix"=>0,"zea"=>0,"oryza"=>0,"nelumbo"=>0,"arabidopsis"=>0,"malus"=>0,"daucus"=>0,"nicotiana"=>0,"rhazya"=>0,"helianthus"=>0);
my %species_length_all = ("Ginkgo"=>0,"Amborella"=>0,"Nymphaea"=>0,"Lirio"=>0,"Saruma"=>0,"Asarum"=>0,"Thottea"=>0,"AriF"=>0,"AriG"=>0,"AriM"=>0,"allium"=>0,"phoenix"=>0,"zea"=>0,"oryza"=>0,"nelumbo"=>0,"arabidopsis"=>0,"malus"=>0,"daucus"=>0,"nicotiana"=>0,"rhazya"=>0,"helianthus"=>0);

open (TAB,"<all_RNA_editing_site.txt");
my @all = <TAB>;
open (ANCESTOR,"<ancestral_RNA_editing_site.txt");
my @ancestor = <ANCESTOR>;
open (DENSITY,">exp.RES.density.txt");

foreach my $sp (@species){
	my %fasta = ();
	my $head;
	foreach my $gene (@gene){
		open (IN,"<../gene_align/$gene.fas") or die();
		while (<IN>) {
			chomp;
			if (/>(\S+)/) {$head = $1} 
			else {
				if ($head =~ /$sp/i) {
					my $tmp = $_;
					$tmp =~ s/-//gi;
					$fasta{$gene} .= $tmp;
				}
			}
		}
	}
	open (DEPTHH,"<../RNA_depth/depth_$sp.txt") or die();
	my @depth_all = <DEPTHH>;
	foreach my $gene (@gene){
		$species_length_all{$sp} += length($fasta{$gene});
		my @depth = grep {$_ =~ /$gene\s+/i} @depth_all;
		for (my $i = 200; $i < scalar(@depth)-200; $i++) {
			my @tmp3 = split /\s+/,$depth[$i];
			$species_exp_all{$sp} += $tmp3[2];
		}
	}
}

foreach my $gene (@gene){count($gene)}

sub count {
	my %species_RES = ("Ginkgo"=>0,"Amborella"=>0,"Nymphaea"=>0,"Lirio"=>0,"Saruma"=>0,"Asarum"=>0,"Thottea"=>0,"AriF"=>0,"AriG"=>0,"AriM"=>0,"allium"=>0,"phoenix"=>0,"zea"=>0,"oryza"=>0,"nelumbo"=>0,"arabidopsis"=>0,"malus"=>0,"daucus"=>0,"nicotiana"=>0,"rhazya"=>0,"helianthus"=>0);
	my %species_A = ("Ginkgo"=>0,"Amborella"=>0,"Nymphaea"=>0,"Lirio"=>0,"Saruma"=>0,"Asarum"=>0,"Thottea"=>0,"AriF"=>0,"AriG"=>0,"AriM"=>0,"allium"=>0,"phoenix"=>0,"zea"=>0,"oryza"=>0,"nelumbo"=>0,"arabidopsis"=>0,"malus"=>0,"daucus"=>0,"nicotiana"=>0,"rhazya"=>0,"helianthus"=>0);
	my %species_exp = ("Ginkgo"=>0,"Amborella"=>0,"Nymphaea"=>0,"Lirio"=>0,"Saruma"=>0,"Asarum"=>0,"Thottea"=>0,"AriF"=>0,"AriG"=>0,"AriM"=>0,"allium"=>0,"phoenix"=>0,"zea"=>0,"oryza"=>0,"nelumbo"=>0,"arabidopsis"=>0,"malus"=>0,"daucus"=>0,"nicotiana"=>0,"rhazya"=>0,"helianthus"=>0);
	my $gene = shift;
	print $gene,"\n";
	open (IN,"<../gene_align/$gene.fas") or die();
	my %fasta = ();
	my $head;
	while (<IN>) {
		chomp;
		if (/>(\S+)/) {$head = $1} else {$fasta{$head} .= $_}
	}
	my @tmp = grep {$_ =~ /$gene\s+/} @all;
	foreach my $i (@tmp){
		my @tmp2 = split /\s+/,$i;
		my $position = $tmp2[1];
		foreach my $sp (@species){
			if ($tmp2[$species{$sp}] eq "R") {
				if (substr($fasta{$sp},$position-1,1) =~ /C/i) {
					$species_RES{$sp}++;
					if (grep {$_ =~ /^$gene\t$position\t1/} @ancestor) {
						$species_A{$sp}++;
					}
				}
			}
		}
	}
	foreach my $sp (@species){$fasta{$sp} =~ s/-//gi}
	foreach my $sp (@species){
		open (DEPTH,"<../RNA_depth/depth_$sp.txt") or die();
		my @depth = grep {$_ =~ /$gene\s+/i} <DEPTH>;
		for (my $i = 200; $i < length($fasta{$sp})+200; $i++) {
			my @tmp3 = split /\s+/,$depth[$i];
			$species_exp{$sp} += $tmp3[2];
		}
		print DENSITY $gene,"\t$sp\t",length($fasta{$sp}),"\t",($species_exp{$sp}/(length($fasta{$sp})+1))/($species_exp_all{$sp}/$species_length_all{$sp}),"\t",$species_RES{$sp}*1000/(length($fasta{$sp})+1),"\t",1-($species_A{$sp}/($species_RES{$sp}+1)),"\n";
	}
}