use strict;

my @gene = ("atp1","atp4","atp6","atp8","atp9","ccmB","ccmC","ccmFC","ccmFN","cob","cox1","cox2","cox3","matR","mttB","nad1","nad2","nad3","nad4","nad4L","nad5","nad6","nad7","nad9","rpl10","rpl16","rpl2","rpl5","rps10","rps11","rps12","rps13","rps14","rps19","rps1","rps2","rps3","rps4","rps7","sdh3","sdh4");
my @species = ("Ginkgo","Amborella","Nymphaea","Lirio","AriG","AriF","AriM","Thottea","Asarum","Saruma","allium","phoenix","zea","oryza","nelumbo","glycine","malus","arabidopsis","gossypium","populus","cucumis","helianthus","nicotiana","rhazya","daucus","salvia");
my %species = ("Ginkgo"=>2,"Amborella"=>3,"Nymphaea"=>4,"Lirio"=>5,"Saruma"=>6,"Asarum"=>7,"Thottea"=>9,"AriF"=>10,"AriG"=>11,"AriM"=>12,"allium"=>13,"phoenix"=>14,"zea"=>15,"oryza"=>16,"nelumbo"=>17,"arabidopsis"=>18,"malus"=>19,"daucus"=>20,"nicotiana"=>21,"rhazya"=>22,"helianthus"=>23);

open (TAB,"<all_RNA_editing_site.txt");
my @all = <TAB>;
open (DENSITY,">RES.density.txt");
print DENSITY "gene\t";
foreach my $sp (@species){print DENSITY $sp,"\t"}
print DENSITY "\n";

foreach my $gene (@gene){count($gene)}

sub count {
	my %species_RES = ("Ginkgo"=>0,"Amborella"=>0,"Nymphaea"=>0,"Lirio"=>0,"Saruma"=>0,"Asarum"=>0,"Thottea"=>0,"AriF"=>0,"AriG"=>0,"AriM"=>0,"allium"=>0,"phoenix"=>0,"zea"=>0,"oryza"=>0,"nelumbo"=>0,"arabidopsis"=>0,"malus"=>0,"daucus"=>0,"nicotiana"=>0,"rhazya"=>0,"helianthus"=>0);
	my $gene = shift;
	print DENSITY "$gene\t";
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
				if (substr($fasta{$sp},$position-1,1) =~ /C/i) {$species_RES{$sp}++}
			}
		}
	}
	foreach my $sp (@species){$fasta{$sp} =~ s/-//gi}
	foreach my $sp (@species){print DENSITY $species_RES{$sp}," (",sprintf("%.4f",$species_RES{$sp}/(length($fasta{$sp})+1)),")\t"}
	print DENSITY "\n";
#	print DENSITY "$gene\t";
#	foreach my $sp (@species){print DENSITY $species_RES{$sp},"\t"}
#	print DENSITY "\n";
}