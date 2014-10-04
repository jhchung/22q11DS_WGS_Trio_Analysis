
$in_file1 = shift; ###tbx1_gene ###
$in_file2 = shift; ###snp_indel ###
$out_file = shift;

die unless open (INFILE1, "$in_file1"); 
die unless open (INFILE2, "$in_file2");
die unless open (OUTFILE, ">$out_file");

undef @query_array;

while ($line=<INFILE1>)
{
	chomp $line;
	$key=$line;
	
	push(@query_array,$key);

	$query_hash{$key}=0;

}


while ($line=<INFILE2>)
{
	chomp $line;
	foreach $gene(@query_array)
	{
		$gene_t="\t".$gene."\t";
		$gene_q1='"'.$gene."(dist";
		$gene_q2=",".$gene."(dist";
		
		#print $gene_q1."\n";
		#print $gene_q2."\n";
		
		if ( $line=~/$gene_t/ || $line=~/\Q$gene_q1/ ||  $line=~/\Q$gene_q2/)
		{
			$query_hash{$gene}++;
		}
	
	}

}

foreach $key (@query_array)
{
	print OUTFILE $key."\t".$query_hash{$key}."\n";	
}

