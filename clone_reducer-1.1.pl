#!/usr/bin/env perl -w
use strict;
use Bio::TreeIO;
use Getopt::Long;
use Pod::Usage;

my $man = 0;
my $help = 0;
my $alignment;
my $tree;
my $perlen = 0.5;
my $bootstrapvalue = 50;
my $consensus;
my $longest;
GetOptions('help|?' => \$help, man => \$man, "alignment=s" => \$alignment, "tree=s" => \$tree, "percent_length=f" => \$perlen, "bootstrap=i" => \$bootstrapvalue, 'consensus' => \$consensus, 'longest' => \$longest) or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


=head1 NAME

clone_reducer

=head1 SYNOPSIS

perl clone_reducer.pl -alignment file -tree file [options] 

=head1 OPTIONS

	-alignment 			Alignment file in FASTA format used to produce gene tree
	-tree 				Newick tree file created from alignment
	-percent_length 		Percent in decimal format of unaligned sequence to alignment length [Default: 0.5]
	-bootstrap 			Minimum bootstrap value for clade to be considered for consensus sequence [Default: 50]
	-consensus 			Flag to output consensus sequence of clones
	-longest 			Flag to output longest sequence of clones
	-help 				Brief help message
	-man 				Full documentation

=cut

my (%alignment, %condensed_clones);
$alignment =~ /^(.+)\.\w+$/;
my $prefix = $1;
open my $bad_seq_file, ">", $prefix . "_clones_condensed.fsa";
open my $logfile, ">", $prefix . "_logfile_condensed_clones.fsa";

&read_fasta($alignment, $perlen);

&clone_remover($tree, $bootstrapvalue);

close $bad_seq_file;
close $logfile;
my $outfile = $prefix . "_clonesremoved.fsa";
open my $cleaned_seq_file, '>', $outfile;
for my $sample_ids (sort keys %alignment){
	my $seq = $alignment{$sample_ids};
	#$seq =~ s/-//g;
	unless ($seq){
		print "Look at $sample_ids\n";
	}
	print $cleaned_seq_file ">$sample_ids\n$seq\n";
}
close $cleaned_seq_file;
	

########################
sub read_fasta {
	open my $file, "<", $_[0];
	my ($seqid, $cloneid, $seq);
	while(<$file>){
		chomp;
		if(/>/){
			if($seq){
				my $original = $seq;
				$original =~ s/-//g;
				if(length($original)/length($seq) >= $_[1]){
					$alignment{$seqid}= $seq;
				}
				else{
					print $bad_seq_file ">$seqid" . "_shortseq\n$original\n";
				}
				$seq=();
			}
			$seqid = substr($_,1);
			#$cloneid = reverse(substr(reverse($_), 0, index(reverse($_), "_"))); 
		}
		else{
			$seq .= $_;
		}
	}
	my $original = $seq;
	$original =~ s/-//g;
	if(length($original)/length($seq) >= $_[1]){
		$alignment{$seqid} = $seq;
	}
	else{
		print $bad_seq_file ">$seqid" . "_shortseq\n$original";
	}
	
}
########################
sub clone_remover{
	my $treeio = new Bio::TreeIO(-format => "newick", -file => "$_[0]", -internal_node_id => 'bootstrap');
	my %used_taxa;
	while( my $tree = $treeio->next_tree ) {
		my @taxa = $tree->get_leaf_nodes;
		TAXA: for my $taxon (@taxa){
			my $taxon_id = $taxon->id;
			if(exists $used_taxa{$taxon_id}){
				next TAXA;
			}
			$used_taxa{$taxon}=1;
			my $sample = substr($taxon_id, 0, index($taxon_id, "_"));
			my (@descendents,@old_descendents);
			my $ancestor = $taxon->ancestor;
			push(@old_descendents, $taxon_id);
			my $bs = $ancestor->bootstrap;
			unless ($bs){
				next TAXA;
			}
			if ($bs >= $_[1]){
				
				my @temp_descendents = $ancestor->get_all_Descendents;
				for my $temp_descend (@temp_descendents){
					if($temp_descend->is_Leaf){
						$temp_descend=$temp_descend->id;
						push(@descendents, $temp_descend);
					}
				}	
				my $same_copy = 0;
				while ($same_copy == 0){
					for my $desc_taxon (@descendents){
							if ($desc_taxon !~ /$sample/){
								if(scalar @old_descendents > 1){
									if($consensus){
										&consensus_seq(\@old_descendents);
									}
									if($longest){
										&longest_seq(\@old_descendents);
									}
									for my $old_taxa (@old_descendents){
										$used_taxa{$old_taxa}=1;
									}
									$same_copy = 1;
									next TAXA;
								}
								else{
									next TAXA;
								}
							
						}
					}
					@old_descendents=@descendents;
					$ancestor = $ancestor->ancestor;
					$bs = $ancestor->bootstrap;
					if ($bs >= $_[1]){
						@temp_descendents = $ancestor->get_all_Descendents;
						@descendents=();
						for my $temp_descend (@temp_descendents){
							if($temp_descend->is_Leaf){
								$temp_descend=$temp_descend->id;
								push(@descendents, $temp_descend);
							}
						}	
					}
					else{
						if($consensus){
							&consensus_seq(\@old_descendents);
						}
						if($longest){
							&longest_seq(\@old_descendents);
						}
						for my $old_taxa (@old_descendents){
							$used_taxa{$old_taxa}=1;
						}
						$same_copy = 1;
						next TAXA;
					}
						
					
				}
			}
			else{
				next TAXA;
			}
		}
	}
}
########################
sub consensus_seq{
	my %consensus;
	my $sample_id = substr(@{$_[0]}[0], 0, index(@{$_[0]}[0], "_")); 
	my $species = reverse(substr(reverse(substr(@{$_[0]}[0], index(@{$_[0]}[0], "_")+1)), index(reverse(substr(@{$_[0]}[0], index(@{$_[0]}[0], "_")+1)), "_")+1));
	print $logfile "Consenus sequence made for:\t$sample_id\nClones used:\t";
	CLONE: for my $clone (@{$_[0]}){
		#RAxML adds an underscore to taxa names.  Must be removed to id match.
		if($clone =~ /_$/){
				$clone = substr($clone, 0, -1);
		}
		if (exists $alignment{$clone}){
			$sample_id = substr($clone, 0, index($clone, "_"));
			my $old_seq = $alignment{$clone};
			$old_seq =~ s/-//g;
			print $bad_seq_file ">$clone" . "_condensedclone\n$old_seq\n";
			my $temp_cloneid = substr($clone, 0, -1);
			my $clone_id_marker= reverse(substr(reverse($temp_cloneid), 0, index(reverse($temp_cloneid), "_")));
			print $logfile "$clone_id_marker\t";
			my $seq = $alignment{$clone};
			my @current = split(//,$seq);
       		my $count = 0;
        	for my $sq (@current) {
        		$count++;
            		$consensus{$clone}{$count} = $sq;
    		}
    		
    		delete $alignment{$clone};
    	}
    }
    print $logfile "\n";
    
    my ($consensus_count, $consensus_id);
    if (exists $condensed_clones{$sample_id}){
    	$consensus_count = scalar keys %{$condensed_clones{$sample_id}};
    	$consensus_count++;
    	$condensed_clones{$sample_id}{$consensus_count}="consensus";
    	$consensus_id = $sample_id . "_$species" . "_$consensus_count" . "consensus";
    	print $logfile "New consenus sequence name:\t$consensus_id\n";
    }
    else{
    	$consensus_count = 1; 
    	$consensus_id = $sample_id . "_$species" . "_$consensus_count" . "consensus";
    	$condensed_clones{$sample_id}{$consensus_count}="consensus";
		print $logfile "New consenus sequence name:\t$consensus_id\n";
    }
    if(scalar keys %consensus == 0){
		return;
    }
    
    my $seqsize;
    
    for my $temp_id (sort keys %consensus){
    	$seqsize = scalar keys %{$consensus{$temp_id}};
    }

	if(scalar keys %consensus == 1){
    	my $con_seq;
		for my $con_id (sort keys %consensus){
			for (my $j = 1; $j <= $seqsize; $j++){
				$con_seq .= $consensus{$con_id}{$j};
			}
			$alignment{$con_id}=$con_seq;
		}
		return;
    }
	
	my $newseq;
	for(my $i = 1; $i <= $seqsize; $i++){
    	my ($A_count, $T_count, $C_count, $G_count, $sp_count) = qw(0 0 0 0 0);
    	for my $con_id (sort keys %consensus){
        	if($consensus{$con_id}{$i} eq "A"){
        		$A_count++;
            }
            if($consensus{$con_id}{$i} eq "T"){
             	$T_count++;
            }
  	       	if($consensus{$con_id}{$i} eq "C"){
            	$C_count++;
            }
        	if($consensus{$con_id}{$i} eq "G"){
            	$G_count++;
            }
  	     	if($consensus{$con_id}{$i} eq "-" || $consensus{$con_id}{$i} eq "N"){
           		$sp_count++;
            }
		}
			
			
        my @counts;
        push(@counts, $A_count);
        push(@counts, $T_count);
        push(@counts, $C_count);
        push(@counts, $G_count);
        push(@counts, $sp_count);
        my $max = 0;
        my $seqval=-1;
        my $j=0;
        for my $counts (@counts){
        	if ($counts > $max){
            	$seqval = $j;
                $max = $counts;
            }
            $j++;
        }
        if ($seqval == 0){
        	$newseq .= "A";
        }
        if ($seqval == 1){
        	$newseq .= "T";
	    }
        if ($seqval == 2){
            $newseq .= "C";
        }
        if ($seqval == 3){
            $newseq .= "G";
        }
        if ($seqval == 4){
        	$newseq .= "-";
        }
    }
    
    $alignment{$consensus_id}=$newseq;
}
		
########################
sub longest_seq {
	my %consensus;
	my $sample_id = substr(@{$_[0]}[0], 0, index(@{$_[0]}[0], "_")); 
	my $species = reverse(substr(reverse(substr(@{$_[0]}[0], index(@{$_[0]}[0], "_")+1)), index(reverse(substr(@{$_[0]}[0], index(@{$_[0]}[0], "_")+1)), "_")+1));
	print $logfile "Largest sequence chosen for:\t$sample_id\nClones used:\t";

	my $longid;
	my $long_len=0;
	my $longseq;
	CLONE: for my $clone (@{$_[0]}){
		#RAxML adds an underscore to taxa names.  Must be removed to id match.
		if($clone =~ /_$/){
				$clone = substr($clone, 0, -1);
		}
		if (exists $alignment{$clone}){
			$sample_id = substr($clone, 0, index($clone, "_"));
			my $old_seq = $alignment{$clone};
			$old_seq =~ s/-//g;
			print $bad_seq_file ">$clone" . "_condensedclone\n$old_seq\n";
			my $temp_cloneid = substr($clone, 0, -1);
			my $clone_id_marker= reverse(substr(reverse($temp_cloneid), 0, index(reverse($temp_cloneid), "_")));
			print $logfile "$clone_id_marker\t";
			my $seq = $alignment{$clone};
			my $cur_len = length($seq);
       		if($cur_len > $long_len){
       				$long_len = $cur_len;
       				$longid = $clone;
       				$longseq = $alignment{$clone};
       		}
        	  	
    		delete $alignment{$clone};
    	}
    }
    print $logfile "\n";
	    	
    print $logfile "Longest Sequence Used:\t$longid\n";
    $alignment{$longid}=$;
}
	
			
