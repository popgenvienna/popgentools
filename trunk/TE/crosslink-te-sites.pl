#!/usr/bin/perl -w
use strict;
use warnings;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use Getopt::Long;
use Pod::Usage;
use Utility;
use TEInsert;
use TEInsertUtility;
use TESiteUtility;
use TEHierarchy;

my $sitefile;
my $polynfile;
my $outputfile;
my $tehierfile;
my $tehiertargetlevel="order";
my $maxdist=250;
my $singlesiteshift=100;
my $mindist=74;
my $help=0;


GetOptions(
    "directional-insertions=s"  =>\$sitefile,
    "poly-n=s"                  =>\$polynfile,
    "output=s"                  =>\$outputfile,
    "te-hierarchy=s"            =>\$tehierfile,
    "te-hier-level=s"           =>\$tehiertargetlevel,
    "max-dist=i"                =>\$maxdist,
    "min-dist=i"                =>\$mindist,
    "single-site-shift=i"       =>\$singlesiteshift,
    "help"                      =>\$help
);
podusage(-verbose=>2) if $help;


print "Loading TE hierarchy...\n";
my $tehierres=get_te_hierarchy_translator($tehierfile,"family",$tehiertargetlevel);
print "Loading poly-N tracts...\n";
my $polynh=Utility::load_poly_n($polynfile);
print "Loading directional insertions...\n";
my $tesites=load_te_sites($sitefile);
print "Identifying overlapping directional insertions..\n";
$tesites=Utility::tag_overlapping_directionalinserts($tesites);
print "Crosslinking directional insertion sites..\n";
my $teinsertions=Utility::crosslink_directional_inserts($tesites,$polynh,$mindist,$maxdist,$singlesiteshift,$tehierres);
print "Writing TE Insertions into file...\n";
write_te_inserts($teinsertions,$outputfile);
print "Done\n";


exit;

{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use TEInsert;

    
    sub crosslink_directional_inserts
    {
        my $tesites = shift;
        my $polynh  = shift;
        my $mindist = shift;
        my $maxdist = shift;
        my $singlesiteshift=shift;
        my $tehierres=shift;
        
        my $makenoise=100;
        
        my $chromosomerep      = _get_chr_tesites($tesites);
        my $teinslist=[];
        while(my($chr,$tesitelist)=each(%$chromosomerep))
        {
            my $tesites=@$tesitelist;
            print "Crosslinking chromosome $chr; Number of TE sites: $tesites\n";
            $tesitelist = [sort {$a->{start}<=>$b->{start}} @$tesitelist];
            my $sitecounter=0;
            
            while(@$tesitelist)
            {
                my $primus=shift @$tesitelist;
                my $teins;
                if($primus->{insdir} eq "F")
                {
                    # start mate search!
                    my $prim_fam=$primus->{teid};
                    my $prim_end=$primus->{end};
                    
                    
                    my $leng=scalar(@$tesitelist);
                    my $truepartnerposition=-1;
                    my $dist_ncorrected=-1;

                    TEINSERT: for my $i (0..($leng-1))
                    {
                        my $candidate=$tesitelist->[$i];
                        unless($candidate){
                            my $lengcheck=@$tesitelist;
                            die "No candidates left for $i in $chr; @$tesitelist" ;
                        }
                        # next te candidate if different family or if not a reverse insertion
                        next unless $candidate->{insdir} eq "R";
                        next unless $candidate->{teid} eq $prim_fam;
                        
                        # get distance
                        my $cand_start=$candidate->{start};
                        $dist_ncorrected=_calculate_n_corrected_distance($chr,$prim_end,$cand_start,$polynh);
                        
                        # check the distance; try next if smaller than minimum distance; if larger than maximum distance than abort
                        next if $dist_ncorrected < $mindist;
                        last TEINSERT if $dist_ncorrected > $maxdist;
                        
                        $truepartnerposition=$i;
                        last TEINSERT;
                    }
                    
                    if($truepartnerposition!=-1)
                    {
                        my $truepartner=splice(@$tesitelist,$truepartnerposition,1);
                        $teins=_annotate_TEInsertion($primus,$truepartner,$singlesiteshift,$tehierres,"ncorrdist=$dist_ncorrected");
                    }
                    
                    unless($teins)
                    {
                        # truepartner position == -1 -> no partner
                        # no mate found -> is single -> store as TE Insert
                        $teins=_annotate_TEInsertion($primus,undef,$singlesiteshift,$tehierres,undef);
                    }
                }
                elsif($primus->{insdir} eq "R")
                {
                    # thats a single -> store as TE insert
                    $teins=_annotate_TEInsertion(undef,$primus,$singlesiteshift,$tehierres,undef);
                }
                else
                {
                    die "impossible";
                }
                
                die "No TEinsertion set"  unless $teins;
                $sitecounter++;
                print "Processed $sitecounter TE inserts \n" unless($sitecounter % $makenoise); 
                push @$teinslist, $teins;
                
            }
        }
        return $teinslist;
    }
    
    
    sub _annotate_TEInsertion
    {
        my $fwd=shift;
        my $rev=shift;
        my $singlesiteshift=shift;
        my $tehierres=shift;
        my $comment=shift;
        
        
        my $teins=undef;
        my($inspos,$teid,$order);
        #SITE: # chr, insdir, teid, support, start, end, comment, overlap
        # TEINsert: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        if($fwd and $rev)
        {
            $inspos=($fwd->{end}+$rev->{start})/2;
            $teid=$fwd->{teid};
            $order=$tehierres->{$teid};
            $teins=TEInsert->new($fwd->{chr},$inspos,$teid,$order,undef,$comment,$fwd->{start},$fwd->{end},undef,undef,$fwd->{overlap},$rev->{start},$rev->{end},undef,undef,$rev->{overlap});
        }
        elsif($fwd)
        {
            $inspos=$fwd->{end}+$singlesiteshift;
            $teid=$fwd->{teid};
            $order=$tehierres->{$teid};
            $teins=TEInsert->new($fwd->{chr},$inspos,$teid,$order,undef,$comment,$fwd->{start},$fwd->{end},undef,undef,$fwd->{overlap},undef,undef,undef,undef,undef);
            
        }
        elsif($rev)
        {
            $inspos=$rev->{start}-$singlesiteshift;
            $teid=$rev->{teid};
            $order=$tehierres->{$teid};
            $teins=TEInsert->new($rev->{chr},$inspos,$teid,$order,undef,$comment,undef,undef,undef,undef,undef,$rev->{start},$rev->{end},undef,undef,$rev->{overlap});
        }
        else
        {
            die "impossible";
        }
        
        return $teins;
    }
    
    sub _calculate_n_corrected_distance
    {
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $polynh=shift;
        
        # speed up
        return $end-$start if($end-$start >1_000_000);

        
        my $dist=0;
        $start++; $end--; # I want the distance in between, without the start and the end position
        for my $i ($start..$end)
        {
            $dist++ unless (exists($polynh->{$chr}{$i}));
        }
        return $dist;
    }
  
    
    sub _get_chr_tesites
    {
        my $tesites=shift;
        # chr, insdir, teid, support, start, end, comment
        
        my $chrrep={};
        foreach my $ts (@$tesites)
        {
            my $chr=$ts->{chr};
            $chrrep->{$chr}||=[];
            push @{$chrrep->{$chr}},$ts;
        }
        return $chrrep;
    }
    
    sub load_poly_n
    {
        my $file = shift;
        open my $ifh, "<",$file or die "Could not open poly-N file";
        
        my $ch={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            # YHet	polyNsearch	polyN	1	148	.	+	.	gene_id "poly_N_1";
            # YHet	polyNsearch	polyN	1731	1810	.	+	.	gene_id "poly_N_2";
            my($chr,undef,undef,$start,$end)=split /\t/,$l;
            for my $i ($start..$end)
            {
                $ch->{$chr}{$i}=1; 
            }
        }
        return $ch;
    }
    
    sub tag_overlapping_directionalinserts
    {
        my $teins=shift;
        
        my $overlapping={};
        my $chrhash={};
        foreach my $t (@$teins)
        {
            my $key=$t->key();
            my $chr=$t->{chr};
            my $insdir=$t->{insdir};
            my $start=$t->{start};
            my $end=$t->{end};
            
            for my $i($start..$end)
            {
               if(exists($chrhash->{$chr}{$insdir}{$i}))
               {
                    my $alreadypresent=$chrhash->{$chr}{$insdir}{$i};
                    $overlapping->{$alreadypresent}=1;
                    $overlapping->{$key}=1;
                
               }
               else
               {
                    $chrhash->{$chr}{$insdir}{$i}=$key;
               }
            }
        }
        
        foreach my $t (@$teins)
        {
            my $key=$t->key();
            my $overlap=0;
            $overlap=1 if(exists($overlapping->{$key}));
            $t->{overlap}=$overlap;
        }
        return $teins;
    }

}


    #"directional-insertions=s"  =>\$sitefile,
    #"poly-n=s"                  =>\$polynfile,
    #"output=s"                  =>\$outputfile,
    #"te-hierarchy=s"            =>\$tehierfile,
    #"te-hier-level=s"           =>\$tehiertargetlevel,
    #"max-dist=i"                =>\$maxdist,
    #"min-dist=i"                =>\$mindist,
    #"single-site-shift=i"       =>\$singlesiteshift
    
=head1 NAME

perl crosslink-te-sites.pl - Crosslinks forward and reverse insertions and outputs transposable element insertions

=head1 SYNOPSIS

 perl crosslink-te-sites.pl --direction-insertions te_insertsites.txt --poly-n polyn.gtf --output te_insertions.txt --te-hierarchy te_hier.txt --te-hier-level family --max-dist 250 --min-dist 70 --single-site-shift 100

=head1 OPTIONS

=over 4

=item B<--directional-insertions>

A list of TE insertion sites as identified by <identify-te-insertsites.pl>; Mandatory

=item B<--output>

The output file; will be a list of TE-insertions;

=item B<--te-hierarchy>

a custom hierarchy of the TE insertions; Must have one entry for every entry present in the sequence database used for repeatmaskin the genome; Mandatory

=item B<--te-hier-level>

the taret level of the TE analysis

=item B<--max-dist>

maximum distance between the forward and reverse insertion; poly-n stretches are not counted

=item B<--min-dist>

the minimum distance between the forward and the reverse insertion; poly-n stretches are not counted

=item B<--single-site-shift>

the exact position of a TE insertion can only be approximated for TE insertions where only the forward or only the reverse insertion was found.
For forward insertions the positon of the TE insertion is calculated as the end of the range + <--single-site-shift>;
For reverse insertions the position of the TE insertion is calculated as the start of the range  - <--single-site-shift>;

=item B<--help>

Display the help

=back



=cut
    


