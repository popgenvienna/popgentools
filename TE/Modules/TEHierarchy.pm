{
    package TEHierarchy;
    use strict;
    use warnings;
    use Data::Dumper;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(get_te_hierarchy_resolver get_TE_hierarchy get_te_hierarchy_translator);
    
    my $verbose=0;
    
    sub get_te_hierarchy_translator
    {
        my $file    =shift;
        my $from    =shift;
        my $to      =shift;
        my $teh=get_TE_hierarchy($file);
        
        my $translator={};
        while(my($key,$value)=each(%$teh))
        {
            my $tfrom   =$value->{$from};
            my $tto     =$value->{$to};
            die "Hierarchy level $from does not exist"unless exists($value->{$from});
            die "Hierarchy level $to does not exist" unless exists($value->{$to});
            if(exists($translator->{$tfrom}))
            {
                my $temp =$translator->{$tfrom};
                warn "Conflicting translation for $tfrom; either $temp or $tto"  unless $temp eq $tto;
            }
            $translator->{$tfrom}=$tto;
        }
        return $translator;
    }
    
    sub get_TE_hierarchy
    {

        
        my $te_hierarchyfile=shift;

        open my $ifh, "<", $te_hierarchyfile or die "Could not open TE file";
        
        # Read the header of the file: Must contain an entry "Insert" in the first column
        my $insertid="insert";
        my $headerline=<$ifh>;
        chomp $headerline;
        my @header=split /\t/,$headerline;
        my %temph=map {($_,1)}@header;
        die "First column of the TE-hierarchy file hast to be the insert; Eg.: column with heading '$insertid' has to be populated by entries like: FBti0018877, FBti0018884.." unless $header[0] eq $insertid;
        
        # Read the rest of the file ie: the individual TE entries
        my $colcount=@header;
        my $teh={};
        while(my $l=<$ifh>)
        {
            chomp $l;
            my @ar=split /\t/,$l;
            die "Length of the header does not fit with length of the entry: header length $colcount - column @ar" unless $colcount == scalar(@ar);
            my $eh={};
            
            for my $i(0..($colcount-1))
            {
                my $headercaption=$header[$i];
                my $entry=$ar[$i];
                $eh->{$headercaption}=$entry;
                
                my $id=$eh->{$insertid};
                $teh->{$id}=$eh;
            }
        }

        return $teh;
    }
    
    
    
    sub get_te_hierarchy_resolver
    {
        my $te_hierarchy_file=shift;
        my $target_level=shift;
        my $teh=get_TE_hierarchy($te_hierarchy_file);
        
        return sub
        {
            my $toresolve=shift;
            return undef unless(exists($teh->{$toresolve}{$target_level}));
            return $teh->{$toresolve}{$target_level};
        }
        
    }
    
}

1;

