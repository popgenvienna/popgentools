{
    package TEInsert;
    use strict;
    use warnings;
    

    # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
    sub new {
        my $class = shift;
        
        
        #chromosome	position	RF/F/R	teid	order	fbid	popfreq		comment
        #rstart rend	popfreq cov	abs	pres	overlap	
        #rstart	rend	popfreq cov     abs	pres	overlap
        
        # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        my $chr=shift;
        my $inspos=shift;
        my $teid=shift;
        my $order=shift;
        my $fbid=shift;
        my $comment=shift;
        my $frstart=shift;
        my $frend=shift;
        my $fpres=shift;
        my $fabs=shift;
        my $foverlap=shift;
        my $rrstart=shift;
        my $rrend=shift;
        my $rpres=shift;
        my $rabs=shift;
        my $roverlap=shift;
        
        my($fcov,$rcov,$sitesupport)=(undef,undef,undef);
        my($popfreq,$rpopfreq,$fpopfreq)=(undef,undef,undef);
        $fcov=$fabs+$fpres if(defined($fabs) and defined($fpres));
        $rcov=$rabs+$rpres if(defined($rabs) and defined($rpres));
        $fpopfreq=$fpres/$fcov if $fcov;
        $rpopfreq=$rpres/$rcov if $rcov;
        
        # BASE: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        # DERIVED: fcov, rcov, sitesupport, popfreq, rpopfreq, fpopfreq
        if($frstart and $rrstart)
        {
            $popfreq=($fpopfreq+$rpopfreq)/2 if(defined($fpopfreq) and defined($rpopfreq));
            $sitesupport="FR";
        }
        elsif($frstart)
        {
            $popfreq=$fpopfreq;
            $sitesupport="F";
            
        }
        elsif($rrstart)
        {
            $popfreq=$rpopfreq;
            $sitesupport="R";
        }
        else
        {
            die "can not creat a TEInsert object withou a valid TESite";
        }
        
        
        my $entry={
            chr         =>$chr,
            inspos      =>$inspos,
            sitesupport =>$sitesupport,
            teid        =>$teid,
            popfreq     =>$popfreq,
            order       =>$order,
            fbid        =>$fbid,
            comment     =>$comment,
            frstart     =>$frstart,
            frend       =>$frend,
            fpopfreq    =>$fpopfreq,
            fcov        =>$fcov,
            fpres       =>$fpres,
            fabs        =>$fabs,
            foverlap    =>$foverlap,
            rrstart     =>$rrstart,
            rrend       =>$rrend,
            rpopfreq    =>$rpopfreq,
            rcov        =>$rcov,
            rpres       =>$rpres,
            rabs        =>$rabs,
            roverlap    =>$roverlap
        };
        # INTERFACE:
        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
        # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
        return bless $entry, __PACKAGE__;
    }
    
    sub tostring
    {
        my $self=shift;
        
        # get variables which may need overwriting;
        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
        # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
        my $popfreq = $self->{popfreq};
        $popfreq="-" unless(defined($popfreq));
        my $order   = $self->{order}    || "-";
        my $fbid    = $self->{fbid}     || "-";
        my $comment = $self->{comment}  || "-";
        my $frstart = $self->{frstart};     my $rrstart     =$self->{rrstart};
        my $frend   = $self->{frend};       my $rrend       =$self->{rrend};
        my $fpopfreq= $self->{fpopfreq};    my $rpopfreq    =$self->{rpopfreq};
        my $fcov    = $self->{fcov};        my $rcov        =$self->{rcov};
        my $fpres   = $self->{fpres};       my $rpres       =$self->{rpres};
        my $fabs    = $self->{fabs};        my $rabs        =$self->{rabs};
        my $foverlap= $self->{foverlap};    my $roverlap    =$self->{roverlap};
        ($fpopfreq,$fcov,$fpres,$fabs)=("-","-","-","-") unless ($fcov);
        ($rpopfreq,$rcov,$rpres,$rabs)=("-","-","-","-") unless ($rcov);
        ($frstart,$frend,$fpopfreq,$fcov,$fpres,$fabs,$foverlap)=("-","-","-","-","-","-","-") unless (defined($frstart));
        ($rrstart,$rrend,$rpopfreq,$rcov,$rpres,$rabs,$roverlap)=("-","-","-","-","-","-","-") unless (defined($rrstart));


        my @ar=($self->{chr},$self->{inspos},$self->{sitesupport},$self->{teid},$popfreq,$order,$fbid,$comment,$frstart,$frend,$fpopfreq,$fcov,$fpres,$fabs,$foverlap,$rrstart,$rrend,$rpopfreq,$rcov,$rpres,$rabs,$roverlap);
        return join ("\t",@ar);
    }
    
    sub clone
    {
        my $self=shift;
        
        # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        my $insert=TEInsert->new($self->{chr},$self->{inspos},$self->{teid},$self->{order},$self->{fbid},$self->{comment},$self->{frstart},$self->{frend},$self->{fpres},$self->{fabs},
                                $self->{foverlap},$self->{rrstart},$self->{rrend},$self->{rpres},$self->{rabs},$self->{roverlap});
        
        return $insert;
    }
    
    sub subclone_polymorphism
    {
        my $self=shift;
        my $fpres=shift;
        my $fabs=shift;
        my $rpres=shift;
        my $rabs=shift;
                # INT: chr, inspos, teid, order, fbid, comment, frstart, frend, fpres, fabs, foverlap, rrstart, rrend, rpres, rabs, roverlap
        my $insert=TEInsert->new($self->{chr},$self->{inspos},$self->{teid},$self->{order},$self->{fbid},$self->{comment},$self->{frstart},$self->{frend},$fpres,$fabs,
                                $self->{foverlap},$self->{rrstart},$self->{rrend},$rpres,$rabs,$self->{roverlap});

        return $insert;
    }
    
    sub subclone_notoverlapping
    {
        my $self=shift;
        # chr, inspos, sitesupport, teid, popfreq, order, fbid, comment
        # frstart, frend, fpopfreq, fcov, fpres, fabs, foverlap
        # rrstart, rrend, rpopfreq, rcov, rpres, rabs, roverlap
        
        # return undef if both the forward and the reverse insertion are overlapping
        return undef if($self->{foverlap} and $self->{roverlap});
        
        if($self->{foverlap})
        {
            # if the TE insert has only forward site, but it is an overlapping one -> return undef
            return undef unless $self->{rrstart};
            
            die "state not allowed"if ($self->{roverlap});
            return TEInsert->new($self->{chr},$self->{inspos},$self->{teid},$self->{order},$self->{fbid},$self->{comment},undef,undef,undef,undef,
                                undef,$self->{rrstart},$self->{rrend},$self->{rpres},$self->{rabs},$self->{roverlap});
        }
        elsif($self->{roverlap})
        {
            # if the TE insert has only one reverse site and it is an overlapping one -> return undef
            return undef unless ($self->{frstart});
            
            die "state not allwed" if ($self->{foverlap});
            return TEInsert->new($self->{chr},$self->{inspos},$self->{teid},$self->{order},$self->{fbid},$self->{comment},$self->{frstart},$self->{frend},$self->{fpres},$self->{fabs},
                                $self->{foverlap},undef,undef,undef,undef,undef);
        }
        else
        {
            return $self; # return the unmodified self if no site is overlapping
        }
    }
    
    sub subclone_mincount
    {
        my $self=shift;
        my $mincount=shift;
        
        my $fcov = $self->{fcov};
        my $rcov = $self->{rcov};
        
        if(defined($fcov) and defined($rcov))
        {
            return undef if($fcov < $mincount and $rcov < $mincount);
            return $self if($fcov >=$mincount and $rcov>=$mincount);
            if($fcov<$mincount)
            {
                return TEInsert->new($self->{chr},$self->{inspos},$self->{teid},$self->{order},$self->{fbid},$self->{comment},undef,undef,undef,undef,
                                undef,$self->{rrstart},$self->{rrend},$self->{rpres},$self->{rabs},$self->{roverlap});
            }
            elsif($rcov<$mincount)
            {
                return TEInsert->new($self->{chr},$self->{inspos},$self->{teid},$self->{order},$self->{fbid},$self->{comment},$self->{frstart},$self->{frend},$self->{fpres},$self->{fabs},
                                $self->{foverlap},undef,undef,undef,undef,undef);
            }
            else
            {
                die "impossible";
            }
            
            
            
        }
        elsif(defined($fcov))
        {
            return undef if $fcov < $mincount;
            return $self;
        }
        elsif(defined($rcov))
        {
            return undef if $rcov <$mincount;
            return $self;
        }
        else
        {
            die "impossible";
        }
        
    }
    
}
1;