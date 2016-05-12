#!usr/bin/perl
#use strict;

my $path = $ARGV[0];                ## Argument: /path/to/file.xmfa ## MAUVE FORMAT ##

my $file;
if ($path !~ /\//) {
    $file = $path;
    $path = "\.";
}else{
    my @path = split /\//, $path;
    $file = pop @path;
    $path = join "\/", @path;
}

print "\nxmfa2fas\.pl\tLSB May 2016\n\n";
print "Program called as\: $ARGV[0]\n\n";
print "Reordering as\: default \(order in XMFA file\)\n";

my @names;                          ## Get names of genomes ##
my $numb = 0;                       ## Get number of genomes ##
my %posits;
my %gennames;
my %posct;
my $check = 0;

my ($name, $seq, $num, $lenfile, $sig);
my %rev;
open XMFA, $ARGV[0];
while(<XMFA>){
    my $l = $_;
    $lenfile .= $l;
    if ($l =~ /^#Sequence\d+File[\t|\s+](.*)/){
        my $n = $1;
        chomp $n;
        if ($n =~ /.*\\(.*)\r/){
            $n = $1;
        }
        $n =~ s/\s/\_/g;
        $n =~ s/\\/\_/;
        $n =~ s/\r//;
        $n =~ s/\.\.\///g;
        $n =~ s/\..*//;
        $numb++;
        $posits{$numb} = $n;
        $gennames{$numb} = $n;
        $posct{$numb} = 0;
    }elsif ($l =~ /^> (.*)\:(\d+\-\d+) ([\+|\-]) .*\n/ && $check == 0){
        $check = 1;
        $name = $2;
        $num = $1;
        $sig = $3;
    }elsif ($l =~ /^> (.*)\:(\d+\-\d+) ([\+|\-]) .*\n/ && $check == 1){
        if ($seq){
            $seq =~ s/\=//g;
            if ($posct{$num} == 0){
                $posits{$num} = {$name => $seq};                    # Create a Hash of Hashes!!! (HoH) #
            }else{
                $posits{$num}{$name} = $seq;                        # Add information to a HoH #
            }
            $posct{$num}++;                                         # Control of the number of breaks #
        }
        undef $name, undef $seq, undef $sig;
        $name = $2;
        $num = $1;
        $sig = $3;
    }elsif ($l !~ /^> ((.*)\:\d+\-\d+ [\+|\-] .*)\n/ && $check == 1) {
        if ($name || $name == 0){
            chomp, $_ =~ s/\r//;
            $seq .= $_;
        }
    }
}
close XMFA;
$seq =~ s/\=//g;                                                    ## Get last sequence ##
if ($posct{$numb} == 0){
    $posits{$num} = {$name => $seq};                    
}else{
    $posits{$num}{$name} = $seq;
}
$posct{$num}++;
undef $name, undef $seq, undef $sig;

my @lenfile = split /\=/, $lenfile;                                 ## Split the whole file by segments ##
pop @lenfile;
my %bone;
my $seg = 0;
my @lenfile2 = @lenfile;
do{
    $seg++;
    my $w = $seg;
    my $x = shift @lenfile2;
    my @x = split /\>/, $x;
    shift @x;
    foreach (@x){
        if (/(\d+\:\d+\-\d+)/){
            $bone{$w} .= "\|$1\|";                                    ## Get num of genome + ini position ##
        }
    }
}until scalar @lenfile2 == 0;

my $c = scalar keys %bone;

print "Nb_segments\: $c\n\n";

my %seqs;

## Segment order as XMFA ##
for (my $j=1; $j<=$seg; $j++){
    my $conts = $bone{$j};
    my @conts = split /\|/, $conts;
    @conts = grep { !/^$/ } @conts;
    my ($ini, $gen);
    my $ch = "";
    my $len;
    foreach(@conts){
        if (/^(\d+)\:(\d+\-\d+)/){
            $gen = $1;
            $ini = $2;
        }
        $ch .= "\_$gen\_";
        foreach (keys %{$posits{$gen}}){
            if ($_ == $ini){
                if(!$seqs{$gen}){
                    $seqs{$gen} = $posits{$gen}->{$ini};
                    $len = $posits{$gen}->{$ini};
                }else{
                    $seqs{$gen} .= $posits{$gen}->{$ini};
                    $len = $posits{$gen}->{$ini};
                }
                last;
            }
        }
    }
        
    my @left;
    for (my $z = 1; $z<=$numb; $z++){                                   ## Get absent genomes in current segment and add '-' ##
        if ($ch !~ /\_$z\_/){
            push @left, $z;
        }
    }
    my $ln = length $len;
    foreach (@left){
        my $lf = $_;
        for (my $x = 1; $x <= $ln; $x++){   
        $seqs{$lf} .= '-';
        }
    }
        
    # Check if different lengths
    #print "$j ";
    #for (my $z = 1; $z <= scalar keys %seqs; $z++){
    #    print length $seqs{$z};
    #    print " ";
    #}
    #print "\n";
}

## Print output ##

$file =~ s/\..*//;
open OUT, ">$path\/$file\_default_order\.fas";
foreach (keys %gennames){
    my $k = $gennames{$_};
    $k =~ s/\s//g;
    print OUT "\>$k\n$seqs{$_}\n";
}
close OUT;

if ($path eq "") {
    print "Output file saved as $file\_default_order\.fas\n\n";
}else{
    print "Output file saved as $path\/$file\_default_order\.fas\n\n";
}


exit;
