#!/usr/bin/env perl

## Author: Aparna Rajpurkar
# quick perl script to handle data from simulation
# intended for use with many repeated simulations
# expects that you know what you are doing with it
# no user validation. Was intended to be quick.

# threshold for transcription
my $PERCENT_M_THRESH = .7;

# get inputs
my ($outfile, $file_list) = @ARGV;
die "usage: <outfilename> <simulation files list>\n" unless @ARGV;

# open file list
open(LI, "<", $file_list) or die "could not open $file_list:$!\n";

# get lines for each file
my @files;
while(<LI>) {
    my $line = $_;
    chomp $line;
    push(@files, $line)
}
close LI;

my %data;

# iterate over all files
foreach my $file (@files) {
    # get important run options from filename
    # which MUST BE IN THIS FORMAT
    my ($t, $d, $rt, $f, $power, $sim_num) = ($file =~ /t(\d+)_d(\d+).+rt(\d+).+f([^_]+)_*(.*)_(\d+)\.txt/);
    warn "t:$t, d:$d, rt:$rt, f:$f, power:$power, simnum:$sim_num\n";

    # open the file
    open(TMP, "<", $file) or die "could not open $file:$!\n";

    # get A and M counts for each line (timestep), calculate totals & percentages
    # calculate gap score and whether the gene is on or not
    # according to the threshold
    my $t_count = 0;
    while(<TMP>) {
        my $line = $_;
        chomp $line;
        @m_count = ($line =~ /M/g);
        @a_count = ($line =~ /A/g);

        my $m = scalar(@m_count);
        my $a = scalar(@a_count);
        my $gap = ($m - $a) / ($m + $a) ;

        my $off = 0;
        if ($m / 60 >= $PERCENT_M_THRESH) {
            $off = 1
        }

        if (!exists($data{$t}) or !exists($data{$t}{$f}) or @{$data{$t}{$f}} < $sim_num) {
            push(@{$data{$t}{$f}}, {$sim_num => {gap =>$gap, off => $off}});
        }
        else {
            $data{$t}{$f}[$t_count]{$sim_num}{gap} = $gap;
            $data{$t}{$f}[$t_count]{$sim_num}{off} = $off;
        }

        $t_count++;
    }

    close TMP;
}

# standard deviation function
sub sd {
    $mean = $_[0];
    %vals = %{$_[1]};

    my $stddev = 0;

    my $count = 0;
    foreach my $val (keys %vals) {
        $stddev += ($vals{$val}{gap} - $mean) ** 2;
        $count++;
    }

    return sqrt($stddev / $count);
}

# print data to outfile
open(OUT, ">", $outfile) or die "could not open $outfile:$!\n";
print OUT "RecruitTime\tFValue\tTimestep\tAvgGapScore\tSDGapScore\tPercOff\n";
foreach my $tot_time (keys %data) {
    foreach my $f_val (keys %{$data{$tot_time}}) {
        for (my $timestep = 0; $timestep < @{$data{$tot_time}{$f_val}}; $timestep++) {
            warn "calculating timestep $timestep\n";
            my $gap = 0;
            my $off = 0;

            foreach my $sim (keys %{$data{$tot_time}{$f_val}[$timestep]}) {
                $gap += $data{$tot_time}{$f_val}[$timestep]{$sim}{gap};
                $off += $data{$tot_time}{$f_val}[$timestep]{$sim}{off};
            }
            my $avg_gap = $gap / 100;
            my $perc_off = $off / 100;
            my $sd = sd($avg_gap, $data{$tot_time}{$f_val}[$timestep]);
            print OUT ($tot_time - 60000) . "\t$f_val\t$timestep\t$avg_gap\t$sd\t$perc_off\n";
        }
    }
}
close OUT;
