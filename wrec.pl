#!/usr/bin/perl

$in = shift @ARGV;
$thresh = shift @ARGV;
$window_size = shift @ARGV;
$max_break = shift @ARGV;

open INPUT, $in or die "Couldn't open $in, dude.\n";
while (<INPUT>){chomp; push @positions, $_;}

for($pos = 0; $pos < scalar(@positions); $pos++){
	#print STDERR "$pos\n";
	$sum = 0;
	for($x = $pos; $x < ($pos + $window_size); $x++){
		if($positions[$x] != 2){
			$sum += $positions[$x];
		}
	}
	$ratio = $sum / $window_size;
	if($ratio >= $thresh){
	#	print STDERR "$pos exceeds threshold\n";
		if(!$in_region){
			$in_region=1;
			if($region_break){
				$break_len = 0;
				$region_break = 0;
			}
			elsif($pos < $window_size){ #left-side genome boundary condition
				push @starts, $pos;
			}else{
				for($i = ($pos - $window_size + 1); $i < $pos; $i++){
					$sum = 0;
					for($x = $i; $x < ($i + $window_size); $x++){	
						if($positions[$x] != 2){
							$sum += $positions[$x];
						}
					}
					$ratio_diff{$i} = $ratio - ($sum / $window_size);
				#	print "$i : $ratio_diff{$i}\n";
				}
				$max = 0;
				for($i = ($pos - $window_size + 1); $i < $pos; $i++){
					if($ratio_diff{$i} >= $max){
						$max = $ratio_diff{$i};
						$last_max = $i;
					}
				}
				push @starts, $last_max + $window_size;
			}
		}
		elsif($in_region){
			$region_break = 0;
			$break_len = 0;
		}
	}else{
		if($in_region){
			if($break_len > $max_break){
				$region_break = 0;
				for($i = ($pos - $max_break); $i < ($pos - $max_break + $window_size); $i++){
					$sum = 0;
					for($x = $i; $x < ($i + $window_size); $x++){
						if($positions[$x] != 2){
							
							$sum += $positions[$x];
						}
					}
					$ratio_diff{$i} = $thresh - ($sum / $window_size);
				#print "$i : $ratio_diff{$i}\n";
				}
				$max = 0;
				for($i = ($pos - $max_break); $i < ($pos - $max_break + $window_size); $i++){
					if($ratio_diff{$i} > $max){
				#		print "$i : $ratio_diff{$i}\n";
						$max = $ratio_diff{$i};
						$first_max = $i;
					}
				}

				push @stops, $first_max;
				$break_len=0;
				$in_region = 0;
			}else{
				$region_break = 1;
				$break_len++;
			}
		}
	}
}

for($i = 0; $i < scalar(@starts); $i++){
	print $starts[$i],"\t",$stops[$i],"\n";
}

#$start = unshift @starts;
#$stop = unshift @stops;


#for($i = 0; $i < scalar(@positions); $i++){
#	if($i == $stop){
#		$start = unshift @starts;
#		$stop = unshift @stops;
#	}
#	if($start <= $i && $i <= $stop){
#		print 1, "\n";
#	}else{
#		print 0, "\n";
#	}
#}
