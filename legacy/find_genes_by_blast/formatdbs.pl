#!/usr/bin/perl


#simple file to copy the contig files from a directory into the base directory

for($i = 0 ; $i <= $#ARGV ; $i++)
{
	system("formatdb -p F -i ".$ARGV[$i]);	
}