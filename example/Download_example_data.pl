#!/usr/bin/perl
#Author:ZHANG Yufeng, yfz96@connect.hku.hk
use strict;
use warnings;
use File::Basename qw(dirname basename);
use Cwd 'abs_path';

$ARGV[0] ||="example";

#if($#ARGV!=0){
#	print STDERR "perl $0 outdir\n";
#	exit 1;
#}

my $outdir=$ARGV[0];#下载目录

$outdir=abs_path($outdir);
&CheckDir("$outdir");


my @a=('saliva_data');

my %hash_path=(
	'saliva_data'=>['',]

	);

my %hash_md5=(
	'saliva_data'=>['1ea2c567daeb21b99efee44bb4734e20']
	);

print STDOUT "Download test data\n";

#download saliva data
for my $i(@a){
        my @tmp=split /\//,$hash_path{$i}[0];
        my $url=join("/",@tmp[0..$#tmp-1]);
        my $name=$tmp[-1];
        my $file_md5;#下载的文件的MD5值
        while(1){
                if(-e "$name"){
                        chomp($file_md5=`md5sum $name`);
                        $file_md5=(split /\s+/,$file_md5)[0];
                }
                if(-e "$name" && $file_md5 eq $hash_md5{$i}[0]){
                        print STDOUT "File $name has been downloaded.\n";
                        `gunzip $name`;
			last;
                }else{
                        `wget -t 0 -O $name $url`;
                }
        }
}

print STDOUT "Congratulations! All test data have been downloaded.\n";

sub CheckDir{
	my $file = shift;
	unless( -d $file ){
		if( -d dirname($file) && -w dirname($file) ){system("mkdir $file");}
		else{print STDERR "$file not exists and cannot be built\n";exit 1;}
	}
	return 1;
}






