#!/usr/bin/perl

use strict;
use warnings;

# 获取命令行参数
my $list_file = shift @ARGV;
my $filter_file = shift @ARGV;
my $output_file = shift @ARGV;

# 打开文件
open(my $list_fh, "<", $list_file) or die "Cannot open $list_file: $!";
open(my $filter_fh, "<", $filter_file) or die "Cannot open $filter_file: $!";
open(my $output_fh, ">", $output_file) or die "Cannot open $output_file: $!";

# 读取第一个文件，存储需要的信息
my %info;
while (my $line = <$list_fh>) {
    chomp $line;
    my @fields = split /\t/, $line;
    $info{$fields[1]} = [$fields[0], $fields[3], $fields[4]];
}

# 处理第二个文件，筛选并输出所需信息
while (my $line = <$filter_fh>) {
    chomp $line;
    my @fields = split /\t/, $line;

    # 如果第2列不等于第一个文件的第2列，则跳过
    next unless ($fields[1] eq $info{$fields[1]}[0]);

    # 如果第3列不在第一个文件的第4列和第5列之间，则跳过
    next unless ($fields[2] > $info{$fields[1]}[1] && $fields[2] < $info{$fields[1]}[2]);

    # 输出符合条件的信息
    print $output_fh join("\t", $info{$fields[1]}[0], @fields[0, 1, 2], $info{$fields[1]}[1], $info{$fields[1]}[2]), "\n";
}

# 关闭文件句柄
close $list_fh;
close $filter_fh;
close $output_fh;