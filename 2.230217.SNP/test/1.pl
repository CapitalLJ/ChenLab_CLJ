#!usr/bin/perl

@N=("A\tT","A\tG","A\tC","T\tA","T\tG","T\tC","G\tA","G\tT","G\tC","C\tA","C\tT","C\tG");
@N1=("A->T","A->G","A->C","T->A","T->G","T->C","G->A","G->T","G->C","C->A","C->T","C->G");

# @NΪ��Ҫ��ȡ�Ĺؼ����б�
while(<>){
    push @T , $_ ;
}

for $a (@N){
    for(@T){
        chomp($_);
        if(/$a/i){
            $number{$a}+=1;
        }
    }
}
for $x (0,1,2,3,4,5,6,7,8,9,10,11,12){
    print "@N1[$x]\t$number{@N[$x]}\n";
    $total=$total+$number{@N[$x]};
}
print $total